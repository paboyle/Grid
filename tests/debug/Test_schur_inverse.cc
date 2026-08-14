/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: Test_schur_inverse.cc

    Copyright (C) 2026

Author: Peter Boyle <pboyle@bnl.gov>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

    See the full license in the file "LICENSE" in the top level distribution directory
*************************************************************************************/
/*  END LEGAL */

//
// Staged regression gate for RecursiveSchurInverse (distributed dense
// inversion by recursive Schur complement) -- the laptop-side certificate
// chain of schur_recursive_inverse_plan.txt section 4B.5.  Runs on a
// CPU-only build (Eigen BLAS backends) under mpirun:
//
//   mpirun -n 1 ./Test_schur_inverse --grid 8.8.8.8  --mpi 1.1.1.1
//   mpirun -n 2 ./Test_schur_inverse --grid 8.8.8.8  --mpi 1.1.1.2
//   mpirun -n 3 ./Test_schur_inverse --grid 8.8.8.12 --mpi 1.1.1.3
//   mpirun -n 4 ./Test_schur_inverse --grid 8.8.8.8  --mpi 1.1.1.4
//
// (n=3 exercises uneven row splits throughout.)  The lattice exists only
// to furnish the communicator; no field is ever constructed.
//
// PRECISION: the inversion runs ENTIRELY in fp64 (decision 2026-08-14,
// superseding the fp32-merge design); certificates are eps64-scaled.
// The single terminal fp32 rounding belongs to the caller (tested at the
// glue level, Test_schur_dense_coarse).
//
// Stages present (cumulative -- earlier tests are never removed):
//   T1a : ownership tables -- CheckRowStart on synthetic uneven partitions,
//         MakeRowStart allgather vs closed form on the live communicator.
//   T1b : STORAGE-CONVENTION PIN -- column-major + ld + window-offset
//         semantics fixed once via identity multiplies through the
//         explicit-ld gemmBatched, on INTEGER-VALUED data so all three
//         cases below are EXACT (values well within the mantissa):
//           (1) alpha=1,beta=0 read from an input column window
//           (2) alpha=-1,beta=1 accumulate (the S-formation case)
//           (3) write INTO an output column window, neighbours untouched
//         No later failure can be a transposition/convention ambiguity.
//   T2  : GatherGemm vs naive fp64 oracle (owner sub-ranges, alpha-beta
//         cases, tiny+huge panels, half-participation call shape).
//   T3  : LeafInvert in-place residual certificate.
//   T4  : full recursive Invert vs Eigen fp64 oracle, growth-scaled
//         certification, adversarial near-singular-A11 family with
//         telemetry-spike assertion.
//
// Hard asserts throughout; thresholds pre-registered in the plan.
//
#include <Grid/Grid.h>
#include <Grid/Grid_Eigen_Dense.h>
#include <Grid/algorithms/multigrid/RecursiveSchurInverse.h>

using namespace std;
using namespace Grid;

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  GridCartesian Comm(GridDefaultLatt(),
                     GridDefaultSimd(Nd,vComplex::Nsimd()),
                     GridDefaultMpi());
  GridBase *grid = &Comm;

  ////////////////////////////////////////////////////////////////
  // T1a : ownership tables
  ////////////////////////////////////////////////////////////////
  {
    // Synthetic partitions of N=97 (prime: every P>1 is uneven)
    const int64_t N = 97;
    for(int P=1; P<=4; P++)
    {
      std::vector<int64_t> table(P+1);
      table[0] = 0;
      for(int r=0; r<P; r++)
      {
        int64_t nr = N/P + ( (r < (int)(N%P)) ? 1 : 0 );
        table[r+1] = table[r] + nr;
      }
      RecursiveSchurInverse::CheckRowStart(table, N);
    }

    // Live allgather: deliberately uneven local counts, closed-form oracle
    int P  = grid->ProcessorCount();
    int me = grid->ThisRank();

    int64_t myNrows = 3 + me;
    std::vector<int64_t> table = RecursiveSchurInverse::MakeRowStart(grid, myNrows);

    std::vector<int64_t> expect(P+1);
    expect[0] = 0;
    for(int r=0; r<P; r++)
    {
      expect[r+1] = expect[r] + (3 + r);
    }
    GRID_ASSERT( (int)table.size() == P+1 );
    for(int r=0; r<=P; r++)
    {
      GRID_ASSERT( table[r] == expect[r] );
    }

    // Constructor smoke: derived ownership matches
    RecursiveSchurInverse RSI(grid, table[P], table, 1024*1024);
    GRID_ASSERT( RSI.P       == P );
    GRID_ASSERT( RSI.me      == me );
    GRID_ASSERT( RSI.myRow0  == expect[me] );
    GRID_ASSERT( RSI.myNrows == myNrows );

    std::cout << GridLogMessage
      << "T1a ownership tables (synthetic P=1..4, live allgather, ctor)   PASS"
      << std::endl;
  }

  ////////////////////////////////////////////////////////////////
  // T1b : storage-convention pin (every rank, local, exact)
  ////////////////////////////////////////////////////////////////
  {
    const int64_t rows = 5;
    const int64_t cols = 13;
    const int64_t col0 = 6;   // input window start
    const int64_t w    = 4;   // window width

    // f(i,j): integer-valued, unique per element
    auto f = [](int64_t i, int64_t j) -> ComplexD
    {
      return ComplexD( (RealD)(1 + i + 10*j), (RealD)(i - j) );
    };

    BlockRows A;
    A.Resize(rows, cols);
    {
      std::vector<ComplexD> Ahost((uint64_t)rows*cols);
      for(int64_t j=0; j<cols; j++)
      {
        for(int64_t i=0; i<rows; i++)
        {
          Ahost[(uint64_t)(i + j*rows)] = f(i,j);
        }
      }
      acceleratorCopyToDevice(&Ahost[0], &A.data[0], (uint64_t)rows*cols*sizeof(ComplexD));
    }

    // Identity I_w, column major
    deviceVector<ComplexD> Idev((uint64_t)w*w);
    {
      std::vector<ComplexD> Ihost((uint64_t)w*w, ComplexD(0.0,0.0));
      for(int64_t d=0; d<w; d++)
      {
        Ihost[(uint64_t)(d + d*w)] = ComplexD(1.0,0.0);
      }
      acceleratorCopyToDevice(&Ihost[0], &Idev[0], (uint64_t)w*w*sizeof(ComplexD));
    }

    GridBLAS BLAS;
    ComplexD one    ( 1.0,0.0);
    ComplexD minus  (-1.0,0.0);
    ComplexD zero   ( 0.0,0.0);

    deviceVector<ComplexD*> Ap(1);
    deviceVector<ComplexD*> Bp(1);
    deviceVector<ComplexD*> Cp(1);
    std::vector<ComplexD*>  ptr_h(1);

    auto setptr = [&](deviceVector<ComplexD*> &d, ComplexD *p)
    {
      ptr_h[0] = p;
      acceleratorCopyToDevice(&ptr_h[0], &d[0], sizeof(ComplexD*));
    };

    ////////////////////////////////////////////////////////////
    // Case 1: C = A(:, col0:col0+w) . I_w      (alpha=1, beta=0)
    ////////////////////////////////////////////////////////////
    {
      deviceVector<ComplexD> Cdev((uint64_t)rows*w);
      setptr(Ap, A.ColumnWindow(col0));
      setptr(Bp, &Idev[0]);
      setptr(Cp, &Cdev[0]);

      BLAS.gemmBatched(GridBLAS_OP_N, GridBLAS_OP_N,
                       (int)rows, (int)w, (int)w,
                       one,  Ap, (int)A.ld,
                             Bp, (int)w,
                       zero, Cp, (int)rows);
      BLAS.synchronise();

      std::vector<ComplexD> Chost((uint64_t)rows*w);
      acceleratorCopyFromDevice(&Cdev[0], &Chost[0], (uint64_t)rows*w*sizeof(ComplexD));
      for(int64_t j=0; j<w; j++)
      {
        for(int64_t i=0; i<rows; i++)
        {
          GRID_ASSERT( Chost[(uint64_t)(i + j*rows)] == f(i, col0+j) );
        }
      }
    }

    ////////////////////////////////////////////////////////////
    // Case 2: C = C0 - A(:, col0:col0+w) . I_w (alpha=-1, beta=1)
    // -- the S-formation accumulate; exact on integer data
    ////////////////////////////////////////////////////////////
    {
      auto g = [](int64_t i, int64_t j) -> ComplexD
      {
        return ComplexD( (RealD)(100 + i + j), (RealD)7 );
      };
      deviceVector<ComplexD> Cdev((uint64_t)rows*w);
      {
        std::vector<ComplexD> Chost((uint64_t)rows*w);
        for(int64_t j=0; j<w; j++)
        {
          for(int64_t i=0; i<rows; i++)
          {
            Chost[(uint64_t)(i + j*rows)] = g(i,j);
          }
        }
        acceleratorCopyToDevice(&Chost[0], &Cdev[0], (uint64_t)rows*w*sizeof(ComplexD));
      }
      setptr(Ap, A.ColumnWindow(col0));
      setptr(Bp, &Idev[0]);
      setptr(Cp, &Cdev[0]);

      BLAS.gemmBatched(GridBLAS_OP_N, GridBLAS_OP_N,
                       (int)rows, (int)w, (int)w,
                       minus, Ap, (int)A.ld,
                              Bp, (int)w,
                       one,   Cp, (int)rows);
      BLAS.synchronise();

      std::vector<ComplexD> Chost((uint64_t)rows*w);
      acceleratorCopyFromDevice(&Cdev[0], &Chost[0], (uint64_t)rows*w*sizeof(ComplexD));
      for(int64_t j=0; j<w; j++)
      {
        for(int64_t i=0; i<rows; i++)
        {
          ComplexD expect = g(i,j) - f(i, col0+j);
          GRID_ASSERT( Chost[(uint64_t)(i + j*rows)] == expect );
        }
      }
    }

    ////////////////////////////////////////////////////////////
    // Case 3: write INTO a column window of a wider C;
    // columns outside the window must be untouched
    ////////////////////////////////////////////////////////////
    {
      const int64_t ccols = 6;
      const int64_t cw0   = 2;   // output window start
      BlockRows C;
      C.Resize(rows, ccols);
      {
        std::vector<ComplexD> Chost((uint64_t)rows*ccols, ComplexD(-999.0, 999.0));
        acceleratorCopyToDevice(&Chost[0], &C.data[0], (uint64_t)rows*ccols*sizeof(ComplexD));
      }
      setptr(Ap, A.ColumnWindow(col0));
      setptr(Bp, &Idev[0]);
      setptr(Cp, C.ColumnWindow(cw0));

      BLAS.gemmBatched(GridBLAS_OP_N, GridBLAS_OP_N,
                       (int)rows, (int)w, (int)w,
                       one,  Ap, (int)A.ld,
                             Bp, (int)w,
                       zero, Cp, (int)C.ld);
      BLAS.synchronise();

      std::vector<ComplexD> Chost((uint64_t)rows*ccols);
      acceleratorCopyFromDevice(&C.data[0], &Chost[0], (uint64_t)rows*ccols*sizeof(ComplexD));
      for(int64_t j=0; j<ccols; j++)
      {
        for(int64_t i=0; i<rows; i++)
        {
          ComplexD got = Chost[(uint64_t)(i + j*rows)];
          if ( (j >= cw0) && (j < cw0+w) )
          {
            GRID_ASSERT( got == f(i, col0 + (j-cw0)) );
          }
          else
          {
            GRID_ASSERT( got == ComplexD(-999.0, 999.0) );
          }
        }
      }
    }

    std::cout << GridLogMessage
      << "T1b storage-convention pin (window read / S-accumulate / window write, exact)   PASS"
      << std::endl;
  }

  ////////////////////////////////////////////////////////////////
  // T2 : GatherGemm vs naive double-precision oracle.
  //
  // Every rank generates the SAME full N x N random fp64 operands
  // from a fixed seed (no comms needed for the oracle), keeps only
  // its own rows in BlockRows form, and after each GatherGemm call
  // checks its output window element-by-element against a plain
  // triple-loop ComplexD accumulation over the same entries.
  //
  // Sweep: N in {8, 96, 97}; owner ranges full/upper-half/single;
  // (alpha,beta) in {(1,0), (-1,1)}; panelBytes tiny (ragged
  // many-chunk gathers) and huge (single panel).  Sentinel columns
  // outside the output window must be untouched.  Finally, a
  // HALF-PARTICIPATION case rehearses the recursion call pattern:
  // lower ranks own B but pass EMPTY A/C (collectives only).
  ////////////////////////////////////////////////////////////////
  {
    int P  = grid->ProcessorCount();
    int me = grid->ThisRank();

    std::mt19937 rng(777);
    std::uniform_real_distribution<double> dist(-1.0,1.0);

    const int64_t nout = 5;   // output width
    const int64_t colB = 3;   // B window offset
    const int64_t colC = 2;   // C window offset

    for(int64_t N : {8L, 96L, 97L})
    {
      // Ownership: uneven for any P not dividing N
      std::vector<int64_t> table(P+1);
      table[0] = 0;
      for(int r=0; r<P; r++)
      {
        int64_t nr = N/P + ( (r < (int)(N%P)) ? 1 : 0 );
        table[r+1] = table[r] + nr;
      }
      int64_t r0   = table[me];
      int64_t myNr = table[me+1] - table[me];

      // Identical full operands on every rank
      std::vector<ComplexD> Aglob((uint64_t)N*N);
      std::vector<ComplexD> Bglob((uint64_t)N*N);
      for(uint64_t i=0; i<(uint64_t)N*N; i++) Aglob[i] = ComplexD(dist(rng),dist(rng));
      for(uint64_t i=0; i<(uint64_t)N*N; i++) Bglob[i] = ComplexD(dist(rng),dist(rng));

      // My rows of a full-matrix operand as a BlockRows
      auto fillRows = [&](BlockRows &X, std::vector<ComplexD> &glob,
                          int64_t row0, int64_t nr)
      {
        X.Resize(nr, N);
        if ( nr == 0 ) return;
        std::vector<ComplexD> h((uint64_t)nr*N);
        for(int64_t j=0; j<N; j++)
        {
          for(int64_t i=0; i<nr; i++)
          {
            h[(uint64_t)(i + j*nr)] = glob[(uint64_t)((row0+i) + j*N)];
          }
        }
        acceleratorCopyToDevice(&h[0], &X.data[0], (uint64_t)nr*N*sizeof(ComplexD));
      };

      // Owner-range cases: full span, upper half, single interior rank
      std::vector<std::pair<int,int> > ranges;
      ranges.push_back(std::make_pair(0, P));
      if ( P > 1 ) ranges.push_back(std::make_pair(P/2, P));
      if ( P > 1 ) ranges.push_back(std::make_pair(1, 2));

      for(auto range : ranges)
      {
        int rB0 = range.first;
        int rB1 = range.second;
        int64_t ka0 = table[rB0];       // A-column window start = B row span
        int64_t k   = table[rB1] - table[rB0];

        for(int acase=0; acase<2; acase++)
        {
          ComplexD alpha = ( acase==0 ) ? ComplexD( 1.0,0.0) : ComplexD(-1.0,0.0);
          ComplexD beta  = ( acase==0 ) ? ComplexD( 0.0,0.0) : ComplexD( 1.0,0.0);

          for(int64_t panelBytes : {64L, 1L<<30})
          {
            RecursiveSchurInverse RSI(grid, N, table, panelBytes);

            BlockRows A;
            BlockRows B;
            BlockRows C;
            fillRows(A, Aglob, r0, myNr);
            fillRows(B, Bglob, r0, myNr);

            // Output: sentinel-filled, window at colC
            const ComplexD sentinel(-999.0, 999.0);
            const int64_t  ccols = colC + nout + 2;
            C.Resize(myNr, ccols);
            std::vector<ComplexD> C0((uint64_t)myNr*ccols, sentinel);
            if ( acase == 1 )
            {
              // beta=1 needs defined window content: g(i,j), integer-valued
              for(int64_t j=0; j<nout; j++)
              {
                for(int64_t i=0; i<myNr; i++)
                {
                  C0[(uint64_t)(i + (colC+j)*myNr)] = ComplexD((RealD)(50+i+j), (RealD)-3);
                }
              }
            }
            if ( myNr > 0 )
            {
              acceleratorCopyToDevice(&C0[0], &C.data[0], (uint64_t)myNr*ccols*sizeof(ComplexD));
            }

            RSI.GatherGemm(alpha, A, ka0, k,
                           rB0, rB1,
                           B, colB, nout,
                           beta, C, colC);

            std::vector<ComplexD> Chost((uint64_t)myNr*ccols);
            if ( myNr > 0 )
            {
              acceleratorCopyFromDevice(&C.data[0], &Chost[0], (uint64_t)myNr*ccols*sizeof(ComplexD));
            }

            double tol = 1.0e-14 * (double)k;
            for(int64_t j=0; j<ccols; j++)
            {
              for(int64_t i=0; i<myNr; i++)
              {
                ComplexD got = Chost[(uint64_t)(i + j*myNr)];
                if ( (j >= colC) && (j < colC+nout) )
                {
                  int64_t jj = j - colC;
                  ComplexD acc(0.0,0.0);
                  if ( acase == 1 )
                  {
                    acc = C0[(uint64_t)(i + j*myNr)];
                  }
                  for(int64_t t=0; t<k; t++)
                  {
                    acc += alpha
                         * Aglob[(uint64_t)((r0+i)  + (ka0+t)*N)]
                         * Bglob[(uint64_t)((ka0+t) + (colB+jj)*N)];
                  }
                  GRID_ASSERT( abs(got - acc) < tol );
                }
                else
                {
                  GRID_ASSERT( got == sentinel );
                }
              }
            }
          }
        }
      }

      ////////////////////////////////////////////////////////////
      // Half-participation: owners = [0,ph) hold B; participants
      // = [ph,P) hold A/C; owners pass EMPTY A/C and column
      // offset 0 (collectives only) -- the recursion call shape.
      ////////////////////////////////////////////////////////////
      if ( P > 1 )
      {
        int ph = ( P+1 ) / 2;
        int64_t ka0 = table[0];
        int64_t k   = table[ph] - table[0];
        int participant = ( me >= ph );

        RecursiveSchurInverse RSI(grid, N, table, 64);

        BlockRows A;
        BlockRows B;
        BlockRows C;
        fillRows(B, Bglob, r0, myNr);
        if ( participant )
        {
          fillRows(A, Aglob, r0, myNr);
          C.Resize(myNr, nout);
        }

        ComplexD one (1.0,0.0);
        ComplexD zero(0.0,0.0);
        int64_t  cA = participant ? ka0 : 0;
        RSI.GatherGemm(one, A, cA, k,
                       0, ph,
                       B, colB, nout,   // owners deposit from their B window
                       zero, C, 0);

        if ( participant )
        {
          std::vector<ComplexD> Chost((uint64_t)myNr*nout);
          acceleratorCopyFromDevice(&C.data[0], &Chost[0], (uint64_t)myNr*nout*sizeof(ComplexD));
          double tol = 1.0e-14 * (double)k;
          for(int64_t j=0; j<nout; j++)
          {
            for(int64_t i=0; i<myNr; i++)
            {
              ComplexD acc(0.0,0.0);
              for(int64_t t=0; t<k; t++)
              {
                acc += Aglob[(uint64_t)((r0+i)  + (ka0+t)*N)]
                     * Bglob[(uint64_t)((ka0+t) + (colB+j)*N)];
              }
              GRID_ASSERT( abs(Chost[(uint64_t)(i + j*myNr)] - acc) < tol );
            }
          }
        }
      }
    }

    std::cout << GridLogMessage
      << "T2  GatherGemm vs oracle (N=8/96/97, 3 owner ranges, 2 alpha-beta, tiny+huge panels, half-participation)   PASS"
      << std::endl;
  }

  ////////////////////////////////////////////////////////////////
  // T3 : LeafInvert -- in-place fp64 inversion of the contiguous
  // leaf window.  Purely local, every rank runs its own
  // uneven-size leaf; residual certificate in ComplexD.
  ////////////////////////////////////////////////////////////////
  {
    int P  = grid->ProcessorCount();
    int me = grid->ThisRank();

    int64_t  w   = 17 + 3*me;
    uint64_t len = (uint64_t)w*w;

    std::vector<int64_t> table = RecursiveSchurInverse::MakeRowStart(grid, w);
    RecursiveSchurInverse RSI(grid, table[P], table, 1<<20);

    // A = w I + R : well conditioned
    std::mt19937 rng(31 + me);
    std::uniform_real_distribution<double> dist(-1.0,1.0);
    std::vector<ComplexD> Ahost(len);
    for(uint64_t i=0; i<len; i++) Ahost[i] = ComplexD(dist(rng),dist(rng));
    for(int64_t d=0; d<w; d++)  Ahost[(uint64_t)(d + d*w)] += ComplexD((RealD)w, 0.0);

    BlockRows Ar;
    Ar.Resize(w, w);
    acceleratorCopyToDevice(&Ahost[0], &Ar.data[0], len*sizeof(ComplexD));
    RSI.LeafInvert(0, w, Ar);
    std::vector<ComplexD> X(len);
    acceleratorCopyFromDevice(&Ar.data[0], &X[0], len*sizeof(ComplexD));

    double maxdev = 0.0;
    for(int64_t j=0; j<w; j++)
    {
      for(int64_t i=0; i<w; i++)
      {
        ComplexD acc(0.0,0.0);
        for(int64_t t=0; t<w; t++)
        {
          acc += Ahost[(uint64_t)(i + t*w)] * X[(uint64_t)(t + j*w)];
        }
        if ( i==j ) acc -= ComplexD(1.0,0.0);
        maxdev = std::max(maxdev, abs(acc));
      }
    }
    GRID_ASSERT( maxdev < 1.0e-13 );

    std::cout << GridLogMessage
      << "T3  LeafInvert in-place fp64 (residual " << maxdev << ")   PASS" << std::endl;
  }

  ////////////////////////////////////////////////////////////////
  // T4 : full recursive Invert vs Eigen fp64 oracle.
  //
  // Every rank builds the SAME N x N fp64 matrix from a fixed seed,
  // keeps its rows, inverts through the full SPMD recursion, then
  // the test gathers the complete inverse (zero-fill GlobalSum) and
  // checks BOTH certificates:
  //     cert1 = || A X - I ||_max          (ComplexD accumulation)
  //     cert2 = max|X - Xref| / max|Xref|  (Xref = Eigen fp64 inverse)
  //
  // Families (eps64-scaled tolerances; the fp32-era growth data
  // rescales by eps64/eps32 ~ 1.9e-9):
  //   kappa-moderate : A = R + 3 sqrt(N) I
  //   kappa-large    : A = R + 0.3 sqrt(N) I
  //   adversarial    : leading block (rank 0's whole leaf) REPLACED by
  //                    1e-2 * (R' + 3 sqrt(b) I) inside a well-conditioned
  //                    A -- the growth spike must REGISTER in telemetry
  //                    (asserted > 10 when P > 1); at fp64 the certificate
  //                    barely notices it: that insensitivity IS the point
  //                    of the fp64 conversion.
  //
  // N=64 runs with panelBytes=128 (ragged many-chunk gathers inside
  // the recursion); larger N with 1 MB panels.
  ////////////////////////////////////////////////////////////////
  {
    int P  = grid->ProcessorCount();
    int me = grid->ThisRank();

    std::mt19937 rng(2026);
    std::uniform_real_distribution<double> dist(-1.0,1.0);

    for(int64_t N : {64L, 200L, 513L})
    {
      std::vector<int64_t> table(P+1);
      table[0] = 0;
      for(int r=0; r<P; r++)
      {
        int64_t nr = N/P + ( (r < (int)(N%P)) ? 1 : 0 );
        table[r+1] = table[r] + nr;
      }
      int64_t r0   = table[me];
      int64_t myNr = table[me+1] - table[me];

      for(int fam=0; fam<3; fam++)
      {
        const char *famname = (fam==0) ? "kappa-moderate" :
                              (fam==1) ? "kappa-large"    : "adversarial-A11";
        double shift = (fam==1) ? 0.3*std::sqrt((double)N) : 3.0*std::sqrt((double)N);
        double tol   = (fam==0) ? 1.0e-12 :
                       (fam==1) ? 1.0e-11 : 5.0e-11;

        // Identical operand on every rank (all draws rank-independent)
        std::vector<ComplexD> Aglob((uint64_t)N*N);
        for(uint64_t i=0; i<(uint64_t)N*N; i++) Aglob[i] = ComplexD(dist(rng),dist(rng));
        for(int64_t d=0; d<N; d++) Aglob[(uint64_t)(d + d*N)] += ComplexD(shift, 0.0);
        if ( fam == 2 )
        {
          // Leading block = rank 0's whole leaf, scaled down 100x but
          // internally well conditioned (shift scales as sqrt(b): a
          // FIXED shift makes A11 itself near-singular at large b).
          int64_t b = ( P > 1 ) ? table[1] : N/4;
          for(int64_t j=0; j<b; j++)
          {
            for(int64_t i=0; i<b; i++)
            {
              Aglob[(uint64_t)(i + j*N)] = ComplexD(0.01,0.0)*ComplexD(dist(rng),dist(rng));
            }
          }
          RealD bshift = (RealD)(0.03*std::sqrt((double)b));
          for(int64_t d=0; d<b; d++) Aglob[(uint64_t)(d + d*N)] += ComplexD(bshift,0.0);
        }

        // Eigen fp64 oracle.  Explicit re/im conversion at the boundary:
        // on HIP builds ComplexD is thrust::complex, which has no
        // operators against Eigen's std::complex.
        auto toStd = [](const ComplexD &z) -> std::complex<double>
        {
          return std::complex<double>(z.real(), z.imag());
        };
        Eigen::MatrixXcd eA(N,N);
        for(int64_t j=0; j<N; j++)
        {
          for(int64_t i=0; i<N; i++)
          {
            eA(i,j) = toStd(Aglob[(uint64_t)(i + j*N)]);
          }
        }
        Eigen::MatrixXcd Xref = eA.inverse();

        // Distribute, invert
        int64_t panelBytes = ( N == 64 ) ? 128 : (1<<20);
        RecursiveSchurInverse RSI(grid, N, table, panelBytes);

        BlockRows Arows;
        Arows.Resize(myNr, N);
        {
          std::vector<ComplexD> h((uint64_t)myNr*N);
          for(int64_t j=0; j<N; j++)
          {
            for(int64_t i=0; i<myNr; i++)
            {
              h[(uint64_t)(i + j*myNr)] = Aglob[(uint64_t)((r0+i) + j*N)];
            }
          }
          acceleratorCopyToDevice(&h[0], &Arows.data[0], (uint64_t)myNr*N*sizeof(ComplexD));
        }

        RSI.Invert(Arows);

        // Gather the full inverse: zero-fill + GlobalSum
        std::vector<ComplexD> Xfull((uint64_t)N*N, ComplexD(0.0,0.0));
        {
          std::vector<ComplexD> h((uint64_t)myNr*N);
          acceleratorCopyFromDevice(&Arows.data[0], &h[0], (uint64_t)myNr*N*sizeof(ComplexD));
          for(int64_t j=0; j<N; j++)
          {
            for(int64_t i=0; i<myNr; i++)
            {
              Xfull[(uint64_t)((r0+i) + j*N)] = h[(uint64_t)(i + j*myNr)];
            }
          }
        }
        grid->GlobalSumVector(&Xfull[0], (int)(N*N));

        // cert1 = ||A X - I||_max
        double cert1 = 0.0;
        for(int64_t j=0; j<N; j++)
        {
          for(int64_t i=0; i<N; i++)
          {
            ComplexD acc(0.0,0.0);
            for(int64_t t=0; t<N; t++)
            {
              acc += Aglob[(uint64_t)(i + t*N)] * Xfull[(uint64_t)(t + j*N)];
            }
            if ( i==j ) acc -= ComplexD(1.0,0.0);
            cert1 = std::max(cert1, abs(acc));
          }
        }

        // cert2 = max|X - Xref| / max|Xref|
        double maxref = 0.0;
        double maxdif = 0.0;
        for(int64_t j=0; j<N; j++)
        {
          for(int64_t i=0; i<N; i++)
          {
            maxref = std::max(maxref, std::abs(Xref(i,j)));
            maxdif = std::max(maxdif, std::abs(toStd(Xfull[(uint64_t)(i + j*N)]) - Xref(i,j)));
          }
        }
        double cert2 = maxdif / maxref;

        double maxNormB = 0.0;
        for(uint64_t i=0; i<RSI.telNormB.size(); i++)
        {
          maxNormB = std::max(maxNormB, RSI.telNormB[i]);
        }

        // GROWTH-SCALED certification, eps64 (the fp32-era model with
        // eps swapped: cert2 ~ (10-12) ||B||_F sqrt(N) eps; threshold =
        // 3x margin, floored at the family tolerance).  ||B||_F capped
        // per family so growth cannot silently excuse a logic error.
        // cert1 remains a loose absolute bound (an O(1) logic error
        // gives cert1 ~ 1e2-1e3; fp64 rounding gives ~1e-10).
        double eps64    = 2.3e-16;
        double tolModel = 30.0 * std::max(1.0, maxNormB) * std::sqrt((double)N) * eps64;
        double tolEff   = std::max(tol, tolModel);
        double capB     = (fam==0) ? 100.0 : (fam==1) ? 2000.0 : 10000.0;

        std::cout << GridLogMessage
          << "T4  N=" << N << " " << famname
          << "  ||AX-I||_max " << cert1
          << "  |X-Xref|/|Xref| " << cert2
          << "  max||B||_F " << maxNormB
          << "  tolEff " << tolEff
          << ( (cert2 < tolEff) && (cert1 < 1.0e-6) ? "   PASS" : "   FAIL" )
          << std::endl;

        GRID_ASSERT( cert2 < tolEff );
        GRID_ASSERT( cert1 < 1.0e-6 );
        GRID_ASSERT( maxNormB < capB );
        if ( (fam == 2) && (P > 1) )
        {
          GRID_ASSERT( maxNormB > 10.0 );   // the spike must REGISTER
        }
      }
    }

    std::cout << GridLogMessage
      << "T4  recursive Invert vs Eigen oracle (N=64/200/513, 3 families)   PASS"
      << std::endl;
  }

  std::cout << GridLogMessage
    << "Test_schur_inverse: ALL STAGES PASS" << std::endl;

  Grid_finalize();
}
