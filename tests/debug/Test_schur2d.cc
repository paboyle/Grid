/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./tests/debug/Test_schur2d.cc

    Copyright (C) 2026

Author: Peter Boyle <pboyle@bnl.gov>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See the full license in the file "LICENSE" in the top level distribution
    directory
*************************************************************************************/
/*  END LEGAL */

//////////////////////////////////////////////////////////////////////////////
// Regression gate for BlockCyclicSchurInverse -- stage 3 of the 2D
// distributed dense inverse.  CPU build under mpirun:
//
//   mpirun -n 1 ./Test_schur2d --grid 16.16.16.32 --mpi 1.1.1.1
//   mpirun -n 2 ./Test_schur2d --grid 16.16.16.32 --mpi 1.1.1.2
//   mpirun -n 3 ./Test_schur2d --grid 16.16.16.48 --mpi 1.1.1.3
//   mpirun -n 4 ./Test_schur2d --grid 16.16.16.32 --mpi 1.1.1.4
//
// Sweeps all process-grid factorisations of P and a battery of (N,nb)
// including ragged trailing blocks, a single-leaf matrix (nblocks==1),
// and nb=3 with many blocks.  The matrices are diagonally dominant --
// the recursion does not pivot, exactly like the 1D implementation, and
// the test respects that contract.
//
//   T1 : certificate  max|A . Ainv - I|  with the product computed by the
//        (independently validated) distributed SUMMA on an untouched copy.
//   T2 : element-wise against a host Gauss-Jordan reference inverse
//        (partial pivoting, fp64).
//   T3 : repeated inversion bitwise identical -- the P2P determinism
//        property, which a collective-reduce implementation cannot offer.
//////////////////////////////////////////////////////////////////////////////

#include <Grid/Grid.h>
#include <Grid/algorithms/multigrid/BlockCyclicSchurInverse.h>

using namespace Grid;

static int failures = 0;

// Portable |z|: ComplexD is std::complex on CPU builds and thrust::complex
// under HIP, where std::abs does not resolve (same trap RecursiveSchurInverse
// documents at FrobNorm2Local).  Member real()/imag() work on both.
static double Cabs(const ComplexD &z)
{
  double re = z.real(), im = z.imag();
  return std::sqrt(re*re + im*im);
}

static void Report(const std::string &name, bool pass, const std::string &detail="")
{
  std::cout << GridLogMessage << "  " << name << (pass ? "  PASS" : "  ** FAIL **");
  if ( detail.size() ) std::cout << "   " << detail;
  std::cout << std::endl;
  if ( !pass ) failures++;
}

static ComplexD Fill(int64_t i, int64_t j, int salt)
{
  double x = std::sin(0.7*i + 1.3*j + 0.31*salt);
  double y = std::cos(1.9*i - 0.4*j + 0.77*salt);
  return ComplexD(x,y);
}

// Diagonally dominant test matrix: no pivoting required at any depth.
static void MakeMatrix(std::vector<ComplexD> &A, int64_t N, int salt)
{
  A.resize((uint64_t)N*N);
  for(int64_t j=0;j<N;j++)
    for(int64_t i=0;i<N;i++)
      A[i+j*N] = Fill(i,j,salt) + ((i==j) ? ComplexD(3.0*N, 0.5) : ComplexD(0.0,0.0));
}

// Host reference inverse: Gauss-Jordan with partial pivoting, fp64.
static void HostInverse(std::vector<ComplexD> A, std::vector<ComplexD> &X, int64_t N)
{
  X.assign((uint64_t)N*N, ComplexD(0.0,0.0));
  for(int64_t i=0;i<N;i++) X[i+i*N] = ComplexD(1.0,0.0);
  for(int64_t c=0;c<N;c++){
    int64_t piv=c; double mx = Cabs(A[c+c*N]);
    for(int64_t r=c+1;r<N;r++)
      if ( Cabs(A[r+c*N]) > mx ){ mx=Cabs(A[r+c*N]); piv=r; }
    GRID_ASSERT( mx > 0.0 );
    if ( piv != c )
      for(int64_t j=0;j<N;j++){
        std::swap(A[c+j*N],A[piv+j*N]);
        std::swap(X[c+j*N],X[piv+j*N]);
      }
    ComplexD d = ComplexD(1.0,0.0)/A[c+c*N];
    for(int64_t j=0;j<N;j++){ A[c+j*N]*=d; X[c+j*N]*=d; }
    for(int64_t r=0;r<N;r++){
      if ( r==c ) continue;
      ComplexD f = A[r+c*N];
      if ( f == ComplexD(0.0,0.0) ) continue;
      for(int64_t j=0;j<N;j++){
        A[r+j*N] -= f*A[c+j*N];
        X[r+j*N] -= f*X[c+j*N];
      }
    }
  }
}

int main(int argc, char **argv)
{
  Grid_init(&argc, &argv);

  GridCartesian *grid = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(),
                                                       GridDefaultSimd(Nd, vComplexD::Nsimd()),
                                                       GridDefaultMpi());
  const int P = grid->ProcessorCount();

  std::vector<std::pair<int,int>> grids;
  for(int r=1;r<=P;r++) if ( P%r==0 ) grids.push_back({r,P/r});

  struct Cfg { int64_t N; int64_t nb; };
  std::vector<Cfg> cfgs = { {24,4}, {26,4}, {17,5}, {30,7}, {8,8}, {33,3}, {40,5} };

  BlockCyclicSchurInverse RSI;
  BlockCyclicSumma        SUMMA;

  std::cout << GridLogMessage << "BlockCyclicSchurInverse regression: P=" << P
            << ", " << grids.size() << " process grids, " << cfgs.size()
            << " layouts" << std::endl;

  ////////////////////////////////////////////////////////////////////////
  // T1 + T2, one pass: invert, certify with distributed SUMMA on an
  // untouched copy, and compare against the host reference inverse.
  ////////////////////////////////////////////////////////////////////////
  {
    bool okC = true, okR = true;
    double worstC = 0.0, worstR = 0.0;
    for(auto &g : grids){
      for(auto &c : cfgs){
        int64_t N = c.N;
        std::vector<ComplexD> Ag, Ref, Ainv, Cert;
        MakeMatrix(Ag, N, 12);
        HostInverse(Ag, Ref, N);

        BlockCyclicMatrix A (grid,N,c.nb,g.first,g.second);
        BlockCyclicMatrix A0(grid,N,c.nb,g.first,g.second);
        BlockCyclicMatrix Ce(grid,N,c.nb,g.first,g.second);
        A.ImportGlobal(Ag);
        A0.ImportGlobal(Ag);

        RSI.Invert(A);                       // in place: A now holds Ainv

        // certificate: Ce = A0 . Ainv, distributed
        SUMMA.Multiply(ComplexD(1.0,0.0),A0,A,ComplexD(0.0,0.0),Ce, 0,N,0,N,0,N);
        Ce.ExportGlobal(Cert);
        double dc = 0.0;
        for(int64_t j=0;j<N;j++)
          for(int64_t i=0;i<N;i++){
            ComplexD id = (i==j) ? ComplexD(1.0,0.0) : ComplexD(0.0,0.0);
            dc = std::max(dc, Cabs(Cert[i+j*N]-id));
          }
        worstC = std::max(worstC,dc);
        if ( dc > 1.0e-10 ) okC = false;

        // reference: element-wise, scaled by the largest inverse entry
        A.ExportGlobal(Ainv);
        double mxref = 0.0, dr = 0.0;
        for(uint64_t i=0;i<Ref.size();i++) mxref = std::max(mxref, Cabs(Ref[i]));
        for(uint64_t i=0;i<Ref.size();i++) dr = std::max(dr, Cabs(Ainv[i]-Ref[i]));
        dr /= mxref;
        worstR = std::max(worstR,dr);
        if ( dr > 1.0e-9 ) okR = false;
      }
    }
    Report("T1  certificate max|A.Ainv - I|, all grids x layouts", okC,
           "worst "+std::to_string(worstC));
    Report("T2  vs host Gauss-Jordan reference (relative)", okR,
           "worst "+std::to_string(worstR));
  }

  ////////////////////////////////////////////////////////////////////////
  // T3 : determinism.  Same import, two inversions, bitwise comparison.
  ////////////////////////////////////////////////////////////////////////
  {
    bool ok = true;
    for(auto &g : grids){
      int64_t N = 30, nb = 7;
      std::vector<ComplexD> Ag, X1, X2;
      MakeMatrix(Ag, N, 13);
      BlockCyclicMatrix A(grid,N,nb,g.first,g.second);
      A.ImportGlobal(Ag); RSI.Invert(A); A.ExportGlobal(X1);
      A.ImportGlobal(Ag); RSI.Invert(A); A.ExportGlobal(X2);
      for(uint64_t i=0;i<X1.size();i++)
        if ( !(X1[i]==X2[i]) ) ok = false;   // BITWISE
    }
    Report("T3  repeated inversion bitwise identical", ok);
  }

  {
    uint64_t f = failures;  grid->GlobalSum(f);
    if ( f && !failures )
      std::cout << GridLogMessage << "  ** failures on OTHER ranks: " << f << " **" << std::endl;
    failures = (int)f;
  }
  std::cout << GridLogMessage << (failures ? "Test_schur2d: FAILURES"
                                           : "Test_schur2d: ALL PASS") << std::endl;
  Grid_finalize();
  return failures ? 1 : 0;
}
