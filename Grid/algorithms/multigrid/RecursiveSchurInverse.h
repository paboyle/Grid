/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: RecursiveSchurInverse.h

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
#pragma once

#include <Grid/algorithms/blas/BatchedBlas.h>
#include <Grid/algorithms/blas/BatchedInverse.h>

NAMESPACE_BEGIN(Grid);

///////////////////////////////////////////////////////////////////////////////
// RecursiveSchurInverse: distributed dense inversion by recursive Schur
// complement over a binary rank-range tree.
//
// CONTRACT: the caller presents an N x N matrix in RANK-MAJOR row ordering,
// distributed by rows -- rank r owns global rows [rowStart[r], rowStart[r+1])
// -- and receives its rows of the INVERSE in the same layout.  This class
// knows nothing of lattices or coarse operators; it consumes a GridBase for
// world collectives, GridBLAS for GEMMs and GridBLASInverse for the leaf
// inversions (all of which have Eigen reference backends, so the whole
// algorithm unit-tests on a CPU-only laptop build under mpirun).
//
// PRECISION (decision 2026-08-14, superseding the fp32-merge design): the
// ENTIRE inversion runs in fp64 (ComplexD).  The apply-side fp32 gain is
// taken where it matters -- inside the iterative process -- by rounding the
// finished inverse ONCE when the caller stores it in the fp32 apply slab.
// Consequences: merge-growth error accumulates in eps64 and the terminal
// rounding gives representation-only ~eps32 accuracy independent of growth;
// the Newton-Schulz refinement and the fp32 escalation ladder are DELETED
// (resurrectable from git history if a future scale forces reduced-precision
// merges).  Setup cost: ~2x panel-gather bytes and ~2x transient memory,
// once per setup; fp64 GEMM runs at fp32 rate on CDNA2/PVC.
//
// EXECUTION MODEL: SPMD full-tree walk.  Every rank executes the identical
// recursion call sequence; participation in DATA is ownership-gated, and
// every collective is a world-communicator zero-fill GlobalSumVector.  No
// sub-communicators exist, so no deadlock surface exists.
//
// STORAGE CONVENTION (pinned by unit test T1b, Test_schur_inverse.cc):
// BlockRows is COLUMN-MAJOR with ld = rows, matching the BLAS world:
// element (i,j) lives at data[ i + j*ld ]; a column window [col0, col0+w)
// is the contiguous slice starting at data[ col0*ld ].
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// My rows of a distributed dense matrix: rows x cols, column major, ld = rows.
///////////////////////////////////////////////////////////////////////////////
class BlockRows
{
public:
  deviceVector<ComplexD> data;
  int64_t                rows;
  int64_t                cols;
  int64_t                ld;

  BlockRows()
  {
    rows = 0;
    cols = 0;
    ld   = 0;
  }
  void Resize(int64_t r, int64_t c)
  {
    rows = r;
    cols = c;
    ld   = r;
    data.resize((uint64_t)r*c);
  }
  ComplexD *ColumnWindow(int64_t col0)
  {
    GRID_ASSERT( col0 >= 0 );
    GRID_ASSERT( col0 <= cols );
    return &data[(uint64_t)col0*ld];
  }
};

class RecursiveSchurInverse
{
public:
  GridBase             *grid;       // world collectives only
  int64_t               N;          // global matrix dimension
  int                   P;          // ranks
  int                   me;         // this rank
  std::vector<int64_t>  rowStart;   // P+1 entries: rank-major row ownership
  int64_t               myRow0;
  int64_t               myNrows;
  int64_t               panelBytes; // gather panel budget (DENSE_PANEL_BYTES)

  GridBLAS              BLAS;
  GridBLASInverse       INV;

  // Growth telemetry (diagnostic, not load-bearing at fp64): one entry per
  // merge node, walk order
  std::vector<double>   telNormB;      // ||B||_F = ||A11inv A12||_F
  std::vector<double>   telSratio;     // ||S||_F / ||A22||_F
  double                telLeafMaxInv; // max |(leaf inverse)_ij| over leaves

  // Phase timers/counters (this rank), accumulated across the whole walk:
  // where does the setup wall actually go?  Reported by ReportTelemetry.
  double                tStage;        // owner B-window device->host
  double                tMemset;       // panel zero-fill (host)
  double                tDeposit;      // owner rows -> panel (host memcpy)
  double                tAllreduce;    // GlobalSumVector on panels
  double                tH2D;          // panel host->device
  double                tGemm;         // strided gemm + synchronise
  double                tLeaf;         // leaf inversions
  double                tARmin;        // fastest single panel collective
  double                tARmax;        // slowest single panel collective
  uint64_t              bytesAllreduce;
  uint64_t              nAllreduce;    // panel collectives
  uint64_t              nGatherGemm;   // GatherGemm calls

  // PERSISTENT comms/staging buffers, grow-only across the whole walk.
  // Fresh per-call host allocations defeat the MPI registration cache
  // (every page unpinned every call => re-registration or bounce-buffer
  // copies inside MPI on ~GB panels) -- the same reason the halo
  // exchange uses allocate-once buffers.  Measured motivation: 0.48 GB/s
  // effective allreduce payload at N=138k with per-call vectors.
  std::vector<ComplexD>  panelBuf;
  deviceVector<ComplexD> dPanelBuf;
  std::vector<ComplexD>  stageBuf;

  ///////////////////////////////////////////////////////////////////////////
  // Ownership-table validation: a proper partition of [0,N).
  // Static and communicator-free so synthetic tables unit-test directly.
  ///////////////////////////////////////////////////////////////////////////
  static void CheckRowStart(const std::vector<int64_t> &table, int64_t N)
  {
    int P = (int)table.size() - 1;
    GRID_ASSERT( P >= 1 );
    GRID_ASSERT( table[0] == 0 );
    GRID_ASSERT( table[P] == N );
    for(int r=0; r<P; r++)
    {
      GRID_ASSERT( table[r+1] >= table[r] );  // zero-row ranks permitted
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // Build the ownership table from each rank's local row count: zero-fill
  // allgather (the standing comms idiom) then prefix sum.  Every rank
  // returns the identical table.
  ///////////////////////////////////////////////////////////////////////////
  static std::vector<int64_t> MakeRowStart(GridBase *g, int64_t myNrows)
  {
    int P  = g->ProcessorCount();
    int me = g->ThisRank();

    std::vector<uint64_t> counts(P, 0);
    counts[me] = (uint64_t)myNrows;
    g->GlobalSumVector(&counts[0], P);

    std::vector<int64_t> table(P+1);
    table[0] = 0;
    for(int r=0; r<P; r++)
    {
      table[r+1] = table[r] + (int64_t)counts[r];
    }
    CheckRowStart(table, table[P]);
    return table;
  }

  RecursiveSchurInverse(GridBase *g,
                        int64_t N_,
                        std::vector<int64_t> &rowStart_,
                        int64_t panelBytes_)
  {
    grid       = g;
    N          = N_;
    P          = g->ProcessorCount();
    me         = g->ThisRank();
    rowStart   = rowStart_;
    panelBytes = panelBytes_;

    GRID_ASSERT( (int)rowStart.size() == P+1 );
    CheckRowStart(rowStart, N);

    myRow0  = rowStart[me];
    myNrows = rowStart[me+1] - rowStart[me];

    telLeafMaxInv = 0.0;
  }

  ///////////////////////////////////////////////////////////////////////////
  // THE communication primitive (plan 3.4 / 4B.3).
  //
  //   C(:, colC : colC+widthB) <-   beta  * C(:, colC : colC+widthB)
  //                               + alpha * A(:, colA : colA+widthA) * Bsub
  //
  // Bsub is the widthA x widthB sub-block of a row-distributed operand
  // owned by ranks [rB0, rB1): owner r contributes its rows of
  // B(:, colB : colB+widthB) at sub-block row offset
  // rowStart[r] - rowStart[rB0].  The sub-block is gathered in panelBytes
  // row-chunks by host zero-fill + world GlobalSumVector.
  //
  // SPMD rules: EVERY rank calls (the collectives are world-wide);
  // non-owners of B add zeros; ranks with A.rows == 0 skip all local
  // compute but still make every collective call.  Column offsets are
  // LOCAL buffer offsets -- non-participants pass 0.
  //
  // Owners stage their whole B window device->host ONCE (ld == rows makes
  // the window contiguous); per-chunk deposits are host memcpy runs.
  ///////////////////////////////////////////////////////////////////////////
  void GatherGemm(ComplexD alpha,
                  BlockRows &A, int64_t colA, int64_t widthA,
                  int rB0, int rB1,
                  BlockRows &B, int64_t colB, int64_t widthB,
                  ComplexD beta,
                  BlockRows &C, int64_t colC)
  {
    GRID_ASSERT( rB0 >= 0 );
    GRID_ASSERT( rB1 >  rB0 );
    GRID_ASSERT( rB1 <= P );

    int64_t k = rowStart[rB1] - rowStart[rB0];
    int64_t m = A.rows;
    int64_t n = widthB;
    GRID_ASSERT( widthA == k );
    GRID_ASSERT( n >= 1 );

    int owner = ( me >= rB0 ) && ( me < rB1 ) && ( B.rows > 0 );
    int64_t myOff = 0;
    if ( owner )
    {
      myOff = rowStart[me] - rowStart[rB0];
    }

    if ( m > 0 )
    {
      GRID_ASSERT( colA + widthA <= A.cols );
      GRID_ASSERT( colC + widthB <= C.cols );
      GRID_ASSERT( C.rows == m );
    }

    nGatherGemm++;
    GRID_TRACE("GatherGemm");

    std::vector<ComplexD> &stage = stageBuf;
    if ( owner )
    {
      GRID_ASSERT( colB + widthB <= B.cols );
      tStage -= usecond();
      if ( stage.size() < (uint64_t)B.rows*n ) stage.resize((uint64_t)B.rows*n);
      acceleratorCopyFromDevice(B.ColumnWindow(colB), &stage[0],
                                (uint64_t)B.rows*n*sizeof(ComplexD));
      tStage += usecond();
    }

    int64_t kc = panelBytes / ( (int64_t)sizeof(ComplexD) * n );
    if ( kc < 1 ) kc = 1;
    if ( kc > k ) kc = k;
    GRID_ASSERT( kc*n < 2147483647L );   // GlobalSumVector count is int

    std::vector<ComplexD>  &panel  = panelBuf;
    deviceVector<ComplexD> &dPanel = dPanelBuf;
    if ( panel.size()  < (uint64_t)kc*n ) panel.resize((uint64_t)kc*n);
    if ( dPanel.size() < (uint64_t)kc*n ) dPanel.resize((uint64_t)kc*n);
    deviceVector<ComplexD*> ap(1);
    deviceVector<ComplexD*> bp(1);
    deviceVector<ComplexD*> cp(1);
    std::vector<ComplexD*>  ptr(1);

    for(int64_t k0=0; k0<k; k0+=kc)
    {
      int64_t kchunk = std::min(kc, k-k0);

      // PLANNED OPTIMISATION (not yet): single-threaded memset zero-fills
      // the WHOLE panel; owners then overwrite their segment.  A threaded
      // zero of only the non-owned rows (thread_for over columns, memset
      // per column run) halves the host traffic and parallelises it.
      // Deliberately deferred until the simple version is proven.
      tMemset -= usecond();
      memset(&panel[0], 0, (uint64_t)kchunk*n*sizeof(ComplexD));
      tMemset += usecond();
      if ( owner )
      {
        int64_t i0 = std::max(k0,        myOff);
        int64_t i1 = std::min(k0+kchunk, myOff+B.rows);
        if ( i1 > i0 )
        {
          int64_t len = i1-i0;
          tDeposit -= usecond();
          thread_for(j, n, {
            memcpy(&panel[(uint64_t)((i0-k0)    + j*kchunk)],
                   &stage[(uint64_t)((i0-myOff) + j*B.rows)],
                   len*sizeof(ComplexD));
          });
          tDeposit += usecond();
        }
      }
      double tar = -usecond();
      grid->GlobalSumVector(&panel[0], (int)(kchunk*n));
      tar += usecond();
      tAllreduce += tar;
      tARmin = std::min(tARmin, tar);
      tARmax = std::max(tARmax, tar);
      bytesAllreduce += (uint64_t)kchunk*n*sizeof(ComplexD);
      nAllreduce++;

      if ( m > 0 )
      {
        tH2D -= usecond();
        acceleratorCopyToDevice(&panel[0], &dPanel[0],
                                (uint64_t)kchunk*n*sizeof(ComplexD));
        tH2D += usecond();

        ComplexD beta_use = ( k0==0 ) ? beta : ComplexD(1.0,0.0);

        ptr[0] = A.ColumnWindow(colA + k0);
        acceleratorCopyToDevice(&ptr[0], &ap[0], sizeof(ComplexD*));
        ptr[0] = &dPanel[0];
        acceleratorCopyToDevice(&ptr[0], &bp[0], sizeof(ComplexD*));
        ptr[0] = C.ColumnWindow(colC);
        acceleratorCopyToDevice(&ptr[0], &cp[0], sizeof(ComplexD*));

        tGemm -= usecond();
        BLAS.gemmBatched(GridBLAS_OP_N, GridBLAS_OP_N,
                         (int)m, (int)n, (int)kchunk,
                         alpha,    ap, (int)A.ld,
                                   bp, (int)kchunk,
                         beta_use, cp, (int)C.ld);
        BLAS.synchronise();
        tGemm += usecond();
      }
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // Local Frobenius norm-squared of a full-height column window.
  // NO comms; callers GlobalSum the result.  Host staging, setup-scale.
  ///////////////////////////////////////////////////////////////////////////
  double FrobNorm2Local(BlockRows &X, int64_t col0, int64_t w)
  {
    if ( X.rows == 0 ) return 0.0;
    GRID_ASSERT( col0 + w <= X.cols );
    uint64_t len = (uint64_t)X.rows*w;
    std::vector<ComplexD> h(len);
    acceleratorCopyFromDevice(X.ColumnWindow(col0), &h[0], len*sizeof(ComplexD));
    // Member real()/imag(): portable across std::complex (CPU) and
    // thrust::complex (HIP), where std::norm does not resolve.
    double s = 0.0;
    for(uint64_t i=0; i<len; i++)
    {
      double re = h[i].real();
      double im = h[i].imag();
      s += re*re + im*im;
    }
    return s;
  }

  ///////////////////////////////////////////////////////////////////////////
  // dst(:, dcol0 : dcol0+w) = - src(:, 0:w).  Both operands have ld == rows
  // so full-height windows are contiguous: flat elementwise device copy.
  ///////////////////////////////////////////////////////////////////////////
  void NegateCopy(BlockRows &src, BlockRows &dst, int64_t dcol0, int64_t w)
  {
    GRID_ASSERT( src.rows == dst.rows );
    GRID_ASSERT( w <= src.cols );
    GRID_ASSERT( dcol0 + w <= dst.cols );
    if ( src.rows == 0 ) return;
    uint64_t len = (uint64_t)src.rows*w;
    ComplexD *s = &src.data[0];
    ComplexD *d = dst.ColumnWindow(dcol0);
    accelerator_for(i, len, 1, {
      d[i] = -s[i];
    });
  }

  ///////////////////////////////////////////////////////////////////////////
  // Leaf inversion.  Purely LOCAL -- the calling rank owns the whole
  // width x width leaf (width == my row count); no collectives, so the
  // SPMD walk stays uniform with other ranks doing nothing.  The window
  // is contiguous (ld == rows == width): invert IN PLACE via
  // GridBLASInverse.  Everything is already fp64; no promote/demote.
  ///////////////////////////////////////////////////////////////////////////
  void LeafInvert(int64_t col0, int64_t width, BlockRows &Arows)
  {
    GRID_TRACE("SchurLeaf");
    GRID_ASSERT( width == Arows.rows );
    GRID_ASSERT( col0 + width <= Arows.cols );
    int64_t  w   = width;
    uint64_t len = (uint64_t)w*w;
    tLeaf -= usecond();

    deviceVector<ComplexD*> bp(1);
    std::vector<ComplexD*>  ptr(1);
    ptr[0] = Arows.ColumnWindow(col0);
    acceleratorCopyToDevice(&ptr[0], &bp[0], sizeof(ComplexD*));
    INV.inverseBatched(w, bp);

    // Telemetry: max |element| of the leaf inverse
    {
      std::vector<ComplexD> h(len);
      acceleratorCopyFromDevice(Arows.ColumnWindow(col0), &h[0], len*sizeof(ComplexD));
      double mx = 0.0;
      for(uint64_t i=0; i<len; i++)
      {
        double re = h[i].real();
        double im = h[i].imag();
        mx = std::max(mx, re*re + im*im);
      }
      telLeafMaxInv = std::max(telLeafMaxInv, std::sqrt(mx));
    }
    tLeaf += usecond();
  }

  ///////////////////////////////////////////////////////////////////////////
  // The recursion (plan 3.5 / 4B.3).  Inverts the diagonal block of the
  // rank-major matrix spanned by ranks [r0, r1), living in every member
  // rank's column window [col0, col0+width) -- IN PLACE.
  //
  // SPMD: every rank calls with IDENTICAL (r0, r1, width) and its own
  // local (col0, Arows); ranks outside [r0, r1) participate in the
  // collectives only (dummy operands, zero contributions).  The collective
  // sequence -- 5 GatherGemm calls + 3 scalar GlobalSums per merge node --
  // is identical on every rank by construction.
  //
  //   I = [r0, mid)   J = [mid, r1)   widths WI, WJ
  //   1. recurse I:      A11 -> A11inv
  //   2. B = A11inv.A12                     (I rows)
  //   3. C = A21.A11inv                     (J rows)
  //   4. S = A22 - A21.B    in place        (J rows)   [alpha=-1, beta=1]
  //   5. recurse J:      S -> Sinv
  //   6. T = Sinv.C                         (J rows)
  //   7. U = B.Sinv                         (I rows)
  //   8. X11 = A11inv + U.C  in place       (I rows)   [beta=1]
  //   9. X12 = -U, X21 = -T  local negates; X22 = Sinv already in place
  ///////////////////////////////////////////////////////////////////////////
  void SchurNode(int r0, int r1, int64_t col0, int64_t width, BlockRows &Arows)
  {
    int span = r1 - r0;
    GRID_ASSERT( span >= 1 );
    GRID_ASSERT( width == rowStart[r1] - rowStart[r0] );

    if ( span == 1 )
    {
      if ( ( me == r0 ) && ( myNrows > 0 ) )
      {
        LeafInvert(col0, width, Arows);
      }
      return;
    }

    int     mid = ( r0 + r1 ) / 2;
    int64_t WI  = rowStart[mid] - rowStart[r0];
    int64_t WJ  = rowStart[r1]  - rowStart[mid];

    // Zero-width child ranges (all ranks of a half owning no rows) are a
    // KNOWN LIMITATION: fail loudly rather than divide mysteriously.
    GRID_ASSERT( WI > 0 );
    GRID_ASSERT( WJ > 0 );

    int inI = ( me >= r0  ) && ( me < mid );
    int inJ = ( me >= mid ) && ( me < r1  );

    ComplexD one  ( 1.0,0.0);
    ComplexD mone (-1.0,0.0);
    ComplexD zero ( 0.0,0.0);

    BlockRows dummy;

    // 1. A11 -> A11inv
    SchurNode(r0, mid, col0, WI, Arows);

    // 2. B = A11inv . A12   (I rows; gather A12 from I owners)
    BlockRows Bbuf;
    if ( inI ) Bbuf.Resize(myNrows, WJ);
    {
      BlockRows &Aop = inI ? Arows : dummy;
      BlockRows &Cop = inI ? Bbuf  : dummy;
      int64_t    cA  = inI ? col0  : 0;
      GatherGemm(one, Aop, cA, WI,
                 r0, mid,
                 Arows, col0+WI, WJ,
                 zero, Cop, 0);
    }
    double nB = FrobNorm2Local(Bbuf, 0, inI ? WJ : 0);
    grid->GlobalSumVector(&nB, 1);
    telNormB.push_back(std::sqrt(nB));

    // 3. C = A21 . A11inv   (J rows; gather A11inv from I owners)
    BlockRows Cbuf;
    if ( inJ ) Cbuf.Resize(myNrows, WI);
    {
      BlockRows &Aop = inJ ? Arows : dummy;
      BlockRows &Cop = inJ ? Cbuf  : dummy;
      int64_t    cA  = inJ ? col0  : 0;
      GatherGemm(one, Aop, cA, WI,
                 r0, mid,
                 Arows, col0, WI,
                 zero, Cop, 0);
    }

    // 4. S = A22 - A21 . B  in place on my A22 window (J rows)
    double nA22 = FrobNorm2Local( inJ ? Arows : dummy, inJ ? col0+WI : 0, inJ ? WJ : 0 );
    grid->GlobalSumVector(&nA22, 1);
    {
      BlockRows &Aop = inJ ? Arows   : dummy;
      BlockRows &Cop = inJ ? Arows   : dummy;
      int64_t    cA  = inJ ? col0    : 0;
      int64_t    cC  = inJ ? col0+WI : 0;
      GatherGemm(mone, Aop, cA, WI,
                 r0, mid,
                 Bbuf, 0, WJ,
                 one, Cop, cC);
    }
    double nS = FrobNorm2Local( inJ ? Arows : dummy, inJ ? col0+WI : 0, inJ ? WJ : 0 );
    grid->GlobalSumVector(&nS, 1);
    telSratio.push_back( std::sqrt(nS) / ( std::sqrt(nA22) + 1.0e-300 ) );

    // 5. S -> Sinv
    SchurNode(mid, r1, col0+WI, WJ, Arows);

    // 6. T = Sinv . C       (J rows; gather C from J owners)
    BlockRows Tbuf;
    if ( inJ ) Tbuf.Resize(myNrows, WI);
    {
      BlockRows &Aop = inJ ? Arows   : dummy;
      BlockRows &Cop = inJ ? Tbuf    : dummy;
      int64_t    cA  = inJ ? col0+WI : 0;
      GatherGemm(one, Aop, cA, WJ,
                 mid, r1,
                 Cbuf, 0, WI,
                 zero, Cop, 0);
    }

    // 7. U = B . Sinv       (I rows; gather Sinv from J owners)
    BlockRows Ubuf;
    if ( inI ) Ubuf.Resize(myNrows, WJ);
    {
      BlockRows &Aop = inI ? Bbuf : dummy;
      BlockRows &Cop = inI ? Ubuf : dummy;
      GatherGemm(one, Aop, 0, WJ,
                 mid, r1,
                 Arows, col0+WI, WJ,
                 zero, Cop, 0);
    }

    // 8. X11 = A11inv + U . C  in place (I rows; gather C from J owners)
    {
      BlockRows &Aop = inI ? Ubuf  : dummy;
      BlockRows &Cop = inI ? Arows : dummy;
      int64_t    cC  = inI ? col0  : 0;
      GatherGemm(one, Aop, 0, WJ,
                 mid, r1,
                 Cbuf, 0, WI,
                 one, Cop, cC);
    }

    // 9. Off-diagonal signs, local
    if ( inI ) NegateCopy(Ubuf, Arows, col0+WI, WJ);
    if ( inJ ) NegateCopy(Tbuf, Arows, col0,    WI);
  }

  ///////////////////////////////////////////////////////////////////////////
  // PUBLIC ENTRY.  Arows: my rows of the rank-major N x N matrix (fp64).
  // On exit Arows holds my rows of the inverse, still fp64; the caller
  // owns the single terminal rounding into its fp32 apply storage.
  ///////////////////////////////////////////////////////////////////////////
  void Invert(BlockRows &Arows)
  {
    GRID_ASSERT( Arows.rows == myNrows );
    GRID_ASSERT( Arows.cols == N );

    telNormB.resize(0);
    telSratio.resize(0);
    telLeafMaxInv = 0.0;

    tStage         = 0.0;
    tMemset        = 0.0;
    tDeposit       = 0.0;
    tAllreduce     = 0.0;
    tH2D           = 0.0;
    tGemm          = 0.0;
    tLeaf          = 0.0;
    tARmin         = 1.0e30;
    tARmax         = 0.0;
    bytesAllreduce = 0;
    nAllreduce     = 0;
    nGatherGemm    = 0;

    SchurNode(0, P, 0, N, Arows);

    RealD mx = telLeafMaxInv;
    grid->GlobalMax(mx);
    telLeafMaxInv = mx;
  }

  // All telemetry values are globally reduced or boss-local; safe to
  // stream on every rank (Grid quiesces stdout to the boss unless
  // --debug-stdout).  NOTE: tAllreduce INCLUDES wait/imbalance -- a rank
  // arriving early books its wait here; the min/max spread across ranks
  // separates true wire time (~min) from skew (max-min).
  void ReportTelemetry(void)
  {
    for(uint64_t i=0; i<telNormB.size(); i++)
    {
      std::cout << GridLogPerformance
                << "SchurNode " << i
                << " ||B||_F "       << telNormB[i]
                << " ||S||/||A22|| " << telSratio[i]
                << std::endl;
    }
    std::cout << GridLogPerformance
              << "Schur leaves max|Ainv| " << telLeafMaxInv
              << std::endl;

    RealD armax =  tAllreduce;
    RealD armin = -tAllreduce;
    grid->GlobalMax(armax);
    grid->GlobalMax(armin);
    armin = -armin;

    std::cout << GridLogMessage << "Schur phases (boss rank, seconds):"
              << "  stage "     << tStage/1.0e6
              << "  memset "    << tMemset/1.0e6
              << "  deposit "   << tDeposit/1.0e6
              << "  allreduce " << tAllreduce/1.0e6
              << "  H2D "       << tH2D/1.0e6
              << "  gemm "      << tGemm/1.0e6
              << "  leaf "      << tLeaf/1.0e6
              << std::endl;
    std::cout << GridLogMessage << "Schur comms:"
              << "  GatherGemm calls " << nGatherGemm
              << "  panel allreduces " << nAllreduce
              << "  allreduce GB "     << bytesAllreduce/1024./1024./1024.
              << "  allreduce s min/max over ranks " << armin/1.0e6
              << " / " << armax/1.0e6
              << std::endl;
    std::cout << GridLogMessage << "Schur comms per-call (boss):"
              << "  min " << tARmin/1.0e3 << " ms"
              << "  avg " << (nAllreduce ? tAllreduce/nAllreduce/1.0e3 : 0.0) << " ms"
              << "  max " << tARmax/1.0e3 << " ms"
              << std::endl;
  }
};

NAMESPACE_END(Grid);
