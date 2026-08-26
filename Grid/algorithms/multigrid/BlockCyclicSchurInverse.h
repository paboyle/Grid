/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./Grid/algorithms/multigrid/BlockCyclicSchurInverse.h

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
#pragma once

#include <Grid/algorithms/blas/BatchedInverse.h>
#include <Grid/algorithms/multigrid/BlockCyclicSumma.h>

NAMESPACE_BEGIN(Grid);

///////////////////////////////////////////////////////////////////////////////
// Stage 3 of the 2D distributed dense inverse: the recursive Schur
// complement on a block-cyclic matrix, in place.
//
// The nine-step algebra is IDENTICAL to RecursiveSchurInverse (1D); what
// changes is the decomposition.  The recursion splits the GLOBAL INDEX
// RANGE at the block boundary nearest the midpoint -- not the rank range --
// so every rank owns part of every sub-block at every depth, and the
// ownership gating (inI/inJ, dummy operands, zero-width rank ranges) of the
// 1D scheme has no analogue here: it is simply gone.
//
//   I = [c0,m)   J = [m,c1)      (block-aligned, m the mid block boundary)
//   1. recurse I  :  A11 -> A11inv                        (in place)
//   2. Bt = A11inv . A12                                  (scratch, I x J)
//   3. Ct = A21 . A11inv                                  (scratch, J x I)
//   4. A22 -= A21 . Bt          == S                      (in place)
//   5. recurse J  :  S -> Sinv                            (in place)
//   6. Tt = Sinv . Ct                                     (scratch, J x I)
//   7. Ut = Bt . Sinv                                     (scratch, I x J)
//   8. A11 += Ut . Ct           == X11                    (in place)
//   9. A12 = -Ut ,  A21 = -Tt                             (window copies)
//
// SCRATCH SHARING.  Four full-size block-cyclic scratch matrices (Bt, Ct,
// Tt, Ut) serve the ENTIRE tree, used through windows.  This is safe at
// every depth because of a window-disjointness invariant:
//
//   * every temporary of a node has its row range in one half of the
//     node's window and its column range in the other (I x J or J x I);
//   * everything any DESCENDANT touches -- its A windows and its own
//     temporaries -- has BOTH ranges inside a single half (I x I during
//     step 1, J x J during step 5).
//
// Hence a descendant window and a live ancestor temporary always differ in
// at least one dimension by disjoint ranges.  Only Bt and Ct are live
// across the step-5 recursion (Tt, Ut are written after it), and both are
// covered by the invariant.
//
// LEAF.  A leaf is a single diagonal block, and block (b,b) of a
// block-cyclic layout lives ENTIRELY on rank (b%Pr, b%Pc).  The leaf
// inversion is therefore purely local -- pack the strided block dense,
// GridBLASInverse, unpack -- with NO communication and no assembly.  The
// leaf-assembly transport question of the 1D scheme does not arise.
// Successive leaves cycle over ranks, so leaf work is naturally spread.
//
// COMMUNICATION.  Every transfer in the whole inversion is a
// SendToRecvFrom inside BlockCyclicSumma's rings: pure point-to-point, no
// collectives on the critical path, deterministic summation order (so
// repeated inversions are bitwise identical).  ReportTelemetry() is the
// one optional exception: it performs reductions, and is only ever called
// explicitly by a caller who wants the numbers.
//
// NUMERICS.  No pivoting, exactly as the 1D scheme: every A11 and every
// Schur complement met on the way down must be non-singular.  The growth
// telemetry stands in for pivoting; note the recursion splits differently
// from the 1D rank-range tree, so DIFFERENT sub-blocks are inverted and
// telemetry values are NOT comparable with the 1D implementation's --
// re-baseline, do not compare.
///////////////////////////////////////////////////////////////////////////////

class BlockCyclicSchurInverse
{
public:
  BlockCyclicSumma        SUMMA;
  GridBLASInverse         INV;

  // Telemetry: accumulated LOCALLY, no comms unless ReportTelemetry().
  double                  telLeafMaxInv;
  uint64_t                nLeaf;
  uint64_t                nNode;
  double                  tLeaf;
  double                  tGemm;      // wall in Multiply calls (comms+gemm)
  double                  tCopy;

  BlockCyclicSchurInverse()
  {
    telLeafMaxInv = 0.0;
    nLeaf = nNode = 0;
    tLeaf = tGemm = tCopy = 0.0;
  }

  ///////////////////////////////////////////////////////////////////////////
  // Window copy-scale:  Dst[i0:i1, j0:j1] = alpha * Src[same window].
  // Both share one layout, so the local bands coincide; pure local kernel.
  ///////////////////////////////////////////////////////////////////////////
  void WindowCopyScale(ComplexD alpha,
                       BlockCyclicMatrix &Src, BlockCyclicMatrix &Dst,
                       int64_t i0, int64_t i1, int64_t j0, int64_t j1)
  {
    BlockCyclicLayout &L = Dst.layout;
    GRID_ASSERT( Src.layout.N==L.N && Src.layout.nb==L.nb );
    GRID_ASSERT( Src.layout.Pr==L.Pr && Src.layout.Pc==L.Pc );
    int64_t li0,li1, lj0,lj1;
    L.RowRange(i0,i1, li0,li1);
    L.ColRange(j0,j1, lj0,lj1);
    int64_t m = li1-li0, n = lj1-lj0;
    if ( !(m && n) ) return;
    ComplexD *src = Src.LocalWindow(li0,lj0);
    ComplexD *dst = Dst.LocalWindow(li0,lj0);
    int64_t ldS = Src.layout.mloc;
    int64_t ldD = L.mloc;
    tCopy -= usecond();
    accelerator_for(idx, (uint64_t)(m*n), 1, {
      int64_t jj = idx / m;
      int64_t ii = idx - jj*m;
      dst[ii + jj*ldD] = alpha*src[ii + jj*ldS];
    });
    tCopy += usecond();
  }

  ///////////////////////////////////////////////////////////////////////////
  // Leaf: single diagonal block (b,b), entirely on rank (b%Pr, b%Pc).
  // Local pack -> dense inverse -> unpack; every other rank does nothing
  // and needs no synchronisation: the next SUMMA's rings pair them up.
  ///////////////////////////////////////////////////////////////////////////
  void Leaf(BlockCyclicMatrix &A, int64_t b)
  {
    GRID_TRACE("SchurLeaf");
    BlockCyclicLayout &L = A.layout;
    nLeaf++;
    if ( (int)(b % L.Pr) != L.prow ) return;
    if ( (int)(b % L.Pc) != L.pcol ) return;

    tLeaf -= usecond();
    int64_t g0 = b*L.nb;
    int64_t g1 = std::min(L.N, g0+L.nb);
    int64_t w  = g1-g0;
    int64_t lr0,lr1, lc0,lc1;
    L.RowRange(g0,g1, lr0,lr1);
    L.ColRange(g0,g1, lc0,lc1);
    GRID_ASSERT( lr1-lr0 == w );
    GRID_ASSERT( lc1-lc0 == w );

    // Pack the strided block dense (inverseBatched assumes lda == w).
    deviceVector<ComplexD> dense((uint64_t)w*w);
    {
      ComplexD *src = A.LocalWindow(lr0,lc0);
      ComplexD *dst = &dense[0];
      int64_t   ld  = L.mloc;
      accelerator_for(idx, (uint64_t)(w*w), 1, {
        int64_t jj = idx / w;
        int64_t ii = idx - jj*w;
        dst[ii + jj*w] = src[ii + jj*ld];
      });
    }
    {
      deviceVector<ComplexD*> bp(1);
      std::vector<ComplexD*>  ptr(1);
      ptr[0] = &dense[0];
      acceleratorCopyToDevice(&ptr[0], &bp[0], sizeof(ComplexD*));
      INV.inverseBatched(w, bp);
    }
    {
      ComplexD *src = &dense[0];
      ComplexD *dst = A.LocalWindow(lr0,lc0);
      int64_t   ld  = L.mloc;
      accelerator_for(idx, (uint64_t)(w*w), 1, {
        int64_t jj = idx / w;
        int64_t ii = idx - jj*w;
        dst[ii + jj*ld] = src[ii + jj*w];
      });
    }
    // Growth telemetry, local only.
    {
      std::vector<ComplexD> h((uint64_t)w*w);
      acceleratorCopyFromDevice(&dense[0], &h[0], h.size()*sizeof(ComplexD));
      double mx = 0.0;
      for(auto &z : h){
        double re=z.real(), im=z.imag();
        mx = std::max(mx, re*re+im*im);
      }
      telLeafMaxInv = std::max(telLeafMaxInv, std::sqrt(mx));
    }
    tLeaf += usecond();
  }

  ///////////////////////////////////////////////////////////////////////////
  // The recursion, on global BLOCK range [b0,b1).  SPMD: every rank calls
  // with identical arguments; there is no ownership gating to get wrong.
  ///////////////////////////////////////////////////////////////////////////
  void SchurNode(BlockCyclicMatrix &A,
                 BlockCyclicMatrix &Bt, BlockCyclicMatrix &Ct,
                 BlockCyclicMatrix &Tt, BlockCyclicMatrix &Ut,
                 int64_t b0, int64_t b1)
  {
    BlockCyclicLayout &L = A.layout;
    int64_t span = b1-b0;
    GRID_ASSERT( span >= 1 );
    if ( span == 1 ) { Leaf(A, b0); return; }
    GRID_TRACE("SchurNode");
    nNode++;

    int64_t bm = b0 + span/2;
    int64_t c0 = b0*L.nb;
    int64_t m  = bm*L.nb;
    int64_t c1 = std::min(L.N, b1*L.nb);
    ComplexD one (1.0,0.0), mone(-1.0,0.0), zero(0.0,0.0);

    // 1. A11 -> A11inv
    SchurNode(A,Bt,Ct,Tt,Ut, b0,bm);

    tGemm -= usecond();
    // 2. Bt[I,J] = A11inv . A12
    SUMMA.Multiply(one,  A,  A,  zero, Bt, c0,m,  m,c1,  c0,m );
    // 3. Ct[J,I] = A21 . A11inv
    SUMMA.Multiply(one,  A,  A,  zero, Ct, m,c1,  c0,m,  c0,m );
    // 4. A22 -= A21 . Bt   (the Schur complement, in place)
    SUMMA.Multiply(mone, A,  Bt, one,  A,  m,c1,  m,c1,  c0,m );
    tGemm += usecond();

    // 5. S -> Sinv          (Bt, Ct live across this call: see invariant)
    SchurNode(A,Bt,Ct,Tt,Ut, bm,b1);

    tGemm -= usecond();
    // 6. Tt[J,I] = Sinv . Ct
    SUMMA.Multiply(one,  A,  Ct, zero, Tt, m,c1,  c0,m,  m,c1 );
    // 7. Ut[I,J] = Bt . Sinv
    SUMMA.Multiply(one,  Bt, A,  zero, Ut, c0,m,  m,c1,  m,c1 );
    // 8. A11 += Ut . Ct
    SUMMA.Multiply(one,  Ut, Ct, one,  A,  c0,m,  c0,m,  m,c1 );
    tGemm += usecond();

    // 9. Off-diagonal signs
    { GRID_TRACE("SchurCopy");
    WindowCopyScale(mone, Ut, A, c0,m,  m,c1);
    WindowCopyScale(mone, Tt, A, m,c1,  c0,m);
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // PUBLIC ENTRY.  In-place inverse of the whole matrix.  Scratch (4x the
  // matrix footprint) is allocated here and released on return.
  ///////////////////////////////////////////////////////////////////////////
  void Invert(BlockCyclicMatrix &A)
  {
    BlockCyclicLayout &L = A.layout;
    GRID_ASSERT( L.N >= 1 );
    int64_t nblocks = (L.N + L.nb - 1)/L.nb;

    BlockCyclicMatrix Bt(A.grid, L.N, L.nb, L.Pr, L.Pc);
    BlockCyclicMatrix Ct(A.grid, L.N, L.nb, L.Pr, L.Pc);
    BlockCyclicMatrix Tt(A.grid, L.N, L.nb, L.Pr, L.Pc);
    BlockCyclicMatrix Ut(A.grid, L.N, L.nb, L.Pr, L.Pc);

    telLeafMaxInv = 0.0;
    nLeaf = nNode = 0;
    tLeaf = tGemm = tCopy = 0.0;

    SchurNode(A, Bt,Ct,Tt,Ut, 0, nblocks);
  }

  ///////////////////////////////////////////////////////////////////////////
  // Optional, and the ONLY place any reduction happens: call it if you
  // want the numbers, never from Invert.
  ///////////////////////////////////////////////////////////////////////////
  void ReportTelemetry(GridBase *grid)
  {
    RealD mx = telLeafMaxInv;
    grid->GlobalMax(mx);
    std::cout << GridLogMessage << "BlockCyclicSchurInverse:"
              << "  nodes " << nNode << "  leaves " << nLeaf
              << "  max|leafinv| " << mx
              << "  (boss secs: gemm+comms " << tGemm/1.0e6
              << "  leaf " << tLeaf/1.0e6
              << "  copy " << tCopy/1.0e6 << ")"
              << std::endl;
    // SUMMA breakdown: boss-rank seconds, plus the ring/gemm time spread over
    // ranks (min/max) -- skew shows as max >> min.
    RealD ring = (SUMMA.tRingA+SUMMA.tRingB)/1.0e6;
    RealD rmin = -ring, rmax = ring; grid->GlobalMax(rmin); grid->GlobalMax(rmax); rmin = -rmin;   // no GlobalMin: max of negation
    RealD gmin = -SUMMA.tGemm/1.0e6, gmax = SUMMA.tGemm/1.0e6; grid->GlobalMax(gmin); grid->GlobalMax(gmax); gmin = -gmin;
    double gb  = SUMMA.bytesRing/1.0e9;
    std::cout << GridLogMessage << "BlockCyclicSumma:"
              << " multiplies " << SUMMA.nMultiply << " gemms " << SUMMA.nGemm << " ring msgs " << SUMMA.nRingMsg
              << " | boss secs: alloc " << SUMMA.tAlloc/1.0e6
              << " pack " << SUMMA.tPack/1.0e6
              << " ringA " << SUMMA.tRingA/1.0e6 << " ringB " << SUMMA.tRingB/1.0e6
              << " gemm " << SUMMA.tGemm/1.0e6
              << " | ring min/max over ranks " << rmin << "/" << rmax
              << " gemm min/max " << gmin << "/" << gmax
              << " | ring bytes/rank " << gb << " GB -> " << (ring>0 ? gb/ring : 0.0) << " GB/s/rank (boss)"
              << std::endl;
  }
};

NAMESPACE_END(Grid);
