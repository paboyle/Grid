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
// The recursion splits the GLOBAL INDEX RANGE at the block boundary
// nearest the midpoint -- not the rank range as the retired 1D
// RecursiveSchurInverse did -- so every rank owns part of every sub-block
// at every depth, and the ownership gating (inI/inJ, dummy operands,
// zero-width rank ranges) a 1D scheme needs has no analogue here.
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
// GridBLASInverse, unpack -- with NO communication and no assembly.
// Successive leaves cycle over ranks, so leaf work is naturally spread.
//
// COMMUNICATION.  Every transfer in the whole inversion is a
// SendToRecvFrom inside BlockCyclicSumma's rings: pure point-to-point, no
// collectives on the critical path, deterministic summation order (so
// repeated inversions are bitwise identical).  ReportTelemetry() is the
// one optional exception: it performs reductions, and is only ever called
// explicitly by a caller who wants the numbers.
//
// NUMERICS.  No pivoting: every A11 and every Schur complement met on the
// way down must be non-singular.  The growth telemetry stands in for
// pivoting.
///////////////////////////////////////////////////////////////////////////////

class BlockCyclicSchurInverse
{
public:
  BlockCyclicSumma        SUMMA;
  GridBLASInverse         INV;

  // Telemetry: accumulated LOCALLY, no comms unless ReportTelemetry().
  double                  telLeafMaxInv;
  uint64_t                nLeaf;
  // BIG LEAVES.  Below span s blocks a sub-block lives on <= s of the Pr
  // process rows / s of the Pc columns; when s is small relative to the grid
  // the SUMMA rings run on a few ranks while the rest block in their next
  // SendToRecvFrom.  Instead: gather the (s*nb)^2 sub-block to one rank,
  // invert locally (inverseLU), scatter back.  Fires only while
  // span < min(Pr,Pc) -- when the whole grid participates the rings are not
  // degenerate and gathering would only concentrate memory (and at the top
  // of the tree would gather the whole matrix).  Larger s concentrates more
  // memory on the root; smaller s loses ring parallelism.
  int                     leafSpan = 9;
  uint64_t                nBigLeaf = 0;  int64_t maxBigW = 0;
  double                  tBigGather = 0, tBigInv = 0, tBigScatter = 0;
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
  void WindowCopyScale(DenseInverseScalar alpha,
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
    DenseInverseScalar *src = Src.LocalWindow(li0,lj0);
    DenseInverseScalar *dst = Dst.LocalWindow(li0,lj0);
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
    deviceVector<DenseInverseScalar> dense((uint64_t)w*w);
    {
      DenseInverseScalar *src = A.LocalWindow(lr0,lc0);
      DenseInverseScalar *dst = &dense[0];
      int64_t   ld  = L.mloc;
      accelerator_for(idx, (uint64_t)(w*w), 1, {
        int64_t jj = idx / w;
        int64_t ii = idx - jj*w;
        dst[ii + jj*w] = src[ii + jj*ld];
      });
    }
    {
      deviceVector<DenseInverseScalar*> bp(1);
      std::vector<DenseInverseScalar*>  ptr(1);
      ptr[0] = &dense[0];
      acceleratorCopyToDevice(&ptr[0], &bp[0], sizeof(DenseInverseScalar*));
      INV.inverseBatched(w, bp);
    }
    {
      DenseInverseScalar *src = &dense[0];
      DenseInverseScalar *dst = A.LocalWindow(lr0,lc0);
      int64_t   ld  = L.mloc;
      accelerator_for(idx, (uint64_t)(w*w), 1, {
        int64_t jj = idx / w;
        int64_t ii = idx - jj*w;
        dst[ii + jj*ld] = src[ii + jj*w];
      });
    }
    // Growth telemetry, local only.
    {
      std::vector<DenseInverseScalar> h((uint64_t)w*w);
      acceleratorCopyFromDevice(&dense[0], &h[0], h.size()*sizeof(DenseInverseScalar));
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
  // BIG LEAF on block range [b0,b1): gather to root = owner of block (b0,b0),
  // invert there, scatter back.  Rank q's piece of the sub-block is the
  // contiguous local window RowRange(c0,c1) x ColRange(c0,c1); its local row
  // ii maps to global row (brq0 + (ii/nb)*Pr)*nb + ii%nb where brq0 is q's
  // first block row >= b0 (closed form: no tables).  Transport is pairwise
  // SendToRecvFrom with root (symmetric byte count: the reverse direction
  // carries a same-size dummy -- a leaf-local cost, accepted for simplicity).
  ///////////////////////////////////////////////////////////////////////////
  static int64_t FirstBlock(int64_t b0, int p, int Pg){ int64_t r = ((b0 % Pg) <= p) ? b0 - (b0 % Pg) + p : b0 - (b0 % Pg) + Pg + p; return r; }
  void BigLeaf(BlockCyclicMatrix &A, int64_t b0, int64_t b1)
  {
    GRID_TRACE("SchurBigLeaf");
    BlockCyclicLayout &L = A.layout;
    GridBase *grid = A.grid;
    const int Pr=L.Pr, Pc=L.Pc, nb=(int)L.nb;
    GRID_ASSERT( L.me == L.prow*Pc + L.pcol );          // rank convention shared with the SUMMA rings
    const int64_t c0 = b0*L.nb, c1 = std::min(L.N, b1*L.nb), W = c1-c0;
    const int root = (int)((b0%Pr)*Pc + (b0%Pc));
    const int me   = L.me;
    nBigLeaf++; maxBigW = std::max(maxBigW, W);

    // my piece
    int64_t lr0,lr1,lc0,lc1; L.RowRange(c0,c1,lr0,lr1); L.ColRange(c0,c1,lc0,lc1);
    const int64_t mq = lr1-lr0, nq = lc1-lc0;

    deviceVector<DenseInverseScalar> dense;            // root only: W x W column major
    deviceVector<DenseInverseScalar> piece, dummy;     // piece: my mq x nq contiguous; dummy: reverse-direction filler
    if ( me == root ) dense.resize((uint64_t)W*W);

    auto pack_piece = [&](DenseInverseScalar *dst, int64_t m, int64_t n, int64_t r0, int64_t cc0){
      DenseInverseScalar *src = A.LocalWindow(r0,cc0); const int64_t ld = L.mloc;
      accelerator_for(idx,(uint64_t)(m*n),1,{ int64_t jj=idx/m, ii=idx-jj*m; dst[ii+jj*m] = src[ii+jj*ld]; });
    };
    auto unpack_piece = [&](DenseInverseScalar *src, int64_t m, int64_t n, int64_t r0, int64_t cc0){
      DenseInverseScalar *dst = A.LocalWindow(r0,cc0); const int64_t ld = L.mloc;
      accelerator_for(idx,(uint64_t)(m*n),1,{ int64_t jj=idx/m, ii=idx-jj*m; dst[ii+jj*ld] = src[ii+jj*m]; });
    };
    // root: piece of rank q  <->  dense, via the closed-form block map
    auto root_place = [&](DenseInverseScalar *pc, int64_t m, int64_t n, int q, int to_dense){
      const int pq=q/Pc, cq=q%Pc;
      const int64_t brq0=FirstBlock(b0,pq,Pr), bcq0=FirstBlock(b0,cq,Pc);
      DenseInverseScalar *dn = &dense[0]; const int64_t WW=W, NB=nb, PR=Pr, PC=Pc, C0=c0;
      accelerator_for(idx,(uint64_t)(m*n),1,{
        int64_t jj=idx/m, ii=idx-jj*m;
        int64_t gr = (brq0 + (ii/NB)*PR)*NB + ii%NB - C0;
        int64_t gc = (bcq0 + (jj/NB)*PC)*NB + jj%NB - C0;
        if ( to_dense ) dn[gr + gc*WW] = pc[ii+jj*m]; else pc[ii+jj*m] = dn[gr + gc*WW];
      });
    };
    auto piece_dims = [&](int q, int64_t &m, int64_t &n){
      const int pq=q/Pc, cq=q%Pc;
      m = BlockCyclicLayout::NumLocal(c1,L.nb,pq,Pr) - BlockCyclicLayout::NumLocal(c0,L.nb,pq,Pr);
      n = BlockCyclicLayout::NumLocal(c1,L.nb,cq,Pc) - BlockCyclicLayout::NumLocal(c0,L.nb,cq,Pc);
    };

    // ---- gather ----
    tBigGather -= usecond();
    if ( mq*nq ) { piece.resize((uint64_t)mq*nq); dummy.resize((uint64_t)mq*nq); pack_piece(&piece[0],mq,nq,lr0,lc0); accelerator_barrier(); }
    if ( me == root ) {
      deviceVector<DenseInverseScalar> stage;
      for(int q=0;q<Pr*Pc;q++){
        int64_t m,n; piece_dims(q,m,n); if ( !(m*n) ) continue;
        if ( q == root ) { root_place(&piece[0],m,n,q,1); continue; }
        if ( stage.size() < (uint64_t)(m*n) ) stage.resize((uint64_t)m*n);
        deviceVector<DenseInverseScalar> junk((uint64_t)m*n);
        grid->SendToRecvFrom((void *)&junk[0], q, (void *)&stage[0], q, (uint64_t)m*n*sizeof(DenseInverseScalar));
        root_place(&stage[0],m,n,q,1);
      }
      accelerator_barrier();
    } else if ( mq*nq ) {
      grid->SendToRecvFrom((void *)&piece[0], root, (void *)&dummy[0], root, (uint64_t)mq*nq*sizeof(DenseInverseScalar));
    }
    tBigGather += usecond();

    // ---- invert on root: blocked getrf_64 + identity getrs_64 ----
    tBigInv -= usecond();
    if ( me == root ) INV.inverseLU(W, &dense[0]);
    tBigInv += usecond();

    // ---- scatter ----
    tBigScatter -= usecond();
    if ( me == root ) {
      deviceVector<DenseInverseScalar> stage;
      for(int q=0;q<Pr*Pc;q++){
        int64_t m,n; piece_dims(q,m,n); if ( !(m*n) ) continue;
        if ( q == root ) { root_place(&piece[0],m,n,q,0); accelerator_barrier(); continue; }
        if ( stage.size() < (uint64_t)(m*n) ) stage.resize((uint64_t)m*n);
        deviceVector<DenseInverseScalar> junk((uint64_t)m*n);
        root_place(&stage[0],m,n,q,0); accelerator_barrier();
        grid->SendToRecvFrom((void *)&stage[0], q, (void *)&junk[0], q, (uint64_t)m*n*sizeof(DenseInverseScalar));
      }
    } else if ( mq*nq ) {
      grid->SendToRecvFrom((void *)&dummy[0], root, (void *)&piece[0], root, (uint64_t)mq*nq*sizeof(DenseInverseScalar));
    }
    if ( mq*nq ) { unpack_piece(&piece[0],mq,nq,lr0,lc0); accelerator_barrier(); }
    tBigScatter += usecond();
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
    if ( span <= leafSpan && span < std::min(L.Pr,L.Pc) ) { BigLeaf(A, b0, b1); return; }
    GRID_TRACE("SchurNode");
    nNode++;

    int64_t bm = b0 + span/2;
    int64_t c0 = b0*L.nb;
    int64_t m  = bm*L.nb;
    int64_t c1 = std::min(L.N, b1*L.nb);
    DenseInverseScalar one (1.0,0.0), mone(-1.0,0.0), zero(0.0,0.0);

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
    if ( nBigLeaf ) {
      RealD ti = tBigInv/1.0e6; grid->GlobalMax(ti);     // inverse runs on the root of each leaf: report the max over ranks
      std::cout << GridLogMessage << "BlockCyclicSchurInverse: BIG LEAVES (span " << leafSpan
                << ", inverseLU): " << nBigLeaf
                << " leaves, max W " << maxBigW
                << "  boss secs: gather " << tBigGather/1.0e6 << " scatter " << tBigScatter/1.0e6
                << "  inverse (max over ranks) " << ti << std::endl;
    }
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
    // Ring time decomposed by message size (boss rank; every rank sends the
    // same sequence of sizes).  Time is wall time inside SendToRecvFrom, so it
    // includes waiting for the partner -- a bucket whose GB/s is far below the
    // probe's for the same size is wait, not wire.
    std::cout << GridLogMessage << "BlockCyclicSumma ring histogram (boss):  size-bucket  msgs  GB  xfer-secs  GB/s  %time" << std::endl;
    std::streamsize oldprec = std::cout.precision();
    for(int b=0;b<SUMMA.NHIST;b++){
      if ( !SUMMA.histN[b] ) continue;
      double sec = SUMMA.histUs[b]/1.0e6, g = SUMMA.histBytes[b]/1.0e9;
      double lo = (double)(1ull<<b);
      char sz[32]; if (lo>=1048576) snprintf(sz,32,"%6.1f MB",lo/1048576.); else if (lo>=1024) snprintf(sz,32,"%6.1f KB",lo/1024.); else snprintf(sz,32,"%6.0f B ",lo);
      std::cout << GridLogMessage << "   >=" << sz
                << std::setw(8) << SUMMA.histN[b]
                << std::setw(10) << std::setprecision(3) << g
                << std::setw(9) << std::setprecision(3) << sec
                << std::setw(9) << std::setprecision(3) << (sec>0 ? g/sec : 0.0)
                << std::setw(8) << std::setprecision(3) << (ring>0 ? 100.0*sec/ring : 0.0);
      std::cout << std::endl;
    }
    std::cout.precision(oldprec);
  }
};

NAMESPACE_END(Grid);
