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
  // BIG LEAVES (SCHUR2D_LEAF_SPAN=s, default 1 = off).  Below span s blocks a
  // sub-block lives on <= s of the Pr process rows / s of the Pc columns; the
  // SUMMA rings then run on a few ranks while the rest block in their next
  // SendToRecvFrom (histogram 2026-08-27: 93% of ring time in the 3.7 MB
  // single-block panels of exactly these levels).  Instead: gather the
  // (s*nb)^2 sub-block to one rank, invert locally, scatter back.
  int                     leafSpan = -1;
  int                     leafLU   = -1;  // SCHUR2D_LEAF_LU=1: big-leaf inverse via GridBLASInverse::inverseLU (blocked getrf_64 + identity getrs_64) instead of getri_batched
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

    deviceVector<ComplexD> dense;            // root only: W x W column major
    deviceVector<ComplexD> piece, dummy;     // piece: my mq x nq contiguous; dummy: reverse-direction filler
    if ( me == root ) dense.resize((uint64_t)W*W);

    auto pack_piece = [&](ComplexD *dst, int64_t m, int64_t n, int64_t r0, int64_t cc0){
      ComplexD *src = A.LocalWindow(r0,cc0); const int64_t ld = L.mloc;
      accelerator_for(idx,(uint64_t)(m*n),1,{ int64_t jj=idx/m, ii=idx-jj*m; dst[ii+jj*m] = src[ii+jj*ld]; });
    };
    auto unpack_piece = [&](ComplexD *src, int64_t m, int64_t n, int64_t r0, int64_t cc0){
      ComplexD *dst = A.LocalWindow(r0,cc0); const int64_t ld = L.mloc;
      accelerator_for(idx,(uint64_t)(m*n),1,{ int64_t jj=idx/m, ii=idx-jj*m; dst[ii+jj*ld] = src[ii+jj*m]; });
    };
    // root: piece of rank q  <->  dense, via the closed-form block map
    auto root_place = [&](ComplexD *pc, int64_t m, int64_t n, int q, int to_dense){
      const int pq=q/Pc, cq=q%Pc;
      const int64_t brq0=FirstBlock(b0,pq,Pr), bcq0=FirstBlock(b0,cq,Pc);
      ComplexD *dn = &dense[0]; const int64_t WW=W, NB=nb, PR=Pr, PC=Pc, C0=c0;
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
      deviceVector<ComplexD> stage;
      for(int q=0;q<Pr*Pc;q++){
        int64_t m,n; piece_dims(q,m,n); if ( !(m*n) ) continue;
        if ( q == root ) { root_place(&piece[0],m,n,q,1); continue; }
        if ( stage.size() < (uint64_t)(m*n) ) stage.resize((uint64_t)m*n);
        deviceVector<ComplexD> junk((uint64_t)m*n);
        grid->SendToRecvFrom((void *)&junk[0], q, (void *)&stage[0], q, (uint64_t)m*n*sizeof(ComplexD));
        root_place(&stage[0],m,n,q,1);
      }
      accelerator_barrier();
    } else if ( mq*nq ) {
      grid->SendToRecvFrom((void *)&piece[0], root, (void *)&dummy[0], root, (uint64_t)mq*nq*sizeof(ComplexD));
    }
    tBigGather += usecond();

    // ---- invert on root ----
    tBigInv -= usecond();
    if ( leafLU < 0 ) leafLU = getenv("SCHUR2D_LEAF_LU") ? atoi(getenv("SCHUR2D_LEAF_LU")) : 0;
    if ( me == root ) {
      if ( leafLU ) {
        INV.inverseLU(W, &dense[0]);
      } else {
        deviceVector<ComplexD*> bp(1); std::vector<ComplexD*> ptr(1); ptr[0] = &dense[0];
        acceleratorCopyToDevice(&ptr[0], &bp[0], sizeof(ComplexD*));
        INV.inverseBatched(W, bp);
      }
    }
    tBigInv += usecond();

    // ---- scatter ----
    tBigScatter -= usecond();
    if ( me == root ) {
      deviceVector<ComplexD> stage;
      for(int q=0;q<Pr*Pc;q++){
        int64_t m,n; piece_dims(q,m,n); if ( !(m*n) ) continue;
        if ( q == root ) { root_place(&piece[0],m,n,q,0); accelerator_barrier(); continue; }
        if ( stage.size() < (uint64_t)(m*n) ) stage.resize((uint64_t)m*n);
        deviceVector<ComplexD> junk((uint64_t)m*n);
        root_place(&stage[0],m,n,q,0); accelerator_barrier();
        grid->SendToRecvFrom((void *)&stage[0], q, (void *)&junk[0], q, (uint64_t)m*n*sizeof(ComplexD));
      }
    } else if ( mq*nq ) {
      grid->SendToRecvFrom((void *)&dummy[0], root, (void *)&piece[0], root, (uint64_t)mq*nq*sizeof(ComplexD));
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
    if ( leafSpan < 0 ) leafSpan = getenv("SCHUR2D_LEAF_SPAN") ? atoi(getenv("SCHUR2D_LEAF_SPAN")) : 1;
    if ( span <= leafSpan ) { BigLeaf(A, b0, b1); return; }
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
  ///////////////////////////////////////////////////////////////////////////
  // POWER - CLOCK - GROUND.  Before the inverse runs, print the preconditions
  // that have differed between fast (62 s) and slow (133-141 s) runs of the
  // SAME inverse, and measure ONE SendToRecvFrom to the ACTUAL ring partners
  // at three sizes, device and host buffers.  Same code in the harness and in
  // production, so the two processes are compared on the identical primitive
  // before any explanation of the SUMMA rings is attempted.
  // SCHUR2D_PROBE=0 disables (costs ~0.1-0.5 s).
  ///////////////////////////////////////////////////////////////////////////
  void Probe(BlockCyclicMatrix &A)
  {
    BlockCyclicLayout &L = A.layout;
    GridBase *grid = A.grid;
    int me = grid->ThisRank();
    // --- banner ---
    int thr = -1;
#ifdef GRID_COMMS_MPI3
    MPI_Query_thread(&thr);
#endif
    const char *omp = getenv("OMP_NUM_THREADS");
    MemoryStatus ms = MemoryManager::GetFootprint();
    std::cout << GridLogMessage << "Schur2D PROBE banner: MPI thread level " << thr
              << " (0 single,1 funneled,2 serialized,3 multiple)  OMP_NUM_THREADS=" << (omp?omp:"unset")
              << "  MemoryManager device bytes " << ms.DeviceBytes/1.0e9 << " GB (LRU " << ms.DeviceLRUBytes/1.0e9
              << " GB, cap " << ms.DeviceMaxBytes/1.0e9 << " GB)"
              << "  grid " << L.Pr << "x" << L.Pc << " nb " << L.nb << std::endl;
#ifdef GRID_HIP
    if ( me==0 ) acceleratorMem();
#endif
    // --- ring partners exactly as SUMMA uses them ---
    int prow=L.prow, pcol=L.pcol, Pr=L.Pr, Pc=L.Pc;
    struct Ring { const char *name; int dest, src; };
    Ring rings[2] = { {"ringA(row, q+-1)", prow*Pc + (pcol+1)%Pc,       prow*Pc + (pcol-1+Pc)%Pc},
                      {"ringB(col, p+-1)", ((prow+1)%Pr)*Pc + pcol,     ((prow-1+Pr)%Pr)*Pc + pcol} };
    // 2/3/4 MB added 2026-08-27: the SUMMA histogram put 93% of ring time in
    // [2,4) MB messages at 0.3 GB/s while >=4 MB ran at 11-13 GB/s.
    const int NSZ = 6;
    uint64_t sizes[NSZ] = { 64ull*1024, 1024ull*1024, 2048ull*1024, 3072ull*1024, 4096ull*1024, 8ull*1024*1024 };
    uint64_t maxb = sizes[NSZ-1];
    deviceVector<char> dsend(maxb), drecv(maxb);
    std::vector<char>  hsend(maxb), hrecv(maxb);
    for(int r=0;r<2;r++){
      if ( (r==0 && Pc==1) || (r==1 && Pr==1) ) continue;
      int off = grid->IsOffNode(rings[r].dest);
      for(int si=0;si<NSZ;si++){
        uint64_t bytes = sizes[si];
        // warm one, time five, both memory spaces
        grid->SendToRecvFrom(&dsend[0], rings[r].dest, &drecv[0], rings[r].src, bytes);
        double t0=usecond();
        for(int i=0;i<5;i++) grid->SendToRecvFrom(&dsend[0], rings[r].dest, &drecv[0], rings[r].src, bytes);
        double td=(usecond()-t0)/5.0;
        grid->SendToRecvFrom(&hsend[0], rings[r].dest, &hrecv[0], rings[r].src, bytes);
        t0=usecond();
        for(int i=0;i<5;i++) grid->SendToRecvFrom(&hsend[0], rings[r].dest, &hrecv[0], rings[r].src, bytes);
        double th=(usecond()-t0)/5.0;
        // spread over ranks
        RealD dmax=td, dmin=-td; grid->GlobalMax(dmax); grid->GlobalMax(dmin); dmin=-dmin;
        std::cout << GridLogMessage << "Schur2D PROBE " << rings[r].name << (off?" OFF-node":" on-node")
                  << " " << bytes/1024 << " KB: device " << td << " us (" << bytes/td/1.0e3 << " GB/s) [min/max over ranks " << dmin << "/" << dmax << " us]"
                  << "  host " << th << " us (" << bytes/th/1.0e3 << " GB/s)" << std::endl;
      }
    }
    ///////////////////////////////////////////////////////////////////////
    // ONE-GCD LEAF MICROBENCHMARK (boss only).  The big-leaf inverse at
    // W=4320 measured 0.53 s with BOTH getri_batched and getrf_64+getrs_64
    // (2026-08-27) -- ~0.4 TF/s on a GCD that runs zgemm at ~15.  Time the
    // three primitives in isolation on a well-conditioned dense matrix so the
    // leaf's rate can be compared with the machine's, and getrf split from
    // getrs.  Sizes: the W of span 4 / 9 / 18 leaves on accelerator builds;
    // tiny on CPU builds (Eigen would take minutes at 4320).
    ///////////////////////////////////////////////////////////////////////
    if ( me == 0 ) {
#if defined(GRID_HIP) || defined(GRID_CUDA) || defined(GRID_SYCL)
      std::vector<int64_t> Ws({1920, 4320, 8640});
#else
      std::vector<int64_t> Ws({240, 480});
#endif
      for(int64_t W : Ws){
        deviceVector<ComplexD> M((uint64_t)W*W), C((uint64_t)W*W);
        { ComplexD *m = &M[0]; const int64_t WW = W;   // diagonally dominant: (i==j ? W : 0) + cos/sin noise
          accelerator_for(idx,(uint64_t)W*W,1,{ int64_t j=idx/WW, i=idx-j*WW; double x=0.37*i+0.61*j;
            m[idx] = ComplexD((i==j)?(double)WW:0.0,0.0) + ComplexD(std::cos(x),std::sin(1.3*x)); });
          accelerator_barrier(); }
        double flopLU = 8.0/3.0*(double)W*W*W;          // complex LU  ~ (4 real flops per complex mult-add) * (2/3 n^3)
        double flopGEMM = 8.0*(double)W*W*W;            // complex GEMM
        // 1. getri_batched (batch 1)
        double tb;
        { deviceVector<ComplexD*> bp(1); std::vector<ComplexD*> ptr(1); ptr[0]=&M[0];
          acceleratorCopyToDevice(&ptr[0],&bp[0],sizeof(ComplexD*));
          double t0=usecond(); INV.inverseBatched(W,bp); tb=usecond()-t0; }
        // 2. inverseLU (getrf_64 + identity getrs_64), timed inside
        double tl; { double t0=usecond(); INV.inverseLU(W,&M[0]); tl=usecond()-t0; }
        // 3. one zgemm W x W x W for the machine rate
        double tg;
        { deviceVector<ComplexD*> ap(1),bp(1),cp(1); std::vector<ComplexD*> ptr(1);
          ptr[0]=&M[0]; acceleratorCopyToDevice(&ptr[0],&ap[0],sizeof(ComplexD*)); acceleratorCopyToDevice(&ptr[0],&bp[0],sizeof(ComplexD*));
          ptr[0]=&C[0]; acceleratorCopyToDevice(&ptr[0],&cp[0],sizeof(ComplexD*));
          double t0=usecond();
          SUMMA.BLAS.gemmBatched(GridBLAS_OP_N,GridBLAS_OP_N,(int)W,(int)W,(int)W,ComplexD(1.0,0.0),ap,(int)W,bp,(int)W,ComplexD(0.0,0.0),cp,(int)W);
          SUMMA.BLAS.synchronise(); tg=usecond()-t0; }
        std::cout << GridLogMessage << "Schur2D PROBE leaf W=" << W
                  << ": getri_batched " << tb/1.0e6 << " s (" << flopLU/tb/1.0e6 << " TF/s)"
                  << "  getrf_64+getrs_64 " << tl/1.0e6 << " s (getrf " << INV.lastGetrfUs/1.0e6 << " getrs " << INV.lastGetrsUs/1.0e6 << ")"
                  << "  zgemm " << tg/1.0e6 << " s (" << flopGEMM/tg/1.0e6 << " TF/s)" << std::endl;
      }
    }
    ///////////////////////////////////////////////////////////////////////
    // The SUMMA's conditions, one at a time, at 8 MB on ring B:
    //  (a) LARGE persistent buffers (the rings use ~0.5 GB Abuf/Bbuf), sending
    //      from offset 0 and from deep inside the region;
    //  (b) a pack kernel + accelerator_barrier immediately before each
    //      message, as the SUMMA does.
    // Isolated 8 MB messages ran at 11-20 GB/s while the SUMMA averaged 2.1;
    // whichever variant drops to ~2 GB/s names the condition.
    ///////////////////////////////////////////////////////////////////////
    {
      int r = (Pr>1) ? 1 : 0;
      uint64_t bytes = sizes[NSZ-1];
      uint64_t big = 512ull*1024*1024;
      deviceVector<char> bsend(big), brecv(big);
      for(int variant=0; variant<3; variant++){
        uint64_t so = (variant==1) ? big-bytes : 0;               // deep offset in the large region
        char *sp=&bsend[so], *rp=&brecv[so];
        grid->SendToRecvFrom(sp, rings[r].dest, rp, rings[r].src, bytes);
        double t0=usecond();
        for(int i=0;i<5;i++){
          if ( variant==2 ) { accelerator_for(k, bytes/8, 1, { ((uint64_t *)sp)[k] = (uint64_t)k; }); accelerator_barrier(); }
          grid->SendToRecvFrom(sp, rings[r].dest, rp, rings[r].src, bytes);
        }
        double t=(usecond()-t0)/5.0;
        RealD tmax=t, tmin=-t; grid->GlobalMax(tmax); grid->GlobalMax(tmin); tmin=-tmin;
        const char *vn[3]={"512MB buffer, offset 0","512MB buffer, offset 504MB","pack kernel + barrier before each send"};
        std::cout << GridLogMessage << "Schur2D PROBE " << rings[r].name << " 8192 KB device, " << vn[variant] << ": "
                  << t << " us (" << bytes/t/1.0e3 << " GB/s) [min/max over ranks " << tmin << "/" << tmax << " us]" << std::endl;
      }
    }
  }

  void Invert(BlockCyclicMatrix &A)
  {
    BlockCyclicLayout &L = A.layout;
    GRID_ASSERT( L.N >= 1 );
    int64_t nblocks = (L.N + L.nb - 1)/L.nb;
    if ( !(getenv("SCHUR2D_PROBE") && atoi(getenv("SCHUR2D_PROBE"))==0) ) Probe(A);

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
      std::cout << GridLogMessage << "BlockCyclicSchurInverse: BIG LEAVES (SCHUR2D_LEAF_SPAN=" << leafSpan
                << (leafLU>0 ? ", SCHUR2D_LEAF_LU=1: getrf_64+getrs_64" : ", getri_batched") << "): " << nBigLeaf
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
    std::cout << GridLogMessage << "BlockCyclicSumma ring histogram (boss):  size-bucket  msgs  GB  xfer-secs  GB/s  %time"
              << (SUMMA.handshake>0 ? "  handshake-secs (partner wait, excluded from xfer)" : "") << std::endl;
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
      if ( SUMMA.handshake>0 ) std::cout << std::setw(12) << std::setprecision(3) << SUMMA.histHsUs[b]/1.0e6;
      std::cout << std::endl;
    }
    std::cout.precision(oldprec);
  }
};

NAMESPACE_END(Grid);
