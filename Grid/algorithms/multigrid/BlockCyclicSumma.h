/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./Grid/algorithms/multigrid/BlockCyclicSumma.h

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

#include <Grid/algorithms/multigrid/BlockCyclic.h>

NAMESPACE_BEGIN(Grid);

///////////////////////////////////////////////////////////////////////////////
// Stage 2 of the 2D distributed dense inverse: the windowed SUMMA product
//
//   C[i0:i1, j0:j1]  <-  beta C[i0:i1, j0:j1]
//                      + alpha A[i0:i1, k0:k1] . B[k0:k1, j0:j1]
//
// on block-cyclic matrices sharing one BlockCyclicLayout.  All ranges are
// BLOCK ALIGNED (multiples of nb, or N itself at the top end): stage 3's
// recursion splits on block boundaries, so nothing else is ever needed, and
// alignment makes every local window a contiguous band of local storage
// (BlockCyclic.h RangeToLocal, Test_blockcyclic T7).
//
// Transport is PURE POINT-TO-POINT: SendToRecvFrom on explicit world ranks
// computed from the row-major rank convention.  No collectives of any kind
// -- no sub-communicators, no broadcast, no allgather -- by design: the
// collective pathologies measured on this machine (MPI_Allgatherv at
// ~0.18 MB/s with skewed counts, mpir_request.h:508 aborts, allreduce at
// 39% of the P2P rate) motivated this implementation.  SendToRecvFrom is
// the most exercised device-buffer path in Grid and the only one never
// implicated.
//
// Algorithm: round-based ring allgather SUMMA.  The k range is processed in
// rounds of Pc consecutive global blocks.  Within a round
//
//   * process column c owns at most one A panel  (blocks s with s%Pc == c);
//     the Pc panels circulate around each process-ROW ring in Pc-1 steps;
//   * process row r owns up to ceil(Pc/Pr) B panels (blocks s%Pr == r);
//     they circulate around each process-COLUMN ring in Pr-1 steps;
//   * every rank then accumulates  Cloc += alpha * Apanel_s . Bpanel_s
//     for each block s of the round, in ascending s: a fixed summation
//     order, so REPEATED RUNS ARE BITWISE IDENTICAL (no reduction, no
//     order ambiguity -- the property the P2P doctrine buys).
//
// Ring chunks are PADDED to a fixed size (full nb panels, fixed
// panels-per-origin): SendToRecvFrom carries one byte count for both
// directions, so symmetric transfers eliminate every variable-size edge
// case at a worst-case ~1/Pc extra volume on ragged rounds.  Padding is
// never read: GEMMs address only the leading nb_s x width of each slot.
//
// Per-rank received volume: (k-extent) * (mloc_i + nloc_j) elements --
// the N^2 (1/Pr + 1/Pc) SUMMA optimum, ~sqrt(P)/2 below the 1D scheme.
///////////////////////////////////////////////////////////////////////////////

class BlockCyclicMatrix
{
public:
  GridBase              *grid;      // borrowed, never owned
  BlockCyclicLayout      layout;
  deviceVector<ComplexD> data;      // column major, ld = layout.mloc

  BlockCyclicMatrix(GridBase *g, int64_t N, int64_t nb, int Pr, int Pc)
    : grid(g),
      layout(N, nb, Pr, Pc, g->ThisRank())
  {
    GRID_ASSERT( Pr*Pc == g->ProcessorCount() );
    uint64_t sz = (uint64_t)layout.mloc*layout.nloc;
    data.resize( sz ? sz : 1 );
  }

  ComplexD *LocalWindow(int64_t li, int64_t lj)
  {
    return &data[0] + li + lj*layout.mloc;
  }

  /////////////////////////////////////////////////////////////////////////
  // TEST-SCALE import/export of a replicated global matrix (host, O(N^2)
  // loops, one collective in Export).  For unit tests and the stage-3
  // oracle only; production data enters through the direct block-cyclic
  // import, never through these.
  /////////////////////////////////////////////////////////////////////////
  void ImportGlobal(const std::vector<ComplexD> &G)
  {
    int64_t N = layout.N;
    GRID_ASSERT( (int64_t)G.size() == N*N );
    std::vector<ComplexD> h((uint64_t)layout.mloc*layout.nloc, ComplexD(0.0,0.0));
    for(int64_t j=0;j<N;j++){
      for(int64_t i=0;i<N;i++){
        if ( layout.Owns(i,j) ) h[layout.LocalOffset(i,j)] = G[i + j*N];
      }
    }
    if ( h.size() )
      acceleratorCopyToDevice(&h[0], &data[0], h.size()*sizeof(ComplexD));
  }
  void ExportGlobal(std::vector<ComplexD> &G)
  {
    int64_t N = layout.N;
    G.assign((uint64_t)N*N, ComplexD(0.0,0.0));
    std::vector<ComplexD> h((uint64_t)layout.mloc*layout.nloc);
    if ( h.size() )
      acceleratorCopyFromDevice(&data[0], &h[0], h.size()*sizeof(ComplexD));
    for(int64_t j=0;j<N;j++){
      for(int64_t i=0;i<N;i++){
        if ( layout.Owns(i,j) ) G[i + j*N] = h[layout.LocalOffset(i,j)];
      }
    }
    if ( N ) grid->GlobalSumVector((ComplexD *)&G[0], (int)(N*N)); // zero-fill: exact
  }
};

class BlockCyclicSumma
{
public:
  GridBLAS BLAS;

  // Per-rank telemetry (boss-rank seconds when printed; no reductions here).
  // Every Multiply is: buffer alloc, pack panels, ring A along the process
  // row, ring B along the process column, then the local GEMMs.  The rings
  // are synchronous SendToRecvFrom, so tRing is time the GPU is idle unless
  // a future version overlaps them with the GEMMs.
  // PERSISTENT ring buffers: allocated once (grow-only) and reused by every
  // Multiply, so the device addresses handed to MPI never change.  Fresh
  // per-call buffers rotated through the caching allocator's blocks, and a
  // GPU-direct RDMA registration cache keyed on address then re-registers
  // per message: measured 1.26 GB/s/rank on 13 MB ring messages (production,
  // GRID_ALLOC_NCACHE_LARGE=64) against 62 s total for the same inverse when
  // hipMalloc returned a stable address.
  deviceVector<ComplexD> Abuf;
  deviceVector<ComplexD> Bbuf;
  double   tAlloc=0, tPack=0, tRingA=0, tRingB=0, tGemm=0;
  uint64_t bytesRing=0, nRingMsg=0, nMultiply=0, nGemm=0;
  // Per-message-size histogram (bucket = floor(log2 bytes)): count, bytes,
  // microseconds -- decomposes the ring time by packet size so a low average
  // GB/s can be attributed (many small latency-bound messages vs slow large
  // ones vs partner-wait).  The 8 MB probe runs at 11-20 GB/s; SUMMA averaged 2.
  static const int NHIST=48;
  uint64_t histN[NHIST]={0}, histBytes[NHIST]={0}; double histUs[NHIST]={0}, histHsUs[NHIST]={0};
  // SUMMA_HANDSHAKE=1: a 4-byte SendToRecvFrom with the same partner
  // immediately before each ring message, timed separately (histHsUs).
  // Handshake time = partner-arrival skew; the remainder = transfer.  Splits
  // the [2,4) MB bucket's 12 ms/msg (2026-08-27) into wait vs wire.
  int handshake = -1; int hsTx=0, hsRx=0;
  void HistAdd(uint64_t bytes, double us, double hs=0.0){ int b=0; while((bytes>>b)>1) b++; histN[b]++; histBytes[b]+=bytes; histUs[b]+=us; histHsUs[b]+=hs; }
  void ResetTelemetry(void){ tAlloc=tPack=tRingA=tRingB=tGemm=0; bytesRing=nRingMsg=nMultiply=nGemm=0; for(int b=0;b<NHIST;b++){histN[b]=histBytes[b]=0; histUs[b]=histHsUs[b]=0;} }

  static int Overlap(int64_t a0,int64_t a1,int64_t b0,int64_t b1)
  { return (a0 < b1) && (b0 < a1); }

  void Multiply(ComplexD alpha,
                BlockCyclicMatrix &A,
                BlockCyclicMatrix &B,
                ComplexD beta,
                BlockCyclicMatrix &C,
                int64_t i0, int64_t i1,
                int64_t j0, int64_t j1,
                int64_t k0, int64_t k1)
  {
    BlockCyclicLayout &L = C.layout;
    GridBase *grid = C.grid;
    const int64_t N  = L.N;
    const int64_t nb = L.nb;
    const int     Pr = L.Pr,   Pc   = L.Pc;
    const int     prow = L.prow, pcol = L.pcol;

    ///////////////////////////////////////////////////////////////////////
    // Conformability: one layout, one communicator, aligned ranges.
    ///////////////////////////////////////////////////////////////////////
    auto same = [&](BlockCyclicLayout &X){
      GRID_ASSERT( X.N==N ); GRID_ASSERT( X.nb==nb );
      GRID_ASSERT( X.Pr==Pr ); GRID_ASSERT( X.Pc==Pc );
      GRID_ASSERT( X.me==L.me );
    };
    same(A.layout); same(B.layout);
    GRID_ASSERT( A.grid==grid ); GRID_ASSERT( B.grid==grid );
    auto aligned = [&](int64_t g0,int64_t g1){
      GRID_ASSERT( 0<=g0 ); GRID_ASSERT( g0<=g1 ); GRID_ASSERT( g1<=N );
      GRID_ASSERT( g0%nb==0 ); GRID_ASSERT( (g1%nb==0)||(g1==N) );
    };
    aligned(i0,i1); aligned(j0,j1); aligned(k0,k1);
    GRID_ASSERT( k1 > k0 );                    // pure scaling not supported here
    // In-place windows are legal only if the written window is disjoint
    // from anything read (stage 3 uses this; make violation loud).
    if ( &C==&A ) GRID_ASSERT( !Overlap(j0,j1,k0,k1) );
    if ( &C==&B ) GRID_ASSERT( !Overlap(i0,i1,k0,k1) );

    ///////////////////////////////////////////////////////////////////////
    // My local bands of the three windows (contiguous: T7).
    ///////////////////////////////////////////////////////////////////////
    int64_t li0,li1, lj0,lj1;
    L.RowRange(i0,i1, li0,li1);
    L.ColRange(j0,j1, lj0,lj1);
    const int64_t mloc_i = li1-li0;            // my rows of the i window
    const int64_t nloc_j = lj1-lj0;            // my cols of the j window

    const int64_t kb0 = k0/nb;
    const int64_t kb1 = (k1+nb-1)/nb;          // block-aligned or ==N: exact
    const int64_t S   = (Pc + Pr - 1)/Pr;      // max B panels per origin row

    ///////////////////////////////////////////////////////////////////////
    // Round buffers, padded to fixed slot sizes (see header comment).
    //   A: Pc slots of mloc_i x nb   (slot c = panel of the block owned by c)
    //   B: Pr slots of S x (nb x nloc_j)
    ///////////////////////////////////////////////////////////////////////
    const uint64_t slotA  = (uint64_t)mloc_i*nb;
    const uint64_t slotB1 = (uint64_t)nb*nloc_j;         // one panel
    const uint64_t slotB  = (uint64_t)S*slotB1;
    nMultiply++;
    if ( handshake < 0 ) handshake = getenv("SUMMA_HANDSHAKE") ? atoi(getenv("SUMMA_HANDSHAKE")) : 0;
    tAlloc -= usecond();
    if ( Abuf.size() < std::max<uint64_t>(slotA*Pc,1) ) Abuf.resize( std::max<uint64_t>(slotA*Pc,1) );
    if ( Bbuf.size() < std::max<uint64_t>(slotB*Pr,1) ) Bbuf.resize( std::max<uint64_t>(slotB*Pr,1) );
    tAlloc += usecond();

    deviceVector<ComplexD *> ap(1), bp(1), cp(1);
    std::vector<ComplexD *>  ptr(1);

    int firstblock = 1;
    for(int64_t r0=kb0; r0<kb1; r0+=Pc){       // rounds of Pc k-blocks
      int64_t r1 = std::min(kb1, r0+Pc);

      /////////////////////////////////////////////////////////////////////
      // Pack MY panels of this round into my origin slots.
      /////////////////////////////////////////////////////////////////////
      tPack -= usecond();
      { GRID_TRACE("SummaPack");
      for(int64_t s=r0; s<r1; s++){
        int64_t nb_s = L.BlockSize(s);
        if ( (int)(s%Pc) == pcol && mloc_i ){  // my A panel: block-column s
          int64_t lc0, lc1;
          L.ColRange(s*nb, std::min(N,(s+1)*nb), lc0, lc1);
          GRID_ASSERT( lc1-lc0 == nb_s );
          ComplexD *src = A.LocalWindow(li0, lc0);
          ComplexD *dst = &Abuf[0] + slotA*pcol;
          int64_t   ld  = A.layout.mloc;
          int64_t   m   = mloc_i;
          accelerator_for(idx, (uint64_t)(m*nb_s), 1, {
            int64_t jj = idx / m;
            int64_t ii = idx - jj*m;
            dst[ii + jj*m] = src[ii + jj*ld];
          });
        }
        if ( (int)(s%Pr) == prow && nloc_j ){  // my B panel: block-row s
          int64_t lr0, lr1;
          L.RowRange(s*nb, std::min(N,(s+1)*nb), lr0, lr1);
          GRID_ASSERT( lr1-lr0 == nb_s );
          int64_t idxs = (s - r0 - ((prow - r0%Pr + Pr) % Pr)) / Pr; // my panel # in round
          GRID_ASSERT( idxs >= 0 ); GRID_ASSERT( idxs < S );
          ComplexD *src = B.LocalWindow(lr0, lj0);
          ComplexD *dst = &Bbuf[0] + slotB*prow + slotB1*idxs;
          int64_t   ld  = B.layout.mloc;
          int64_t   nn  = nloc_j;
          accelerator_for(idx, (uint64_t)(nb_s*nn), 1, {
            int64_t jj = idx / nb_s;
            int64_t ii = idx - jj*nb_s;
            dst[ii + jj*nb_s] = src[ii + jj*ld];
          });
        }
      }
      accelerator_barrier();
      }
      tPack += usecond();

      /////////////////////////////////////////////////////////////////////
      // Ring allgather along my process ROW: Pc-1 symmetric steps.  At
      // step t send the slot of origin (pcol-t+1), receive origin (pcol-t).
      /////////////////////////////////////////////////////////////////////
      if ( Pc > 1 && slotA ){
        GRID_TRACE("SummaRingA");
        tRingA -= usecond();
        int dest = prow*Pc + (pcol+1)%Pc;
        int src  = prow*Pc + (pcol-1+Pc)%Pc;
        for(int t=1;t<Pc;t++){
          int cs = (pcol - t + 1 + Pc*Pc) % Pc;
          int cr = (pcol - t     + Pc*Pc) % Pc;
          double ths = 0.0, tm = usecond();
          if ( handshake ) { grid->SendToRecvFrom((void *)&hsTx, dest, (void *)&hsRx, src, sizeof(int)); ths = usecond()-tm; tm = usecond(); }
          grid->SendToRecvFrom((void *)(&Abuf[0]+slotA*cs), dest,
                               (void *)(&Abuf[0]+slotA*cr), src,
                               slotA*sizeof(ComplexD));
          HistAdd(slotA*sizeof(ComplexD), usecond()-tm, ths);
          bytesRing += slotA*sizeof(ComplexD); nRingMsg++;
        }
        tRingA += usecond();
      }
      /////////////////////////////////////////////////////////////////////
      // Ring allgather along my process COLUMN: Pr-1 symmetric steps.
      /////////////////////////////////////////////////////////////////////
      if ( Pr > 1 && slotB ){
        GRID_TRACE("SummaRingB");
        tRingB -= usecond();
        int dest = ((prow+1)%Pr)*Pc + pcol;
        int src  = ((prow-1+Pr)%Pr)*Pc + pcol;
        for(int t=1;t<Pr;t++){
          int rs = (prow - t + 1 + Pr*Pr) % Pr;
          int rr = (prow - t     + Pr*Pr) % Pr;
          double ths = 0.0, tm = usecond();
          if ( handshake ) { grid->SendToRecvFrom((void *)&hsTx, dest, (void *)&hsRx, src, sizeof(int)); ths = usecond()-tm; tm = usecond(); }
          grid->SendToRecvFrom((void *)(&Bbuf[0]+slotB*rs), dest,
                               (void *)(&Bbuf[0]+slotB*rr), src,
                               slotB*sizeof(ComplexD));
          HistAdd(slotB*sizeof(ComplexD), usecond()-tm, ths);
          bytesRing += slotB*sizeof(ComplexD); nRingMsg++;
        }
        tRingB += usecond();
      }

      /////////////////////////////////////////////////////////////////////
      // Local update, ascending s: fixed order, bitwise-reproducible.
      /////////////////////////////////////////////////////////////////////
      tGemm -= usecond();
      { GRID_TRACE("SummaGEMM");
      for(int64_t s=r0; s<r1; s++){
        int64_t nb_s = L.BlockSize(s);
        if ( !(mloc_i && nloc_j && nb_s) ) { firstblock = 0; continue; }
        int cA = (int)(s%Pc);
        int rB = (int)(s%Pr);
        int64_t idxs = (s - r0 - ((rB - r0%Pr + Pr) % Pr)) / Pr;
        ComplexD beta_use = firstblock ? beta : ComplexD(1.0,0.0);
        firstblock = 0;

        ptr[0] = &Abuf[0] + slotA*cA;
        acceleratorCopyToDevice(&ptr[0], &ap[0], sizeof(ComplexD *));
        ptr[0] = &Bbuf[0] + slotB*rB + slotB1*idxs;
        acceleratorCopyToDevice(&ptr[0], &bp[0], sizeof(ComplexD *));
        ptr[0] = C.LocalWindow(li0, lj0);
        acceleratorCopyToDevice(&ptr[0], &cp[0], sizeof(ComplexD *));

        BLAS.gemmBatched(GridBLAS_OP_N, GridBLAS_OP_N,
                         (int)mloc_i, (int)nloc_j, (int)nb_s,
                         alpha,    ap, (int)mloc_i,
                                   bp, (int)nb_s,
                         beta_use, cp, (int)C.layout.mloc);
        BLAS.synchronise();
        nGemm++;
      }
      }
      tGemm += usecond();
    }
  }
};

NAMESPACE_END(Grid);
