/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./Grid/algorithms/multigrid/BlockCyclicRedistribute.h

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

#include <Grid/algorithms/multigrid/BlockCyclicSumma.h>

NAMESPACE_BEGIN(Grid);

///////////////////////////////////////////////////////////////////////////////
// Stage 4 of the 2D distributed dense inverse: redistribution between the
// 1D rank-major row layout (BlockRows: rank r owns contiguous global rows
// [rowStart[r], rowStart[r+1]) of an N x N matrix, stored rows x N column
// major with ld = rows) and the 2D block-cyclic layout.
//
// This is what lets the EXISTING stencil->dense import, its certificate,
// the fp32 slab conversion and the apply path all remain byte-for-byte
// untouched: the 2D inverse slots between them as
//
//    RowsToCyclic -> BlockCyclicSchurInverse::Invert -> CyclicToRows
//
// Volume is one matrix pass each way -- N^2/P elements per rank (~1 GB at
// production), trivial against the inversion itself.
//
// Transport: PURE POINT-TO-POINT, like everything else in this stack.
// Ranks exchange in a round-robin TOURNAMENT (the circle method, on an odd
// modulus M so it covers every pair exactly once for any P, with byes):
// at round r, ranks x and y are partners iff x+y == r (mod M).  Each
// meeting handles both directed edges of the pair in ONE SendToRecvFrom,
// padded to the larger of the two edge sizes -- SendToRecvFrom carries a
// single byte count for both directions, and both endpoints compute the
// same max from the shared descriptors, so there is no asymmetric-size
// case and no zero-count shape.  Pairs with nothing to exchange skip the
// round, decided identically at both ends.
//
// Element enumeration within an edge is canonical -- ascending global
// column outer, ascending global row inner -- and each endpoint builds its
// OWN local offset tables from the shared descriptors, so no index data is
// ever transmitted.  The round trip is BITWISE exact (pure data movement,
// no arithmetic): Test_schur2d_redist proves it.
///////////////////////////////////////////////////////////////////////////////

class BlockCyclicRedistribute
{
public:
  /////////////////////////////////////////////////////////////////////////
  // The directed edge (1D rank r1, 2D rank r2): global rows of r1's range
  // whose row-block coordinate is r2's prow; ALL global columns whose
  // column-block coordinate is r2's pcol.  Every rank can enumerate any
  // edge from (rowStart, layout) alone.
  /////////////////////////////////////////////////////////////////////////
  static void EdgeRows(const std::vector<int64_t> &rowStart, int r1,
                       const BlockCyclicLayout &L, int r2,
                       std::vector<int64_t> &rows)
  {
    rows.clear();
    int p = r2 / L.Pc;                       // row-major rank convention
    for(int64_t i=rowStart[r1]; i<rowStart[r1+1]; i++)
      if ( (int)((i/L.nb) % L.Pr) == p ) rows.push_back(i);
  }
  static void EdgeCols(const BlockCyclicLayout &L, int r2,
                       std::vector<int64_t> &cols)
  {
    cols.clear();
    int q = r2 % L.Pc;
    for(int64_t b=0; b*L.nb<L.N; b++){
      if ( (int)(b % L.Pc) != q ) continue;
      int64_t g0=b*L.nb, g1=std::min(L.N,(b+1)*L.nb);
      for(int64_t j=g0;j<g1;j++) cols.push_back(j);
    }
  }

  /////////////////////////////////////////////////////////////////////////
  // Gather/scatter one edge between a matrix (device, column major, ld)
  // and a dense edge buffer, through device offset tables.
  //   buffer(a,b) = elem(rows[a], cols[b]),  a fastest.
  /////////////////////////////////////////////////////////////////////////
  static void MoveEdge(int toBuffer,
                       ComplexD *mat, int64_t ld,
                       const std::vector<int64_t> &roff,  // per-row offset in mat
                       const std::vector<int64_t> &coff,  // per-col offset in mat
                       ComplexD *buf)
  {
    int64_t nr = roff.size();
    int64_t nc = coff.size();
    if ( !(nr && nc) ) return;
    deviceVector<int64_t> dro(nr), dco(nc);
    acceleratorCopyToDevice((void *)&roff[0], (void *)&dro[0], nr*sizeof(int64_t));
    acceleratorCopyToDevice((void *)&coff[0], (void *)&dco[0], nc*sizeof(int64_t));
    int64_t *ro = &dro[0];
    int64_t *co = &dco[0];
    if ( toBuffer ) {
      accelerator_for(idx, (uint64_t)(nr*nc), 1, {
        int64_t b = idx / nr;
        int64_t a = idx - b*nr;
        buf[a + b*nr] = mat[ ro[a] + co[b]*ld ];
      });
    } else {
      accelerator_for(idx, (uint64_t)(nr*nc), 1, {
        int64_t b = idx / nr;
        int64_t a = idx - b*nr;
        mat[ ro[a] + co[b]*ld ] = buf[a + b*nr];
      });
    }
  }

  /////////////////////////////////////////////////////////////////////////
  // My offset tables for an edge, on whichever side I am.
  /////////////////////////////////////////////////////////////////////////
  static void Offsets1D(const std::vector<int64_t> &rows,
                        const std::vector<int64_t> &cols,
                        int64_t row0,
                        std::vector<int64_t> &roff, std::vector<int64_t> &coff)
  {
    roff.resize(rows.size());  coff.resize(cols.size());
    for(uint64_t a=0;a<rows.size();a++) roff[a] = rows[a]-row0;   // local row
    for(uint64_t b=0;b<cols.size();b++) coff[b] = cols[b];        // global col
  }
  static void Offsets2D(const BlockCyclicLayout &L,
                        const std::vector<int64_t> &rows,
                        const std::vector<int64_t> &cols,
                        std::vector<int64_t> &roff, std::vector<int64_t> &coff)
  {
    roff.resize(rows.size());  coff.resize(cols.size());
    for(uint64_t a=0;a<rows.size();a++){
      int p; int64_t l;
      BlockCyclicLayout::GlobalToLocal(rows[a], L.nb, L.Pr, p, l);
      GRID_ASSERT( p == L.prow );
      roff[a] = l;
    }
    for(uint64_t b=0;b<cols.size();b++){
      int q; int64_t l;
      BlockCyclicLayout::GlobalToLocal(cols[b], L.nb, L.Pc, q, l);
      GRID_ASSERT( q == L.pcol );
      coff[b] = l;
    }
  }

  /////////////////////////////////////////////////////////////////////////
  // The worker.  dir=+1 : 1D rows -> block cyclic ;  dir=-1 : back.
  /////////////////////////////////////////////////////////////////////////
  static void Redistribute(int dir, GridBase *grid,
                           const std::vector<int64_t> &rowStart,
                           ComplexD *rows1d, int64_t myrows,
                           BlockCyclicMatrix &A)
  {
    BlockCyclicLayout &L = A.layout;
    int P  = grid->ProcessorCount();
    int me = grid->ThisRank();
    GRID_ASSERT( (int)rowStart.size() == P+1 );
    GRID_ASSERT( rowStart[P] == L.N );
    GRID_ASSERT( rowStart[me+1]-rowStart[me] == myrows );
    int64_t row0 = rowStart[me];
    int64_t ld1  = myrows ? myrows : 1;

    std::vector<int64_t> rows, cols, roff, coff;
    deviceVector<ComplexD> sbuf(1), rbuf(1);

    ///////////////////////////////////////////////////////////////////////
    // Self edge: purely local, via a bounce buffer (shares all the code).
    ///////////////////////////////////////////////////////////////////////
    EdgeRows(rowStart, me, L, me, rows);
    EdgeCols(L, me, cols);
    if ( rows.size() && cols.size() ){
      uint64_t ne = rows.size()*cols.size();
      if ( sbuf.size() < ne ) sbuf.resize(ne);
      std::vector<int64_t> roff2, coff2;
      Offsets1D(rows, cols, row0, roff, coff);
      Offsets2D(L,  rows, cols, roff2, coff2);
      if ( dir > 0 ) {
        MoveEdge(1, rows1d,     ld1,     roff,  coff,  &sbuf[0]);
        MoveEdge(0, &A.data[0], L.mloc,  roff2, coff2, &sbuf[0]);
      } else {
        MoveEdge(1, &A.data[0], L.mloc,  roff2, coff2, &sbuf[0]);
        MoveEdge(0, rows1d,     ld1,     roff,  coff,  &sbuf[0]);
      }
    }

    ///////////////////////////////////////////////////////////////////////
    // Tournament over all pairs: odd modulus M, partner = (r - me) mod M.
    // Every unordered pair meets exactly once; partner==me or >=P is a bye.
    ///////////////////////////////////////////////////////////////////////
    int M = (P%2) ? P : P+1;
    for(int r=0;r<M;r++){
      int partner = (int)(((int64_t)r - me + 2L*M) % M);
      if ( partner == me || partner >= P ) continue;

      // outbound edge: my (dir>0 ? 1D rows : 2D data) -> partner
      // inbound  edge: partner -> my (dir>0 ? 2D data : 1D rows)
      std::vector<int64_t> orow, ocol, irow, icol;
      if ( dir > 0 ) { EdgeRows(rowStart, me,      L, partner, orow); EdgeCols(L, partner, ocol);
                       EdgeRows(rowStart, partner, L, me,      irow); EdgeCols(L, me,      icol); }
      else           { EdgeRows(rowStart, partner, L, me,      orow); EdgeCols(L, me,      ocol);
                       EdgeRows(rowStart, me,      L, partner, irow); EdgeCols(L, partner, icol); }

      uint64_t nout = orow.size()*ocol.size();
      uint64_t nin  = irow.size()*icol.size();
      if ( !(nout || nin) ) continue;              // both ends compute this identically

      uint64_t nmax = std::max(nout,nin);          // symmetric padded transfer
      if ( sbuf.size() < nmax ) sbuf.resize(nmax);
      if ( rbuf.size() < nmax ) rbuf.resize(nmax);

      if ( nout ){
        if ( dir > 0 ) { Offsets1D(orow, ocol, row0, roff, coff);
                         MoveEdge(1, rows1d, ld1, roff, coff, &sbuf[0]); }
        else           { Offsets2D(L, orow, ocol, roff, coff);
                         MoveEdge(1, &A.data[0], L.mloc, roff, coff, &sbuf[0]); }
      }
      grid->SendToRecvFrom((void *)&sbuf[0], partner,
                           (void *)&rbuf[0], partner,
                           nmax*sizeof(ComplexD));
      if ( nin ){
        if ( dir > 0 ) { Offsets2D(L, irow, icol, roff, coff);
                         MoveEdge(0, &A.data[0], L.mloc, roff, coff, &rbuf[0]); }
        else           { Offsets1D(irow, icol, row0, roff, coff);
                         MoveEdge(0, rows1d, ld1, roff, coff, &rbuf[0]); }
      }
    }
  }

  static void RowsToCyclic(GridBase *grid, const std::vector<int64_t> &rowStart,
                           ComplexD *rows1d, int64_t myrows, BlockCyclicMatrix &A)
  { Redistribute(+1, grid, rowStart, rows1d, myrows, A); }

  static void CyclicToRows(GridBase *grid, const std::vector<int64_t> &rowStart,
                           BlockCyclicMatrix &A, ComplexD *rows1d, int64_t myrows)
  { Redistribute(-1, grid, rowStart, rows1d, myrows, A); }
};

NAMESPACE_END(Grid);
