/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./Grid/algorithms/multigrid/BlockCyclic.h

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

NAMESPACE_BEGIN(Grid);

///////////////////////////////////////////////////////////////////////////////
// BlockCyclicLayout: the index arithmetic of a 2D block-cyclic distribution
// of an N x N matrix over a Pr x Pc logical process grid with block size nb.
//
// This is stage 1 of the 2D distributed dense inverse
// (documentation/DistributedDenseInverse2D.tex).  It is deliberately
// COMMUNICATOR-FREE: every mapping is a static pure function of
// (N, nb, Pr, Pc), so the whole layout is exhaustively unit-testable on one
// rank with no MPI in the loop (Test_blockcyclic).  A thin instance layer
// binds a world rank to a grid coordinate and caches local extents.
//
// Conventions (fixed here, relied on by every later stage):
//
//  * Global block b of a dimension with Pg processes is owned by process
//    coordinate  b % Pg  (ScaLAPACK csrc=0), and is that process's local
//    block  b / Pg.
//  * Rank <-> grid coordinate is ROW MAJOR over the process grid:
//        rank = p*Pc + q ,  p = rank/Pc ,  q = rank%Pc .
//    The eventual ring transport must construct its neighbour tables with
//    the same convention.
//  * Local storage is COLUMN MAJOR with ld = mloc, matching BlockRows:
//    local element (i,j) lives at  data[i + j*mloc].
//  * The trailing partial block (N % nb != 0) belongs to the owner of the
//    last full-size block position; only that one block is short.
//
// Element (gi,gj) therefore lives on grid coordinate
//    ( (gi/nb) % Pr , (gj/nb) % Pc )
// at local coordinate
//    ( ((gi/nb)/Pr)*nb + gi%nb , ((gj/nb)/Pc)*nb + gj%nb ).
//
// Everything here is host-side integer arithmetic; nothing allocates.
///////////////////////////////////////////////////////////////////////////////

class BlockCyclicLayout
{
public:
  ///////////////////////////////////////////////////////////////////////////
  // Closest-to-square factorisation Pr*Pc == P with Pr <= Pc.
  // For fixed P the per-rank SUMMA volume  N^2 (1/Pr + 1/Pc)  is minimised
  // at the most square grid.  P=288 -> 16 x 18.
  ///////////////////////////////////////////////////////////////////////////
  static void ChooseProcessGrid(int P, int &Pr, int &Pc)
  {
    GRID_ASSERT(P >= 1);
    Pr = 1;
    for(int r=1; (int64_t)r*r <= (int64_t)P; r++)
      if ( P % r == 0 ) Pr = r;
    Pc = P / Pr;
  }

  ///////////////////////////////////////////////////////////////////////////
  // Number of rows (or columns) of a dimension of global extent N, block nb,
  // owned by process coordinate p of Pg.  ScaLAPACK "numroc", csrc=0.
  ///////////////////////////////////////////////////////////////////////////
  static int64_t NumLocal(int64_t N, int64_t nb, int p, int Pg)
  {
    GRID_ASSERT(N  >= 0);
    GRID_ASSERT(nb >= 1);
    GRID_ASSERT(p  >= 0);
    GRID_ASSERT(p  <  Pg);
    int64_t nblocks = N / nb;            // full blocks
    int64_t extra   = N % nb;            // trailing partial block
    int64_t full    = nblocks / Pg;      // full blocks everyone owns
    int64_t rem     = nblocks % Pg;      // coords [0,rem) own one more
    int64_t n = full*nb;
    if ( p <  (int)rem ) n += nb;        // an extra full block
    if ( p == (int)rem ) n += extra;     // the partial block, if any
    return n;
  }

  ///////////////////////////////////////////////////////////////////////////
  // Global index -> (owner coordinate, local index) in one dimension.
  ///////////////////////////////////////////////////////////////////////////
  static void GlobalToLocal(int64_t g, int64_t nb, int Pg,
                            int &owner, int64_t &loc)
  {
    GRID_ASSERT(g >= 0);
    int64_t b = g / nb;                  // global block
    owner = (int)(b % Pg);
    loc   = (b / Pg)*nb + (g % nb);
  }

  ///////////////////////////////////////////////////////////////////////////
  // (process coordinate, local index) -> global index in one dimension.
  // Inverse of GlobalToLocal on the owned set.
  ///////////////////////////////////////////////////////////////////////////
  static int64_t LocalToGlobal(int64_t l, int64_t nb, int p, int Pg)
  {
    GRID_ASSERT(l >= 0);
    int64_t lb = l / nb;                 // local block
    int64_t b  = lb*Pg + p;              // global block
    return b*nb + (l % nb);
  }

  ///////////////////////////////////////////////////////////////////////////
  // Instance layer: bind a rank of a Pr x Pc grid.
  ///////////////////////////////////////////////////////////////////////////
  int64_t N;         // global matrix dimension (square)
  int64_t nb;        // block size
  int     Pr, Pc;    // process grid
  int     me;        // world rank within the grid, row major
  int     prow, pcol;// my grid coordinate
  int64_t mloc,nloc; // my local extents; storage column major, ld = mloc

  BlockCyclicLayout(int64_t N_, int64_t nb_, int Pr_, int Pc_, int me_)
  {
    N  = N_;
    nb = nb_;
    Pr = Pr_;
    Pc = Pc_;
    me = me_;
    GRID_ASSERT( N  >= 0 );
    GRID_ASSERT( nb >= 1 );
    GRID_ASSERT( Pr >= 1 );
    GRID_ASSERT( Pc >= 1 );
    GRID_ASSERT( me >= 0 );
    GRID_ASSERT( me <  Pr*Pc );
    prow = me / Pc;                      // ROW MAJOR rank convention
    pcol = me % Pc;
    mloc = NumLocal(N, nb, prow, Pr);
    nloc = NumLocal(N, nb, pcol, Pc);
  }

  // Owning rank of global element (gi,gj), row-major rank convention.
  int OwnerRank(int64_t gi, int64_t gj) const
  {
    int pr,pc; int64_t li,lj;
    GlobalToLocal(gi, nb, Pr, pr, li);
    GlobalToLocal(gj, nb, Pc, pc, lj);
    return pr*Pc + pc;
  }

  // My local storage offset of global element (gi,gj).
  // The caller must know I own it; asserted, not assumed.
  int64_t LocalOffset(int64_t gi, int64_t gj) const
  {
    int pr,pc; int64_t li,lj;
    GlobalToLocal(gi, nb, Pr, pr, li);
    GlobalToLocal(gj, nb, Pc, pc, lj);
    GRID_ASSERT( pr == prow );
    GRID_ASSERT( pc == pcol );
    return li + lj*mloc;                 // column major, ld = mloc
  }

  // Do I own global element (gi,gj)?
  int Owns(int64_t gi, int64_t gj) const
  {
    return OwnerRank(gi,gj) == me;
  }

  ///////////////////////////////////////////////////////////////////////////
  // Block-aligned global range [g0,g1) -> my contiguous local range [l0,l1).
  //
  // For fixed owner p the local index is monotone in the global index, so a
  // coordinate's owned elements of ANY global range are contiguous in local
  // storage; and for a BLOCK-ALIGNED range the bounds are exactly
  // NumLocal(g0) and NumLocal(g1), because NumLocal(g,...) counts the owned
  // elements below g.  This is what lets a windowed product view the local
  // sub-matrix of a global window as &data[l0 + c0*mloc] with the SAME ld --
  // no gather, no copy.  Verified exhaustively in Test_blockcyclic T7.
  //
  // g0 must be a block multiple; g1 a block multiple or N itself.
  ///////////////////////////////////////////////////////////////////////////
  static void RangeToLocal(int64_t g0, int64_t g1,
                           int64_t N, int64_t nb, int p, int Pg,
                           int64_t &l0, int64_t &l1)
  {
    GRID_ASSERT( 0 <= g0 );
    GRID_ASSERT( g0 <= g1 );
    GRID_ASSERT( g1 <= N );
    GRID_ASSERT( g0 % nb == 0 );
    GRID_ASSERT( (g1 % nb == 0) || (g1 == N) );
    l0 = NumLocal(g0, nb, p, Pg);
    l1 = NumLocal(g1, nb, p, Pg);
  }

  // Instance forms, rows and columns of my own coordinate.
  void RowRange(int64_t g0, int64_t g1, int64_t &l0, int64_t &l1) const
  { RangeToLocal(g0,g1,N,nb,prow,Pr,l0,l1); }
  void ColRange(int64_t g0, int64_t g1, int64_t &l0, int64_t &l1) const
  { RangeToLocal(g0,g1,N,nb,pcol,Pc,l0,l1); }

  // Size of global block b (the trailing block may be short).
  int64_t BlockSize(int64_t b) const
  {
    int64_t lo = b*nb;
    GRID_ASSERT( lo < N );
    return std::min(N, lo+nb) - lo;
  }
};

NAMESPACE_END(Grid);
