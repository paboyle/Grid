/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./tests/debug/Test_blockcyclic.cc

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
// Regression gate for BlockCyclicLayout -- stage 1 of the 2D distributed
// dense inverse (documentation/DistributedDenseInverse2D.tex).
//
// The layout is pure index arithmetic, so this test is EXHAUSTIVE rather
// than statistical: every stage sweeps a battery of (N, nb, Pr, Pc)
// configurations chosen for their edge cases -- nb=1, nb>N, N%nb!=0,
// more processes than blocks, prime N, the production shape -- and checks
// every global index (T2/T3) or every element of an outer-product grid (T4)
// against a brute-force reference built by walking blocks.
//
// No communication: correct at mpirun -n 1, and identical at any rank count.
//
//   T1 : NumLocal partitions N (sums over coords; against brute force).
//   T2 : GlobalToLocal / LocalToGlobal round trip, every index.
//   T3 : local indices are dense [0,mloc): a bijection, not just a cover.
//   T4 : 2D ownership: every element has exactly one owner; per-rank counts
//        equal mloc*nloc; LocalOffset is a bijection onto [0,mloc*nloc).
//   T5 : block contiguity: within any owned global block, consecutive
//        global rows are consecutive local rows (what SUMMA panels rely on).
//   T6 : ChooseProcessGrid: exact factorisation, Pr<=Pc, most-square.
//////////////////////////////////////////////////////////////////////////////

#include <Grid/Grid.h>
#include <Grid/algorithms/multigrid/BlockCyclic.h>

using namespace Grid;

static int failures = 0;

static void Report(const std::string &name, bool pass, const std::string &detail="")
{
  std::cout << GridLogMessage << "  " << name << (pass ? "  PASS" : "  ** FAIL **");
  if ( detail.size() ) std::cout << "   " << detail;
  std::cout << std::endl;
  if ( !pass ) failures++;
}

// One-dimensional configurations: {N, nb, Pg}
struct Cfg1d { int64_t N; int64_t nb; int Pg; };
static const std::vector<Cfg1d> configs1d = {
  {  0,  1, 1},          // empty matrix
  {  1,  1, 1},          // minimal
  { 16,  4, 4},          // exact tiling
  { 17,  4, 4},          // trailing partial block
  { 16,  1, 4},          // nb=1: pure cyclic
  { 16, 32, 4},          // nb>N: one short block, most coords empty
  { 16,  4, 7},          // more processes than blocks
  { 97, 13, 5},          // prime N, awkward everything
  {138240, 480, 16},     // production N, row dimension of 16x18
  {138240, 480, 18},     // production N, column dimension of 16x18
  {138240, 512, 16},     // production N, non-dividing block
};

int main(int argc, char **argv)
{
  Grid_init(&argc, &argv);

  std::cout << GridLogMessage << "BlockCyclicLayout regression (communicator-free index arithmetic)" << std::endl;

  ////////////////////////////////////////////////////////////////////////
  // T1 : NumLocal partitions N, and agrees with brute force block-walking.
  ////////////////////////////////////////////////////////////////////////
  {
    bool ok = true;
    for(auto &c : configs1d){
      // brute force: walk global blocks, count elements per coordinate
      std::vector<int64_t> counts(c.Pg, 0);
      for(int64_t b=0; b*c.nb < c.N || (c.N==0 && b<0); b++){
        if ( b*c.nb >= c.N ) break;
        int64_t lo = b*c.nb;
        int64_t hi = std::min(c.N, lo+c.nb);
        counts[b % c.Pg] += hi-lo;
      }
      int64_t sum = 0;
      for(int p=0;p<c.Pg;p++){
        int64_t n = BlockCyclicLayout::NumLocal(c.N, c.nb, p, c.Pg);
        if ( n != counts[p] ) ok = false;
        sum += n;
      }
      if ( sum != c.N ) ok = false;
    }
    Report("T1  NumLocal partitions N, == brute force", ok);
  }

  ////////////////////////////////////////////////////////////////////////
  // T2 : round trip, every global index of every configuration.
  ////////////////////////////////////////////////////////////////////////
  {
    bool ok = true;
    for(auto &c : configs1d){
      int64_t sweep = std::min<int64_t>(c.N, 200000); // full N except production
      for(int64_t g=0; g<sweep; g++){
        int p; int64_t l;
        BlockCyclicLayout::GlobalToLocal(g, c.nb, c.Pg, p, l);
        if ( BlockCyclicLayout::LocalToGlobal(l, c.nb, p, c.Pg) != g ) ok = false;
      }
      // and the production tail, where the arithmetic could overflow or drift
      for(int64_t g=std::max<int64_t>(0,c.N-1000); g<c.N; g++){
        int p; int64_t l;
        BlockCyclicLayout::GlobalToLocal(g, c.nb, c.Pg, p, l);
        if ( BlockCyclicLayout::LocalToGlobal(l, c.nb, p, c.Pg) != g ) ok = false;
      }
    }
    Report("T2  GlobalToLocal <-> LocalToGlobal round trip", ok);
  }

  ////////////////////////////////////////////////////////////////////////
  // T3 : per coordinate, the local indices hit [0, NumLocal) exactly once.
  ////////////////////////////////////////////////////////////////////////
  {
    bool ok = true;
    for(auto &c : configs1d){
      if ( c.N > 4096 ) continue;        // dense bitmap check: small cfgs only
      for(int p=0;p<c.Pg;p++){
        int64_t n = BlockCyclicLayout::NumLocal(c.N, c.nb, p, c.Pg);
        std::vector<int> hit(n, 0);
        for(int64_t g=0; g<c.N; g++){
          int pp; int64_t l;
          BlockCyclicLayout::GlobalToLocal(g, c.nb, c.Pg, pp, l);
          if ( pp != p ) continue;
          if ( l < 0 || l >= n ) { ok = false; continue; }
          hit[l]++;
        }
        for(int64_t l=0;l<n;l++) if ( hit[l] != 1 ) ok = false;
      }
    }
    Report("T3  local indices dense and unique on [0,NumLocal)", ok);
  }

  ////////////////////////////////////////////////////////////////////////
  // T4 : 2D ownership.  Small grids, every element: exactly one owner,
  // owner counts equal mloc*nloc, LocalOffset bijective onto storage.
  ////////////////////////////////////////////////////////////////////////
  {
    bool ok = true;
    struct Cfg2d { int64_t N; int64_t nb; int Pr; int Pc; };
    std::vector<Cfg2d> cfgs = {
      { 16, 4, 2, 2},
      { 17, 4, 2, 3},                    // partial block, non-square grid
      { 23, 5, 3, 2},
      { 12, 2, 3, 4},                    // 12 ranks
      { 30, 7, 4, 2},
    };
    for(auto &c : cfgs){
      int P = c.Pr*c.Pc;
      std::vector<BlockCyclicLayout> L;
      for(int r=0;r<P;r++) L.push_back(BlockCyclicLayout(c.N,c.nb,c.Pr,c.Pc,r));
      // ownership count per rank, and per-rank storage bitmap
      std::vector<int64_t> owned(P,0);
      std::vector<std::vector<int>> slot(P);
      for(int r=0;r<P;r++) slot[r].assign(L[r].mloc*L[r].nloc, 0);
      for(int64_t i=0;i<c.N;i++){
        for(int64_t j=0;j<c.N;j++){
          int r = L[0].OwnerRank(i,j);
          if ( r < 0 || r >= P ) { ok=false; continue; }
          // every layout instance must agree on the owner
          if ( L[r].OwnerRank(i,j) != r ) ok = false;
          if ( !L[r].Owns(i,j) )          ok = false;
          owned[r]++;
          int64_t off = L[r].LocalOffset(i,j);
          if ( off < 0 || off >= (int64_t)slot[r].size() ) { ok=false; continue; }
          slot[r][off]++;
        }
      }
      int64_t tot=0;
      for(int r=0;r<P;r++){
        if ( owned[r] != L[r].mloc*L[r].nloc ) ok = false;
        for(auto h : slot[r]) if ( h != 1 ) ok = false;
        tot += owned[r];
      }
      if ( tot != c.N*c.N ) ok = false;
    }
    Report("T4  2D ownership: unique owner, counts == mloc*nloc, offsets bijective", ok);
  }

  ////////////////////////////////////////////////////////////////////////
  // T5 : block contiguity.  For every owned global block, consecutive
  // global rows are consecutive local rows: g and g+1 in the same block
  // must give l and l+1 on the same coordinate.  SUMMA panel extraction
  // (stage 2) sends whole local blocks as contiguous strides through the
  // column-major store; that only works if this holds.
  ////////////////////////////////////////////////////////////////////////
  {
    bool ok = true;
    for(auto &c : configs1d){
      if ( c.N == 0 ) continue;
      int64_t sweep = std::min<int64_t>(c.N-1, 100000);
      for(int64_t g=0; g<sweep; g++){
        if ( (g+1) % c.nb == 0 ) continue;   // block boundary: owner may change
        int p0,p1; int64_t l0,l1;
        BlockCyclicLayout::GlobalToLocal(g,   c.nb, c.Pg, p0, l0);
        BlockCyclicLayout::GlobalToLocal(g+1, c.nb, c.Pg, p1, l1);
        if ( p1 != p0     ) ok = false;
        if ( l1 != l0 + 1 ) ok = false;
      }
    }
    Report("T5  intra-block contiguity (SUMMA panel precondition)", ok);
  }

  ////////////////////////////////////////////////////////////////////////
  // T6 : ChooseProcessGrid.
  ////////////////////////////////////////////////////////////////////////
  {
    bool ok = true;
    for(int P : {1,2,3,4,6,8,12,16,17,64,96,144,256,288,512}){
      int Pr,Pc;
      BlockCyclicLayout::ChooseProcessGrid(P,Pr,Pc);
      if ( Pr*Pc != P ) ok = false;
      if ( Pr > Pc )    ok = false;
      // most-square: no divisor r with Pr < r <= sqrt(P)
      for(int r=Pr+1; (int64_t)r*r <= (int64_t)P; r++)
        if ( P % r == 0 ) ok = false;
    }
    int Pr,Pc;
    BlockCyclicLayout::ChooseProcessGrid(288,Pr,Pc);
    if ( !(Pr==16 && Pc==18) ) ok = false;
    Report("T6  ChooseProcessGrid exact, Pr<=Pc, most-square (288 -> 16x18)", ok);
  }

  ////////////////////////////////////////////////////////////////////////
  // T7 : RangeToLocal.  For every block-aligned range of every small
  // configuration: the owned global indices of [g0,g1) map exactly onto
  // local [l0,l1), contiguously and in order (brute force).
  ////////////////////////////////////////////////////////////////////////
  {
    bool ok = true;
    for(auto &c : configs1d){
      if ( c.N == 0 || c.N > 4096 ) continue;
      int64_t nblocks = (c.N + c.nb - 1)/c.nb;
      for(int p=0;p<c.Pg;p++){
        for(int64_t b0=0;b0<=nblocks;b0++){
          for(int64_t b1=b0;b1<=nblocks;b1++){
            int64_t g0 = b0*c.nb;
            int64_t g1 = std::min(c.N, b1*c.nb);
            if ( g0 > c.N ) continue;
            int64_t l0,l1;
            BlockCyclicLayout::RangeToLocal(g0,g1,c.N,c.nb,p,c.Pg,l0,l1);
            // brute force: owned globals in [g0,g1) in ascending order
            int64_t expect = l0;
            for(int64_t g=g0; g<g1; g++){
              int pp; int64_t l;
              BlockCyclicLayout::GlobalToLocal(g, c.nb, c.Pg, pp, l);
              if ( pp != p ) continue;
              if ( l != expect ) ok = false;   // contiguous, in order
              expect++;
            }
            if ( expect != l1 ) ok = false;    // count matches the bounds
          }
        }
      }
    }
    Report("T7  RangeToLocal contiguous, ordered, exact bounds", ok);
  }

  std::cout << GridLogMessage << (failures ? "Test_blockcyclic: FAILURES"
                                           : "Test_blockcyclic: ALL PASS") << std::endl;
  Grid_finalize();
  return failures ? 1 : 0;
}
