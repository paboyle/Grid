/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/lattice/Lattice_peekpoke.h

    Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: Peter Boyle <peterboyle@Peters-MacBook-Pro-2.local>

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
#ifndef GRID_LATTICE_PEEK_H
#define GRID_LATTICE_PEEK_H

///////////////////////////////////////////////
// Peeking and poking around
///////////////////////////////////////////////

NAMESPACE_BEGIN(Grid);


// FIXME accelerator_loop and accelerator_inline these
////////////////////////////////////////////////////////////////////////////////////////////////////
// Peek internal indices of a Lattice object
////////////////////////////////////////////////////////////////////////////////////////////////////
template<int Index,class vobj> 
auto PeekIndex(const Lattice<vobj> &lhs,int i) -> Lattice<decltype(peekIndex<Index>(vobj(),i))>
{
  Lattice<decltype(peekIndex<Index>(vobj(),i))> ret(lhs.Grid());
  ret.Checkerboard()=lhs.Checkerboard();
  autoView( ret_v, ret, AcceleratorWrite);
  autoView( lhs_v, lhs, AcceleratorRead);
  accelerator_for( ss, lhs_v.size(), 1, {
    ret_v[ss] = peekIndex<Index>(lhs_v[ss],i);
  });
  return ret;
};
template<int Index,class vobj> 
auto PeekIndex(const Lattice<vobj> &lhs,int i,int j) -> Lattice<decltype(peekIndex<Index>(vobj(),i,j))>
{
  Lattice<decltype(peekIndex<Index>(vobj(),i,j))> ret(lhs.Grid());
  ret.Checkerboard()=lhs.Checkerboard();
  autoView( ret_v, ret, AcceleratorWrite);
  autoView( lhs_v, lhs, AcceleratorRead);
  accelerator_for( ss, lhs_v.size(), 1, {
    ret_v[ss] = peekIndex<Index>(lhs_v[ss],i,j);
  });
  return ret;
};

////////////////////////////////////////////////////////////////////////////////////////////////////
// Poke internal indices of a Lattice object
////////////////////////////////////////////////////////////////////////////////////////////////////
template<int Index,class vobj>  
void PokeIndex(Lattice<vobj> &lhs,const Lattice<decltype(peekIndex<Index>(vobj(),0))> & rhs,int i)
{
  autoView( rhs_v, rhs, AcceleratorRead);
  autoView( lhs_v, lhs, AcceleratorWrite);
  accelerator_for( ss, lhs_v.size(), 1, {
    pokeIndex<Index>(lhs_v[ss],rhs_v[ss],i);
  });
}
template<int Index,class vobj> 
void PokeIndex(Lattice<vobj> &lhs,const Lattice<decltype(peekIndex<Index>(vobj(),0,0))> & rhs,int i,int j)
{
  autoView( rhs_v, rhs, AcceleratorRead);
  autoView( lhs_v, lhs, AcceleratorWrite);
  accelerator_for( ss, lhs_v.size(), 1, {
    pokeIndex<Index>(lhs_v[ss],rhs_v[ss],i,j);
  });
}

//////////////////////////////////////////////////////
// Poke a scalar object into the SIMD array
//////////////////////////////////////////////////////
template<class vobj,class sobj> 
void pokeSite(const sobj &s,Lattice<vobj> &l,const Coordinate &site){

  GridBase *grid=l.Grid();

  int Nsimd = grid->Nsimd();

  assert( l.Checkerboard()== l.Grid()->CheckerBoard(site));
  assert( sizeof(sobj)*Nsimd == sizeof(vobj));

  int rank,odx,idx;
  // Optional to broadcast from node 0.
  grid->GlobalCoorToRankIndex(rank,odx,idx,site);
  grid->Broadcast(grid->BossRank(),s);

  // extract-modify-merge cycle is easiest way and this is not perf critical
  ExtractBuffer<sobj> buf(Nsimd);
  autoView( l_v , l, CpuWrite);
  if ( rank == grid->ThisRank() ) {
    extract(l_v[odx],buf);
    buf[idx] = s;
    merge(l_v[odx],buf);
  }

  return;
};


//////////////////////////////////////////////////////////
// Peek a scalar object from the SIMD array
//////////////////////////////////////////////////////////
template<class vobj>
typename vobj::scalar_object peekSite(const Lattice<vobj> &l,const Coordinate &site){
  typename vobj::scalar_object s;
  peekSite(s,l,site);
  return s;
}        
template<class vobj,class sobj>
void peekSite(sobj &s,const Lattice<vobj> &l,const Coordinate &site){
        
  GridBase *grid=l.Grid();

  int Nsimd = grid->Nsimd();

  assert( l.Checkerboard() == l.Grid()->CheckerBoard(site));

  int rank,odx,idx;
  grid->GlobalCoorToRankIndex(rank,odx,idx,site);

  ExtractBuffer<sobj> buf(Nsimd);
  autoView( l_v , l, CpuWrite);
  extract(l_v[odx],buf);

  s = buf[idx];

  grid->Broadcast(rank,s);

  return;
};

//////////////////////////////////////////////////////////
// Peek a scalar object from the SIMD array
//////////////////////////////////////////////////////////
// Must be CPU read view
template<class vobj,class sobj>
inline void peekLocalSite(sobj &s,const LatticeView<vobj> &l,Coordinate &site)
{
  GridBase *grid = l.getGrid();
  assert(l.mode==CpuRead);
  typedef typename vobj::scalar_type scalar_type;
  typedef typename vobj::vector_type vector_type;

  int Nsimd = grid->Nsimd();

  //  assert( l.Checkerboard()== grid->CheckerBoard(site));
  assert( sizeof(sobj)*Nsimd == sizeof(vobj));

  static const int words=sizeof(vobj)/sizeof(vector_type);
  int odx,idx;
  idx= grid->iIndex(site);
  odx= grid->oIndex(site);
  
  const vector_type *vp = (const vector_type *) &l[odx];
  scalar_type * pt = (scalar_type *)&s;
      
  for(int w=0;w<words;w++){
    pt[w] = getlane(vp[w],idx);
  }
  //  std::cout << "peekLocalSite "<<site<<" "<<odx<<","<<idx<<" "<<s<<std::endl;
  return;
};
template<class vobj,class sobj>
inline void peekLocalSite(sobj &s,const Lattice<vobj> &l,Coordinate &site)
{
  autoView(lv,l,CpuRead);
  peekLocalSite(s,lv,site);
  return;
};

// Must be CPU write view
template<class vobj,class sobj>
inline void pokeLocalSite(const sobj &s,LatticeView<vobj> &l,Coordinate &site)
{
  GridBase *grid=l.getGrid();
  assert(l.mode==CpuWrite);

  typedef typename vobj::scalar_type scalar_type;
  typedef typename vobj::vector_type vector_type;

  int Nsimd = grid->Nsimd();

  //  assert( l.Checkerboard()== grid->CheckerBoard(site));
  assert( sizeof(sobj)*Nsimd == sizeof(vobj));

  static const int words=sizeof(vobj)/sizeof(vector_type);
  int odx,idx;
  idx= grid->iIndex(site);
  odx= grid->oIndex(site);

  vector_type * vp = (vector_type *)&l[odx];
  scalar_type * pt = (scalar_type *)&s;
  for(int w=0;w<words;w++){
    putlane(vp[w],pt[w],idx);
  }
  return;
};

template<class vobj,class sobj>
inline void pokeLocalSite(const sobj &s, Lattice<vobj> &l,Coordinate &site)
{
  autoView(lv,l,CpuWrite);
  pokeLocalSite(s,lv,site);
  return;
};

NAMESPACE_END(Grid);
#endif

