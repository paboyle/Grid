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
auto PeekIndex(const Lattice<vobj> &lhs,int i) -> Lattice<decltype(peekIndex<Index>(lhs[0],i))>
{
  Lattice<decltype(peekIndex<Index>(lhs[0],i))> ret(lhs._grid);
  ret.Checkerboard()=lhs.Checkerboard();
  cpu_loop( ss, lhs, {
    ret[ss] = peekIndex<Index>(lhs[ss],i);
  });
  return ret;
};
template<int Index,class vobj> 
auto PeekIndex(const Lattice<vobj> &lhs,int i,int j) -> Lattice<decltype(peekIndex<Index>(lhs[0],i,j))>
{
  Lattice<decltype(peekIndex<Index>(lhs[0],i,j))> ret(lhs._grid);
  ret.Checkerboard()=lhs.Checkerboard();
  cpu_loop( ss, lhs, {
    ret[ss] = peekIndex<Index>(lhs[ss],i,j);
  });
  return ret;
};

////////////////////////////////////////////////////////////////////////////////////////////////////
// Poke internal indices of a Lattice object
////////////////////////////////////////////////////////////////////////////////////////////////////
template<int Index,class vobj>  
void PokeIndex(Lattice<vobj> &lhs,const Lattice<decltype(peekIndex<Index>(lhs[0],0))> & rhs,int i)
{
  cpu_loop( ss, lhs, {
    pokeIndex<Index>(lhs[ss],rhs[ss],i);
  });
}
template<int Index,class vobj> 
void PokeIndex(Lattice<vobj> &lhs,const Lattice<decltype(peekIndex<Index>(lhs[0],0,0))> & rhs,int i,int j)
{
  cpu_loop( ss, lhs, {
    pokeIndex<Index>(lhs[ss],rhs[ss],i,j);
  });
}

//////////////////////////////////////////////////////
// Poke a scalar object into the SIMD array
//////////////////////////////////////////////////////
template<class vobj,class sobj> 
void pokeSite(const sobj &s,Lattice<vobj> &l,const std::vector<int> &site){

  GridBase *grid=l._grid;

  typedef typename vobj::scalar_type scalar_type;
  typedef typename vobj::vector_type vector_type;

  int Nsimd = grid->Nsimd();

  assert( l.Checkerboard()== l._grid->CheckerBoard(site));
  assert( sizeof(sobj)*Nsimd == sizeof(vobj));

  int rank,odx,idx;
  // Optional to broadcast from node 0.
  grid->GlobalCoorToRankIndex(rank,odx,idx,site);
  grid->Broadcast(grid->BossRank(),s);

  std::vector<sobj> buf(Nsimd);

  // extract-modify-merge cycle is easiest way and this is not perf critical
  if ( rank == grid->ThisRank() ) {
    extract(l[odx],buf);
    buf[idx] = s;
    merge(l[odx],buf);
  }

  return;
};


//////////////////////////////////////////////////////////
// Peek a scalar object from the SIMD array
//////////////////////////////////////////////////////////
template<class vobj,class sobj>
void peekSite(sobj &s,const Lattice<vobj> &l,const std::vector<int> &site){
        
  GridBase *grid=l._grid;

  typedef typename vobj::scalar_type scalar_type;
  typedef typename vobj::vector_type vector_type;

  int Nsimd = grid->Nsimd();

  assert( l.Checkerboard() == l._grid->CheckerBoard(site));

  int rank,odx,idx;
  grid->GlobalCoorToRankIndex(rank,odx,idx,site);

  std::vector<sobj> buf(Nsimd);
  extract(l[odx],buf);

  s = buf[idx];

  grid->Broadcast(rank,s);

  return;
};


//////////////////////////////////////////////////////////
// Peek a scalar object from the SIMD array
//////////////////////////////////////////////////////////
template<class vobj,class sobj>
void peekLocalSite(sobj &s,const Lattice<vobj> &l,std::vector<int> &site){
        
  GridBase *grid = l._grid;

  typedef typename vobj::scalar_type scalar_type;
  typedef typename vobj::vector_type vector_type;

  int Nsimd = grid->Nsimd();

  assert( l.Checkerboard()== l._grid->CheckerBoard(site));
  assert( sizeof(sobj)*Nsimd == sizeof(vobj));

  static const int words=sizeof(vobj)/sizeof(vector_type);
  int odx,idx;
  idx= grid->iIndex(site);
  odx= grid->oIndex(site);

  scalar_type * vp = (scalar_type *)&l[odx];
  scalar_type * pt = (scalar_type *)&s;
      
  for(int w=0;w<words;w++){
    pt[w] = vp[idx+w*Nsimd];
  }
      
  return;
};

template<class vobj,class sobj>
void pokeLocalSite(const sobj &s,Lattice<vobj> &l,std::vector<int> &site){

  GridBase *grid=l._grid;

  typedef typename vobj::scalar_type scalar_type;
  typedef typename vobj::vector_type vector_type;

  int Nsimd = grid->Nsimd();

  assert( l.Checkerboard()== l._grid->CheckerBoard(site));
  assert( sizeof(sobj)*Nsimd == sizeof(vobj));

  static const int words=sizeof(vobj)/sizeof(vector_type);
  int odx,idx;
  idx= grid->iIndex(site);
  odx= grid->oIndex(site);

  scalar_type * vp = (scalar_type *)&l[odx];
  scalar_type * pt = (scalar_type *)&s;
      
  for(int w=0;w<words;w++){
    vp[idx+w*Nsimd] = pt[w];
  }

  return;
};

NAMESPACE_END(Grid);
#endif

