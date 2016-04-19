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

namespace Grid {

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Peek internal indices of a Lattice object
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    template<int Index,class vobj>
       auto PeekIndex(const Lattice<vobj> &lhs,int i) -> Lattice<decltype(peekIndex<Index>(lhs._odata[0],i))>
    {
      Lattice<decltype(peekIndex<Index>(lhs._odata[0],i))> ret(lhs._grid);
      ret.checkerboard=lhs.checkerboard;
PARALLEL_FOR_LOOP
        for(int ss=0;ss<lhs._grid->oSites();ss++){
	  ret._odata[ss] = peekIndex<Index>(lhs._odata[ss],i);
        }
        return ret;
    };
    template<int Index,class vobj>
       auto PeekIndex(const Lattice<vobj> &lhs,int i,int j) -> Lattice<decltype(peekIndex<Index>(lhs._odata[0],i,j))>
    {
      Lattice<decltype(peekIndex<Index>(lhs._odata[0],i,j))> ret(lhs._grid);
      ret.checkerboard=lhs.checkerboard;
PARALLEL_FOR_LOOP
        for(int ss=0;ss<lhs._grid->oSites();ss++){
	  ret._odata[ss] = peekIndex<Index>(lhs._odata[ss],i,j);
        }
        return ret;
    };

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Poke internal indices of a Lattice object
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    template<int Index,class vobj> 
    void PokeIndex(Lattice<vobj> &lhs,const Lattice<decltype(peekIndex<Index>(lhs._odata[0],0))> & rhs,int i)
    {
PARALLEL_FOR_LOOP
        for(int ss=0;ss<lhs._grid->oSites();ss++){
	  pokeIndex<Index>(lhs._odata[ss],rhs._odata[ss],i);
	}      
    }
    template<int Index,class vobj>
      void PokeIndex(Lattice<vobj> &lhs,const Lattice<decltype(peekIndex<Index>(lhs._odata[0],0,0))> & rhs,int i,int j)
    {
PARALLEL_FOR_LOOP
        for(int ss=0;ss<lhs._grid->oSites();ss++){
	  pokeIndex<Index>(lhs._odata[ss],rhs._odata[ss],i,j);
	}      
    }

    //////////////////////////////////////////////////////
    // Poke a scalar object into the SIMD array
    //////////////////////////////////////////////////////
    template<class vobj,class sobj>
    void pokeSite(const sobj &s,Lattice<vobj> &l,std::vector<int> &site){

      GridBase *grid=l._grid;

      typedef typename vobj::scalar_type scalar_type;
      typedef typename vobj::vector_type vector_type;

      int Nsimd = grid->Nsimd();

      assert( l.checkerboard== l._grid->CheckerBoard(site));
      assert( sizeof(sobj)*Nsimd == sizeof(vobj));

      int rank,odx,idx;
      // Optional to broadcast from node 0.
      grid->GlobalCoorToRankIndex(rank,odx,idx,site);
      grid->Broadcast(grid->BossRank(),s);

      std::vector<sobj> buf(Nsimd);

      // extract-modify-merge cycle is easiest way and this is not perf critical
      if ( rank == grid->ThisRank() ) {
	extract(l._odata[odx],buf);
	buf[idx] = s;
	merge(l._odata[odx],buf);
      }

      return;
    };


    //////////////////////////////////////////////////////////
    // Peek a scalar object from the SIMD array
    //////////////////////////////////////////////////////////
    template<class vobj,class sobj>
      void peekSite(sobj &s,const Lattice<vobj> &l,std::vector<int> &site){
        
      GridBase *grid=l._grid;

      typedef typename vobj::scalar_type scalar_type;
      typedef typename vobj::vector_type vector_type;

      int Nsimd = grid->Nsimd();

      assert( l.checkerboard == l._grid->CheckerBoard(site));

      // FIXME
      //      assert( sizeof(sobj)*Nsimd == sizeof(vobj));

      int rank,odx,idx;
      grid->GlobalCoorToRankIndex(rank,odx,idx,site);

      std::vector<sobj> buf(Nsimd);
      extract(l._odata[odx],buf);

      s = buf[idx];

      grid->Broadcast(rank,s);

      return;
    };


    //////////////////////////////////////////////////////////
    // Peek a scalar object from the SIMD array
    //////////////////////////////////////////////////////////
    template<class vobj,class sobj>
    void peekLocalSite(sobj &s,const Lattice<vobj> &l,std::vector<int> &site){
        
      GridBase *grid=l._grid;

      typedef typename vobj::scalar_type scalar_type;
      typedef typename vobj::vector_type vector_type;

      int Nsimd = grid->Nsimd();

      assert( l.checkerboard== l._grid->CheckerBoard(site));
      assert( sizeof(sobj)*Nsimd == sizeof(vobj));

      int odx,idx;
      idx= grid->iIndex(site);
      odx= grid->oIndex(site);

      std::vector<sobj> buf(Nsimd);

      extract(l._odata[odx],buf);
      
      s = buf[idx];

      return;
    };

    template<class vobj,class sobj>
    void pokeLocalSite(const sobj &s,Lattice<vobj> &l,std::vector<int> &site){

      GridBase *grid=l._grid;

      typedef typename vobj::scalar_type scalar_type;
      typedef typename vobj::vector_type vector_type;

      int Nsimd = grid->Nsimd();

      assert( l.checkerboard== l._grid->CheckerBoard(site));
      assert( sizeof(sobj)*Nsimd == sizeof(vobj));

      int odx,idx;
      idx= grid->iIndex(site);
      odx= grid->oIndex(site);

      std::vector<sobj> buf(Nsimd);

      // extract-modify-merge cycle is easiest way and this is not perf critical
      extract(l._odata[odx],buf);
      
      buf[idx] = s;

      merge(l._odata[odx],buf);

      return;
    };

}
#endif

