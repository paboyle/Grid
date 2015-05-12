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
    inline auto peekIndex(const Lattice<vobj> &lhs)
      -> Lattice<decltype(peekIndex<Index>(lhs._odata[0]))>
    {
      Lattice<decltype(peekIndex<Index>(lhs._odata[0]))> ret(lhs._grid);
PARALLEL_FOR_LOOP
        for(int ss=0;ss<lhs._grid->oSites();ss++){
            ret._odata[ss] = peekIndex<Index>(lhs._odata[ss]);
        }
        return ret;
    };
    template<int Index,class vobj>
      inline auto peekIndex(const Lattice<vobj> &lhs,int i)
      -> Lattice<decltype(peekIndex<Index>(lhs._odata[0],i))>
    {
      Lattice<decltype(peekIndex<Index>(lhs._odata[0],i))> ret(lhs._grid);
PARALLEL_FOR_LOOP
        for(int ss=0;ss<lhs._grid->oSites();ss++){
	  ret._odata[ss] = peekIndex<Index>(lhs._odata[ss],i);
        }
        return ret;
    };
    template<int Index,class vobj>
      inline auto peekIndex(const Lattice<vobj> &lhs,int i,int j)
      -> Lattice<decltype(peekIndex<Index>(lhs._odata[0],i,j))>
    {
      Lattice<decltype(peekIndex<Index>(lhs._odata[0],i,j))> ret(lhs._grid);
PARALLEL_FOR_LOOP
        for(int ss=0;ss<lhs._grid->oSites();ss++){
	  ret._odata[ss] = peekIndex<Index>(lhs._odata[ss],i,j);
        }
        return ret;
    };

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Poke internal indices of a Lattice object
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    template<int Index,class vobj> inline
    void pokeIndex(Lattice<vobj> &lhs,const Lattice<decltype(peekIndex<Index>(lhs._odata[0]))> & rhs)
    {
PARALLEL_FOR_LOOP
        for(int ss=0;ss<lhs._grid->oSites();ss++){
	  pokeIndex<Index>(lhs._odata[ss],rhs._odata[ss]);
	}      
    }
    template<int Index,class vobj> inline
    void pokeIndex(Lattice<vobj> &lhs,const Lattice<decltype(peekIndex<Index>(lhs._odata[0],0))> & rhs,int i)
    {
PARALLEL_FOR_LOOP
        for(int ss=0;ss<lhs._grid->oSites();ss++){
	  pokeIndex<Index>(lhs._odata[ss],rhs._odata[ss],i);
	}      
    }
    template<int Index,class vobj> inline
    void pokeIndex(Lattice<vobj> &lhs,const Lattice<decltype(peekIndex<Index>(lhs._odata[0],0,0))> & rhs,int i,int j)
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
      void peekSite(sobj &s,Lattice<vobj> &l,std::vector<int> &site){
        
      GridBase *grid=l._grid;

      typedef typename vobj::scalar_type scalar_type;
      typedef typename vobj::vector_type vector_type;

      int Nsimd = grid->Nsimd();

      assert( l.checkerboard== l._grid->CheckerBoard(site));
      assert( sizeof(sobj)*Nsimd == sizeof(vobj));

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
    void peekLocalSite(sobj &s,Lattice<vobj> &l,std::vector<int> &site){
        
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

