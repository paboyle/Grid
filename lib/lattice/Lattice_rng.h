    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/lattice/Lattice_rng.h

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: paboyle <paboyle@ph.ed.ac.uk>

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
#ifndef GRID_LATTICE_RNG_H
#define GRID_LATTICE_RNG_H

#include <random>
#include <Grid/sitmo_rng/sitmo_prng_engine.hpp>

namespace Grid {

  //http://nvlpubs.nist.gov/nistpubs/SpecialPublications/NIST.SP.800-90Ar1.pdf ?

  //////////////////////////////////////////////////////////////
  // Allow the RNG state to be less dense than the fine grid
  //////////////////////////////////////////////////////////////
  inline int RNGfillable(GridBase *coarse,GridBase *fine)
  {

    int rngdims = coarse->_ndimension;

    // trivially extended in higher dims, with locality guaranteeing RNG state is local to node
    int lowerdims   = fine->_ndimension - coarse->_ndimension;
    assert(lowerdims >= 0);
    for(int d=0;d<lowerdims;d++){
      assert(fine->_simd_layout[d]==1);
      assert(fine->_processors[d]==1);
    }

    int multiplicity=1;
    for(int d=0;d<lowerdims;d++){
      multiplicity=multiplicity*fine->_rdimensions[d];
    }
    // local and global volumes subdivide cleanly after SIMDization
    for(int d=0;d<rngdims;d++){
      int fd= d+lowerdims;
      assert(coarse->_processors[d]  == fine->_processors[fd]);
      assert(coarse->_simd_layout[d] == fine->_simd_layout[fd]);
      assert(((fine->_rdimensions[fd] / coarse->_rdimensions[d])* coarse->_rdimensions[d])==fine->_rdimensions[fd]); 

      multiplicity = multiplicity *fine->_rdimensions[fd] / coarse->_rdimensions[d]; 
    }

    return multiplicity;
  }

  // real scalars are one component
  template<class scalar,class distribution,class generator> void fillScalar(scalar &s,distribution &dist,generator & gen)
  {
    s=dist(gen);
  }
  template<class distribution,class generator> void fillScalar(ComplexF &s,distribution &dist, generator &gen)
  {
    s=ComplexF(dist(gen),dist(gen));
  }
  template<class distribution,class generator> void fillScalar(ComplexD &s,distribution &dist,generator &gen)
  {
    s=ComplexD(dist(gen),dist(gen));
  }
  
  class GridRNGbase {
  public:
    // One generator per site.
    // Uniform and Gaussian distributions from these generators.
#ifdef RNG_RANLUX
    typedef std::ranlux48 RngEngine;
    typedef uint64_t      RngStateType;
    static const int RngStateCount = 15;
#endif 
#ifdef RNG_MT19937 
    typedef std::mt19937 RngEngine;
    typedef uint32_t     RngStateType;
    static const int     RngStateCount = std::mt19937::state_size;
#endif
#ifdef RNG_SITMO
    typedef sitmo::prng_engine 	RngEngine;
    typedef uint64_t    	RngStateType;
    static const int    	RngStateCount = 4;
#endif
    ///////////////////////
    // support for parallel init
    ///////////////////////
#ifdef RNG_SITMO
    static void Skip(RngEngine &eng)
    {
      uint64_t skip = 0x1; skip = skip<<40;
      eng.discard(skip);
    } 
#endif
    static RngEngine Reseed(RngEngine &eng)
    {
      const int reseeds=4;
      std::uniform_int_distribution<uint32_t> uid;
      std::vector<uint32_t> newseed(reseeds);
      for(int i=0;i<reseeds;i++){
	newseed[i] = uid(eng);
      }
      std::seed_seq sseq(newseed.begin(),newseed.end());
      return RngEngine(sseq);
    }    

    std::vector<RngEngine>                             _generators;
    std::vector<std::uniform_real_distribution<RealD> > _uniform;
    std::vector<std::normal_distribution<RealD> >       _gaussian;
    std::vector<std::discrete_distribution<int32_t> >   _bernoulli;
    std::vector<std::uniform_int_distribution<uint32_t> > _uid;

    void GetState(std::vector<RngStateType> & saved,RngEngine &eng) {
      saved.resize(RngStateCount);
      std::stringstream ss;
      ss<<eng;
      ss.seekg(0,ss.beg);
      for(int i=0;i<RngStateCount;i++){
	ss>>saved[i];
      }
    }
    void GetState(std::vector<RngStateType> & saved,int gen) {
      GetState(saved,_generators[gen]);
    }
    void SetState(std::vector<RngStateType> & saved,RngEngine &eng){
      assert(saved.size()==RngStateCount);
      std::stringstream ss;
      for(int i=0;i<RngStateCount;i++){
	ss<< saved[i]<<" ";
      }
      ss.seekg(0,ss.beg);
      ss>>eng;
    }
    void SetState(std::vector<RngStateType> & saved,int gen){
      SetState(saved,_generators[gen]);
    }
    void SetEngine(RngEngine &Eng, int gen){
      _generators[gen]=Eng;
    }
    void GetEngine(RngEngine &Eng, int gen){
      Eng=_generators[gen];
    }

    template<class source> void Seed(source &src, int gen)
    {
      _generators[gen] = RngEngine(src);
    }    
  };

  class GridSerialRNG : public GridRNGbase {
  public:

    GridSerialRNG() : GridRNGbase() {
      _generators.resize(1);
      _uniform.resize(1,std::uniform_real_distribution<RealD>{0,1});
      _gaussian.resize(1,std::normal_distribution<RealD>(0.0,1.0) );
      _bernoulli.resize(1,std::discrete_distribution<int32_t>{1,1});
      _uid.resize(1,std::uniform_int_distribution<uint32_t>() );
    }

    template <class sobj,class distribution> inline void fill(sobj &l,std::vector<distribution> &dist){

      typedef typename sobj::scalar_type scalar_type;
 
      int words = sizeof(sobj)/sizeof(scalar_type);

      scalar_type *buf = (scalar_type *) & l;

      dist[0].reset();
      for(int idx=0;idx<words;idx++){
	fillScalar(buf[idx],dist[0],_generators[0]);
      }
      
      CartesianCommunicator::BroadcastWorld(0,(void *)&l,sizeof(l));

    };
    template <class distribution>  inline void fill(ComplexF &l,std::vector<distribution> &dist){
      dist[0].reset();
      fillScalar(l,dist[0],_generators[0]);
      CartesianCommunicator::BroadcastWorld(0,(void *)&l,sizeof(l));
    }
    template <class distribution>  inline void fill(ComplexD &l,std::vector<distribution> &dist){
      dist[0].reset();
      fillScalar(l,dist[0],_generators[0]);
      CartesianCommunicator::BroadcastWorld(0,(void *)&l,sizeof(l));
    }
    template <class distribution>  inline void fill(RealF &l,std::vector<distribution> &dist){
      dist[0].reset();
      fillScalar(l,dist[0],_generators[0]);
      CartesianCommunicator::BroadcastWorld(0,(void *)&l,sizeof(l));
    }
    template <class distribution>  inline void fill(RealD &l,std::vector<distribution> &dist){
      dist[0].reset();
      fillScalar(l,dist[0],_generators[0]);
      CartesianCommunicator::BroadcastWorld(0,(void *)&l,sizeof(l));
    }
    // vector fill
    template <class distribution>  inline void fill(vComplexF &l,std::vector<distribution> &dist){
      RealF *pointer=(RealF *)&l;
      dist[0].reset();
      for(int i=0;i<2*vComplexF::Nsimd();i++){
	fillScalar(pointer[i],dist[0],_generators[0]);
      }
      CartesianCommunicator::BroadcastWorld(0,(void *)&l,sizeof(l));
    }
    template <class distribution>  inline void fill(vComplexD &l,std::vector<distribution> &dist){
      RealD *pointer=(RealD *)&l;
      dist[0].reset();
      for(int i=0;i<2*vComplexD::Nsimd();i++){
	fillScalar(pointer[i],dist[0],_generators[0]);
      }
      CartesianCommunicator::BroadcastWorld(0,(void *)&l,sizeof(l));
    }
    template <class distribution>  inline void fill(vRealF &l,std::vector<distribution> &dist){
      RealF *pointer=(RealF *)&l;
      dist[0].reset();
      for(int i=0;i<vRealF::Nsimd();i++){
	fillScalar(pointer[i],dist[0],_generators[0]);
      }
      CartesianCommunicator::BroadcastWorld(0,(void *)&l,sizeof(l));
    }
    template <class distribution>  inline void fill(vRealD &l,std::vector<distribution> &dist){
      RealD *pointer=(RealD *)&l;
      dist[0].reset();
      for(int i=0;i<vRealD::Nsimd();i++){
	fillScalar(pointer[i],dist[0],_generators[0]);
      }
      CartesianCommunicator::BroadcastWorld(0,(void *)&l,sizeof(l));
    }

    void SeedFixedIntegers(const std::vector<int> &seeds){
      CartesianCommunicator::BroadcastWorld(0,(void *)&seeds[0],sizeof(int)*seeds.size());
      std::seed_seq src(seeds.begin(),seeds.end());
      Seed(src,0);
    }

  };

  class GridParallelRNG : public GridRNGbase {
  public:

    GridBase *_grid;
    int _vol;

    int generator_idx(int os,int is){
      return is*_grid->oSites()+os;
    }

    GridParallelRNG(GridBase *grid) : GridRNGbase() {
      _grid=grid;
      _vol =_grid->iSites()*_grid->oSites();

      _generators.resize(_vol);
      _uniform.resize(_vol,std::uniform_real_distribution<RealD>{0,1});
      _gaussian.resize(_vol,std::normal_distribution<RealD>(0.0,1.0) );
      _bernoulli.resize(_vol,std::discrete_distribution<int32_t>{1,1});
      _uid.resize(_vol,std::uniform_int_distribution<uint32_t>() );
    }

    template <class vobj,class distribution> inline void fill(Lattice<vobj> &l,std::vector<distribution> &dist){

      typedef typename vobj::scalar_object scalar_object;
      typedef typename vobj::scalar_type scalar_type;
      typedef typename vobj::vector_type vector_type;
      
      int multiplicity = RNGfillable(_grid,l._grid);

      int     Nsimd =_grid->Nsimd();
      int     osites=_grid->oSites();
      int words=sizeof(scalar_object)/sizeof(scalar_type);


      parallel_for(int ss=0;ss<osites;ss++){

	std::vector<scalar_object> buf(Nsimd);
	for(int m=0;m<multiplicity;m++) {// Draw from same generator multiplicity times

	  int sm=multiplicity*ss+m;      // Maps the generator site to the fine site

	  for(int si=0;si<Nsimd;si++){
	    int gdx = generator_idx(ss,si); // index of generator state
	    scalar_type *pointer = (scalar_type *)&buf[si];
	    dist[gdx].reset();
	    for(int idx=0;idx<words;idx++){
	      fillScalar(pointer[idx],dist[gdx],_generators[gdx]);
	    }
	  }

	  // merge into SIMD lanes
	  merge(l._odata[sm],buf);
	}
      }
    };

    void SeedFixedIntegers(const std::vector<int> &seeds){

      CartesianCommunicator::BroadcastWorld(0,(void *)&seeds[0],sizeof(int)*seeds.size());

      std::seed_seq source(seeds.begin(),seeds.end());

      RngEngine master_engine(source);

#ifdef RNG_SITMO
      std::vector<int> gcoor;

      for(int gidx=0;gidx<_grid->_gsites;gidx++){

	Skip(master_engine); // advance the state; does this work?

	int rank,o_idx,i_idx;
	_grid->GlobalIndexToGlobalCoor(gidx,gcoor);
	_grid->GlobalCoorToRankIndex(rank,o_idx,i_idx,gcoor);

	if( rank == _grid->ThisRank() ){
	  int l_idx=generator_idx(o_idx,i_idx);
	  _generators[l_idx] = master_engine;
	}

      }
#else 
      // Machine and thread decomposition dependent seeding
      // is efficient and maximally parallel; but not
      // reproducible from machine to machine. Not ideal, but fast.
      // Different seed for each node.
      {
	int Nproc = _grid->ProcessorCount();
	int me= _grid->ThisRank();
	std::vector<RngEngine> seeders(Nproc);
	
	for(int p=0;p<Nproc;p++){
	  seeders[p] = Reseed(master_engine);
	}
	master_engine = seeders[me];
      }

      // Different seed for each thread
      {
	int Nthread = GridThread::GetThreads();
	std::vector<RngEngine> seeders(Nthread);
	for(int t=0;t<Nthread;t++){
	  seeders[t] = Reseed(master_engine);
	}

	parallel_for(int t=0;t<Nthread;t++) {
	  master_engine = seeders[t];
	  for(int l=0;l<_grid->lSites();l++) {
	    if ( (l%Nthread)==t ) {
	      _generators[l] = Reseed(master_engine);
	    }
	  }
	}
      }
#endif
    }

    uint32_t GlobalU01(int gsite){

      std::vector<int> gcoor;
      _grid->GlobalIndexToGlobalCoor(gsite,gcoor);

      int rank,o_idx,i_idx;
      _grid->GlobalCoorToRankIndex(rank,o_idx,i_idx,gcoor);

      int l_idx=generator_idx(o_idx,i_idx);

      uint32_t the_number;
      if( rank == _grid->ThisRank() ){
	the_number = _uid[l_idx](_generators[l_idx]);
      }

      _grid->Broadcast(rank,(void *)&the_number,sizeof(the_number));

      return the_number;
    }


  };

  template <class vobj> inline void random(GridParallelRNG &rng,Lattice<vobj> &l){
    rng.fill(l,rng._uniform);
  }

  template <class vobj> inline void gaussian(GridParallelRNG &rng,Lattice<vobj> &l){
    rng.fill(l,rng._gaussian);
  }
  
  template <class vobj> inline void bernoulli(GridParallelRNG &rng,Lattice<vobj> &l){
    rng.fill(l,rng._bernoulli);
  }

  template <class sobj> inline void random(GridSerialRNG &rng,sobj &l){
    rng.fill(l,rng._uniform);
  }
  
  template <class sobj> inline void gaussian(GridSerialRNG &rng,sobj &l){
    rng.fill(l,rng._gaussian);
  }
  
  template <class sobj> inline void bernoulli(GridSerialRNG &rng,sobj &l){
    rng.fill(l,rng._bernoulli);
  }

}
#endif
