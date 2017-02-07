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

  // Wrap seed_seq to give common interface with random_device
  // Should rather wrap random_device and have a generate
  class fixedSeed {
  public:

    typedef std::seed_seq::result_type result_type;

    std::seed_seq src;
    
    template<class int_type> fixedSeed(const std::vector<int_type> &seeds) : src(seeds.begin(),seeds.end()) {};

    template< class RandomIt > void generate( RandomIt begin, RandomIt end ) {
      src.generate(begin,end);
    }

  };


  class deviceSeed {
  public:

    std::random_device rd;

    typedef std::random_device::result_type result_type;
    
    deviceSeed(void) : rd(){};

    template< class RandomIt > void generate( RandomIt begin, RandomIt end ) {
      for(RandomIt it=begin; it!=end;it++){
	*it = rd();
      }
    }
  };

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

    int _seeded;
    // One generator per site.
    // Uniform and Gaussian distributions from these generators.
#ifdef RNG_RANLUX
    typedef uint64_t      RngStateType;
    typedef std::ranlux48 RngEngine;
    static const int RngStateCount = 15;
#else
    typedef std::mt19937 RngEngine;
    typedef uint32_t     RngStateType;
    static const int     RngStateCount = std::mt19937::state_size;
#endif
    std::vector<RngEngine>                             _generators;
    std::vector<std::uniform_real_distribution<RealD>> _uniform;
    std::vector<std::normal_distribution<RealD>>       _gaussian;
    std::vector<std::discrete_distribution<int32_t>>   _bernoulli;

    void GetState(std::vector<RngStateType> & saved,int gen) {
      saved.resize(RngStateCount);
      std::stringstream ss;
      ss<<_generators[gen];
      ss.seekg(0,ss.beg);
      for(int i=0;i<RngStateCount;i++){
	ss>>saved[i];
      }
    }
    void SetState(std::vector<RngStateType> & saved,int gen){
      assert(saved.size()==RngStateCount);
      std::stringstream ss;
      for(int i=0;i<RngStateCount;i++){
	ss<< saved[i]<<" ";
      }
      ss.seekg(0,ss.beg);
      ss>>_generators[gen];
    }
  };

  class GridSerialRNG : public GridRNGbase {
  public:

    // FIXME ... do we require lockstep draws of randoms 
    // from all nodes keeping seeds consistent.
    // place a barrier/broadcast in the fill routine

    GridSerialRNG() : GridRNGbase() {
      _generators.resize(1);
      _uniform.resize(1,std::uniform_real_distribution<RealD>{0,1});
      _gaussian.resize(1,std::normal_distribution<RealD>(0.0,1.0) );
      _bernoulli.resize(1,std::discrete_distribution<int32_t>{1,1});
      _seeded=0;
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

    template<class source> void Seed(source &src)
    {
      _generators[0] = RngEngine(src);
      _seeded=1;
    }    
    void SeedRandomDevice(void){
      deviceSeed src;
      Seed(src);
    }
    void SeedFixedIntegers(const std::vector<int> &seeds){
      CartesianCommunicator::BroadcastWorld(0,(void *)&seeds[0],sizeof(int)*seeds.size());
      fixedSeed src(seeds);
      Seed(src);
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
      _seeded=0;
    }



    //FIXME implement generic IO and create state save/restore
    //void SaveState(const std::string<char> &file);
    //void LoadState(const std::string<char> &file);

    template <class vobj,class distribution> inline void fill(Lattice<vobj> &l,std::vector<distribution> &dist){

      typedef typename vobj::scalar_object scalar_object;
      typedef typename vobj::scalar_type scalar_type;
      typedef typename vobj::vector_type vector_type;
      
      int multiplicity = RNGfillable(_grid,l._grid);

      int     Nsimd =_grid->Nsimd();
      int     osites=_grid->oSites();
      int words=sizeof(scalar_object)/sizeof(scalar_type);


PARALLEL_FOR_LOOP
      for(int ss=0;ss<osites;ss++){

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

    // This loop could be made faster to avoid the Ahmdahl by
    // i)  seed generators on each timeslice, for x=y=z=0;
    // ii) seed generators on each z for x=y=0
    // iii)seed generators on each y,z for x=0
    // iv) seed generators on each y,z,x 
    // made possible by physical indexing.
    template<class source> void Seed(source &src)
    {

      typedef typename source::result_type seed_t;
      std::uniform_int_distribution<seed_t> uid;

      int numseed=4;
      int gsites = _grid->_gsites;
      std::vector<seed_t> site_init(numseed);
      std::vector<int> gcoor;


      // Master RngEngine
      std::vector<seed_t> master_init(numseed);  src.generate(master_init.begin(),master_init.end());
      _grid->Broadcast(0,(void *)&master_init[0],sizeof(seed_t)*numseed);
      fixedSeed master_seed(master_init);
      RngEngine master_engine(master_seed);

      // Per node RngEngine
      std::vector<seed_t> node_init(numseed);
      for(int r=0;r<_grid->ProcessorCount();r++) {

	std::vector<seed_t> rank_init(numseed);
	for(int i=0;i<numseed;i++) rank_init[i] = uid(master_engine);

	std::cout << GridLogMessage << "SeedSeq for rank "<<r;
	for(int i=0;i<numseed;i++) std::cout<<" "<<rank_init[i];
	std::cout <<std::endl;

	if ( r==_grid->ThisRank() ) { 
	  for(int i=0;i<numseed;i++) node_init[i] = rank_init[i];
	}

      }

      ////////////////////////////////////////////////////
      // Set up a seed_seq wrapper with these 8 words
      // and draw for each site within node.
      ////////////////////////////////////////////////////
      fixedSeed node_seed(node_init);
      RngEngine node_engine(node_seed);

      for(int gidx=0;gidx<gsites;gidx++){
	int rank,o_idx,i_idx;

	_grid->GlobalIndexToGlobalCoor(gidx,gcoor);
	_grid->GlobalCoorToRankIndex(rank,o_idx,i_idx,gcoor);

	if( rank == _grid->ThisRank() ){
	  int l_idx=generator_idx(o_idx,i_idx);
	  for(int i=0;i<numseed;i++)  site_init[i] = uid(node_engine);
	  fixedSeed site_seed(site_init);
	  _generators[l_idx] = RngEngine(site_seed);
	}
      }
      _seeded=1;
    }    
    void SeedRandomDevice(void){
      deviceSeed src;
      Seed(src);
    }
    void SeedFixedIntegers(const std::vector<int> &seeds){
      CartesianCommunicator::BroadcastWorld(0,(void *)&seeds[0],sizeof(int)*seeds.size());
      fixedSeed src(seeds);
      Seed(src);
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
