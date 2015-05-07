#ifndef GRID_LATTICE_RNG_H
#define GRID_LATTICE_RNG_H

#include <random>

namespace Grid {

  // Wrap seed_seq to give common interface with random_device
  class fixedSeed {
  public:

    typedef std::seed_seq::result_type result_type;

    std::seed_seq src;
    
    fixedSeed(std::vector<int> &seeds) : src(seeds.begin(),seeds.end()) {};

    result_type operator () (void){

      std::vector<result_type> list(1);

      src.generate(list.begin(),list.end());

      return list[0];

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

   GridRNGbase() : _uniform{0,1}, _gaussian(0.0,1.0) {};

    int _seeded;
    // One generator per site.
    // Uniform and Gaussian distributions from these generators.
    std::vector<std::ranlux48>             _generators;
    std::uniform_real_distribution<double> _uniform;
    std::normal_distribution<double>       _gaussian;


  };

  class GridSerialRNG : public GridRNGbase {
  public:

    // FIXME ... do we require lockstep draws of randoms 
    // from all nodes keeping seeds consistent.
    // place a barrier/broadcast in the fill routine
    template<class source> void Seed(source &src)
    {
      typename source::result_type init = src();
      CartesianCommunicator::BroadcastWorld(0,(void *)&init,sizeof(init));
      _generators[0] = std::ranlux48(init);
      _seeded=1;
    }    

    GridSerialRNG() : GridRNGbase() {
      _generators.resize(1);
      _seeded=0;
    }



    template <class sobj,class distribution> inline void fill(sobj &l,distribution &dist){

      typedef typename sobj::scalar_type scalar_type;
 
      int words = sizeof(sobj)/sizeof(scalar_type);

      scalar_type *buf = (scalar_type *) & l;

      for(int idx=0;idx<words;idx++){
	fillScalar(buf[idx],dist,_generators[0]);
      }
      
      CartesianCommunicator::BroadcastWorld(0,(void *)&l,sizeof(l));

    };

    template <class distribution>  inline void fill(ComplexF &l,distribution &dist){
      fillScalar(l,dist,_generators[0]);
      CartesianCommunicator::BroadcastWorld(0,(void *)&l,sizeof(l));
    }
    template <class distribution>  inline void fill(ComplexD &l,distribution &dist){
      fillScalar(l,dist,_generators[0]);
      CartesianCommunicator::BroadcastWorld(0,(void *)&l,sizeof(l));
    }
    template <class distribution>  inline void fill(RealF &l,distribution &dist){
      fillScalar(l,dist,_generators[0]);
      CartesianCommunicator::BroadcastWorld(0,(void *)&l,sizeof(l));
    }
    template <class distribution>  inline void fill(RealD &l,distribution &dist){
      fillScalar(l,dist,_generators[0]);
      CartesianCommunicator::BroadcastWorld(0,(void *)&l,sizeof(l));
    }
    // vector fill
    template <class distribution>  inline void fill(vComplexF &l,distribution &dist){
      RealF *pointer=(RealF *)&l;
      for(int i=0;i<2*vComplexF::Nsimd();i++){
	fillScalar(pointer[i],dist,_generators[0]);
      }
      CartesianCommunicator::BroadcastWorld(0,(void *)&l,sizeof(l));
    }
    template <class distribution>  inline void fill(vComplexD &l,distribution &dist){
      RealD *pointer=(RealD *)&l;
      for(int i=0;i<2*vComplexD::Nsimd();i++){
	fillScalar(pointer[i],dist,_generators[0]);
      }
      CartesianCommunicator::BroadcastWorld(0,(void *)&l,sizeof(l));
    }
    template <class distribution>  inline void fill(vRealF &l,distribution &dist){
      RealF *pointer=(RealF *)&l;
      for(int i=0;i<vRealF::Nsimd();i++){
	fillScalar(pointer[i],dist,_generators[0]);
      }
      CartesianCommunicator::BroadcastWorld(0,(void *)&l,sizeof(l));
    }
    template <class distribution>  inline void fill(vRealD &l,distribution &dist){
      RealD *pointer=(RealD *)&l;
      for(int i=0;i<vRealD::Nsimd();i++){
	fillScalar(pointer[i],dist,_generators[0]);
      }
      CartesianCommunicator::BroadcastWorld(0,(void *)&l,sizeof(l));
    }


    void SeedRandomDevice(void){
      std::random_device rd;
      Seed(rd);
    }
    void SeedFixedIntegers(std::vector<int> &seeds){
      fixedSeed src(seeds);
      Seed(src);
    }

  };

  class GridParallelRNG : public GridRNGbase {
  public:
    // One generator per site.
    std::vector<std::ranlux48>             _generators;
    
    // Uniform and Gaussian distributions from these generators.
    std::uniform_real_distribution<double> _uniform;
    std::normal_distribution<double>       _gaussian;

    GridBase *_grid;
    int _vol;

    int generator_idx(int os,int is){
      return is*_grid->oSites()+os;
    }

    GridParallelRNG(GridBase *grid) : GridRNGbase() {
      _grid=grid;
      _vol =_grid->iSites()*_grid->oSites();
      _generators.resize(_vol);
      _seeded=0;
    }


    // This loop could be made faster to avoid the Ahmdahl by
    // i)  seed generators on each timeslice, for x=y=z=0;
    // ii) seed generators on each z for x=y=0
    // iii)seed generators on each y,z for x=0
    // iv) seed generators on each y,z,x 
    // made possible by physical indexing.
    template<class source> void Seed(source &src)
    {
      std::vector<int> gcoor;

      int gsites = _grid->_gsites;

      typename source::result_type init = src();
      std::ranlux48 pseeder(init);
      std::uniform_int_distribution<uint64_t> ui;

      for(int gidx=0;gidx<gsites;gidx++){

	int rank,o_idx,i_idx;
	_grid->GlobalIndexToGlobalCoor(gidx,gcoor);
	_grid->GlobalCoorToRankIndex(rank,o_idx,i_idx,gcoor);

	int l_idx=generator_idx(o_idx,i_idx);
	
	std::vector<int> site_seeds(4);
	for(int i=0;i<4;i++){
	  site_seeds[i]= ui(pseeder);
	}

	_grid->Broadcast(0,(void *)&site_seeds[0],sizeof(int)*site_seeds.size());

	if( rank == _grid->ThisRank() ){
	  fixedSeed ssrc(site_seeds);
	  typename source::result_type sinit = ssrc();
	  _generators[l_idx] = std::ranlux48(sinit);
	}
      }
      _seeded=1;
    }    

    //FIXME implement generic IO and create state save/restore
    //void SaveState(const std::string<char> &file);
    //void LoadState(const std::string<char> &file);

    template <class vobj,class distribution> inline void fill(Lattice<vobj> &l,distribution &dist){

      typedef typename vobj::scalar_object scalar_object;
      typedef typename vobj::scalar_type scalar_type;
      typedef typename vobj::vector_type vector_type;
      
      conformable(_grid,l._grid);

      int     Nsimd =_grid->Nsimd();
      int     osites=_grid->oSites();
      int words=sizeof(scalar_object)/sizeof(scalar_type);

      std::vector<scalar_object> buf(Nsimd);
      
      for(int ss=0;ss<osites;ss++){
	for(int si=0;si<Nsimd;si++){

	  int gdx = generator_idx(ss,si); // index of generator state
	  scalar_type *pointer = (scalar_type *)&buf[si];
	  for(int idx=0;idx<words;idx++){
	    fillScalar(pointer[idx],dist,_generators[gdx]);
	  }

	}
	// merge into SIMD lanes
	merge(l._odata[ss],buf);
      }
    };

    void SeedRandomDevice(void){
      std::random_device rd;
      Seed(rd);
    }
    void SeedFixedIntegers(std::vector<int> &seeds){
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


  template <class sobj> inline void random(GridSerialRNG &rng,sobj &l){
    rng.fill(l,rng._uniform);
  }
  template <class sobj> inline void gaussian(GridSerialRNG &rng,sobj &l){
    rng.fill(l,rng._gaussian);
  }

}
#endif
