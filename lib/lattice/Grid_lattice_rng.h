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
	buf[idx] = dist(_generators[0]);
      }
      
      CartesianCommunicator::BroadcastWorld(0,(void *)&l,sizeof(l));

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

      for(int gidx=0;gidx<_grid->_gsites;gidx++){
	int rank,o_idx,i_idx;
	_grid->GlobalIndexToGlobalCoor(gidx,gcoor);
	_grid->GlobalCoorToRankIndex(rank,o_idx,i_idx,gcoor);

	int l_idx=generator_idx(o_idx,i_idx);

	typename source::result_type init = src();

	_grid->Broadcast(0,(void *)&init,sizeof(init));
	if( rank == _grid->ThisRank() ){
	  _generators[l_idx] = std::ranlux48(init);
	}
      }
      _seeded=1;
    }    

    //FIXME implement generic IO and create state save/restore
    //void SaveState(const std::string<char> &file);
    //void LoadState(const std::string<char> &file);

    template <class vobj,class distribution> inline void fill(Lattice<vobj> &l,distribution &dist){

      typedef typename vobj::scalar_type scalar_type;
      typedef typename vobj::vector_type vector_type;
      
      conformable(_grid,l._grid);

      int     Nsimd =_grid->Nsimd();
      int     osites=_grid->oSites();

      int words = sizeof(vobj)/sizeof(vector_type);
      std::vector<std::vector<scalar_type> > buf(Nsimd,std::vector<scalar_type>(words));
      std::vector<scalar_type *> pointers(Nsimd);  
      
      for(int ss=0;ss<osites;ss++){

	for(int si=0;si<Nsimd;si++){

	  int gdx = generator_idx(ss,si); // index of generator state

	  pointers[si] = (scalar_type *)&buf[si][0];
	  for(int idx=0;idx<words;idx++){
	    pointers[si][idx] = dist(_generators[gdx]);
	  }

	}
	// merge into SIMD lanes
	merge(l._odata[ss],pointers);
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
