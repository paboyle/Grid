/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/lattice/Lattice_rng.h

    Copyright (C) 2015

    Author: Peter Boyle <paboyle@ph.ed.ac.uk>
    Author: Guido Cossu <guido.cossu@ed.ac.uk>

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

#ifdef RNG_SITMO
#include <Grid/sitmo_rng/sitmo_prng_engine.hpp>
#endif 

#if defined(RNG_SITMO)
#define RNG_FAST_DISCARD
#else 
#undef  RNG_FAST_DISCARD
#endif

NAMESPACE_BEGIN(Grid);

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

  
// merge of April 11 2017
// this function is necessary for the LS vectorised field
inline int RNGfillable_general(GridBase *coarse,GridBase *fine)
{
  int rngdims = coarse->_ndimension;
    
  // trivially extended in higher dims, with locality guaranteeing RNG state is local to node
  int lowerdims   = fine->_ndimension - coarse->_ndimension;  assert(lowerdims >= 0);
  // assumes that the higher dimensions are not using more processors
  // all further divisions are local
  for(int d=0;d<lowerdims;d++) assert(fine->_processors[d]==1);
  for(int d=0;d<rngdims;d++) assert(coarse->_processors[d] == fine->_processors[d+lowerdims]);

  // then divide the number of local sites
  // check that the total number of sims agree, meanse the iSites are the same
  assert(fine->Nsimd() == coarse->Nsimd());

  // check that the two grids divide cleanly
  assert( (fine->lSites() / coarse->lSites() ) * coarse->lSites() == fine->lSites() );

  return fine->lSites() / coarse->lSites();
}
  
// real scalars are one component
template<class scalar,class distribution,class generator> 
void fillScalar(scalar &s,distribution &dist,generator & gen)
{
  s=dist(gen);
}
template<class distribution,class generator> 
void fillScalar(ComplexF &s,distribution &dist, generator &gen)
{
  //  s=ComplexF(dist(gen),dist(gen));
  s.real(dist(gen));
  s.imag(dist(gen));
}
template<class distribution,class generator> 
void fillScalar(ComplexD &s,distribution &dist,generator &gen)
{
  //  s=ComplexD(dist(gen),dist(gen));
  s.real(dist(gen));
  s.imag(dist(gen));
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
  static const int    	RngStateCount = 13;
#endif

  std::vector<RngEngine>                             _generators;
  std::vector<std::uniform_real_distribution<RealD> > _uniform;
  std::vector<std::normal_distribution<RealD> >       _gaussian;
  std::vector<std::discrete_distribution<int32_t> >   _bernoulli;
  std::vector<std::uniform_int_distribution<uint32_t> > _uid;

  ///////////////////////
  // support for parallel init
  ///////////////////////
#ifdef RNG_FAST_DISCARD
  static void Skip(RngEngine &eng,uint64_t site)
  {
#if 0
    /////////////////////////////////////////////////////////////////////////////////////
    // Skip by 2^40 elements between successive lattice sites
    // This goes by 10^12.
    // Consider quenched updating; likely never exceeding rate of 1000 sweeps
    // per second on any machine. This gives us of order 10^9 seconds, or 100 years
    // skip ahead.
    // For HMC unlikely to go at faster than a solve per second, and 
    // tens of seconds per trajectory so this is clean in all reasonable cases,
    // and margin of safety is orders of magnitude.
    // We could hack Sitmo to skip in the higher order words of state if necessary
    //
    // Replace with 2^30 ; avoid problem on large volumes
    //
    /////////////////////////////////////////////////////////////////////////////////////
    //      uint64_t skip = site+1;  //   Old init Skipped then drew.  Checked compat with faster init
    const int shift = 30;

    ////////////////////////////////////////////////////////////////////
    // Weird compiler bug in Intel 2018.1 under O3 was generating 32bit and not 64 bit left shift.
    ////////////////////////////////////////////////////////////////////
    volatile uint64_t skip = site;

    skip = skip<<shift;

    assert((skip >> shift)==site); // check for overflow

    eng.discard(skip);
#else
    eng.discardhi(site);
#endif
    //      std::cout << " Engine  " <<site << " state " <<eng<<std::endl;
  } 
#endif
  static RngEngine Reseed(RngEngine &eng)
  {
    std::vector<uint32_t> newseed;
    std::uniform_int_distribution<uint32_t> uid;
    return Reseed(eng,newseed,uid);
  }
  static RngEngine Reseed(RngEngine &eng,std::vector<uint32_t> & newseed,
			  std::uniform_int_distribution<uint32_t> &uid)
  {
    const int reseeds=4;
      
    newseed.resize(reseeds);
    for(int i=0;i<reseeds;i++){
      newseed[i] = uid(eng);
    }
    std::seed_seq sseq(newseed.begin(),newseed.end());
    return RngEngine(sseq);
  }    

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

  }

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

    void SeedUniqueString(const std::string &s){
      std::vector<int> seeds;
      std::stringstream sha;
      seeds = GridChecksum::sha256_seeds(s);
      for(int i=0;i<seeds.size();i++) { 
        sha << std::hex << seeds[i];
      }
      std::cout << GridLogMessage << "Intialising serial RNG with unique string '" 
                << s << "'" << std::endl;
      std::cout << GridLogMessage << "Seed SHA256: " << sha.str() << std::endl;
      SeedFixedIntegers(seeds);
    }
};

class GridParallelRNG : public GridRNGbase {
private:
  double _time_counter;
  GridBase *_grid;
  unsigned int _vol;

public:
  GridBase *Grid(void) const { return _grid; }
  int generator_idx(int os,int is) {
    return is*_grid->oSites()+os;
  }

  GridParallelRNG(GridBase *grid) : GridRNGbase() {
    _grid = grid;
    _vol  =_grid->iSites()*_grid->oSites();

    _generators.resize(_vol);
    _uniform.resize(_vol,std::uniform_real_distribution<RealD>{0,1});
    _gaussian.resize(_vol,std::normal_distribution<RealD>(0.0,1.0) );
    _bernoulli.resize(_vol,std::discrete_distribution<int32_t>{1,1});
    _uid.resize(_vol,std::uniform_int_distribution<uint32_t>() );
  }
  template <class vobj,class distribution> inline void fill(Lattice<vobj> &l,std::vector<distribution> &dist)
  {
    if ( l.Grid()->_isCheckerBoarded ) {
      Lattice<vobj> tmp(_grid);
      fill(tmp,dist);
      pickCheckerboard(l.Checkerboard(),l,tmp);
      return;
    }
    typedef typename vobj::scalar_object scalar_object;
    typedef typename vobj::scalar_type scalar_type;
    typedef typename vobj::vector_type vector_type;

    double inner_time_counter = usecond();

    int multiplicity = RNGfillable_general(_grid, l.Grid()); // l has finer or same grid
    int Nsimd  = _grid->Nsimd();  // guaranteed to be the same for l.Grid() too
    int osites = _grid->oSites();  // guaranteed to be <= l.Grid()->oSites() by a factor multiplicity
    int words  = sizeof(scalar_object) / sizeof(scalar_type);

    autoView(l_v, l, CpuWrite);
    thread_for( ss, osites, {
      ExtractBuffer<scalar_object> buf(Nsimd);
      for (int m = 0; m < multiplicity; m++) {  // Draw from same generator multiplicity times

	int sm = multiplicity * ss + m;  // Maps the generator site to the fine site

	for (int si = 0; si < Nsimd; si++) {
            
	  int gdx = generator_idx(ss, si);  // index of generator state
	  scalar_type *pointer = (scalar_type *)&buf[si];
	  dist[gdx].reset();
	  for (int idx = 0; idx < words; idx++) 
	    fillScalar(pointer[idx], dist[gdx], _generators[gdx]);
	}
	// merge into SIMD lanes, FIXME suboptimal implementation
	merge(l_v[sm], buf);
      }
      });
    //    });

    _time_counter += usecond()- inner_time_counter;
  }

    void SeedUniqueString(const std::string &s){
      std::vector<int> seeds;
      seeds = GridChecksum::sha256_seeds(s);
      std::cout << GridLogMessage << "Intialising parallel RNG with unique string '" 
                << s << "'" << std::endl;
      std::cout << GridLogMessage << "Seed SHA256: " << GridChecksum::sha256_string(seeds) << std::endl;
      SeedFixedIntegers(seeds);
    }
  void SeedFixedIntegers(const std::vector<int> &seeds, int britney=0){

    // Everyone generates the same seed_seq based on input seeds
    CartesianCommunicator::BroadcastWorld(0,(void *)&seeds[0],sizeof(int)*seeds.size());

    std::seed_seq source(seeds.begin(),seeds.end());

    RngEngine master_engine(source);

#ifdef RNG_FAST_DISCARD
    ////////////////////////////////////////////////
    // Skip ahead through a single stream.
    // Applicable to SITMO and other has based/crypto RNGs
    // Should be applicable to Mersenne Twister, but the C++11
    // MT implementation does not implement fast discard even though
    // in principle this is possible
    ////////////////////////////////////////////////
    thread_for( lidx, _grid->lSites(), {

	int64_t gidx;
	int o_idx;
	int i_idx;
	int rank;
	Coordinate pcoor;
	Coordinate lcoor;
	Coordinate gcoor;
	_grid->LocalIndexToLocalCoor(lidx,lcoor);
	pcoor=_grid->ThisProcessorCoor();
	_grid->ProcessorCoorLocalCoorToGlobalCoor(pcoor,lcoor,gcoor);
	_grid->GlobalCoorToGlobalIndex(gcoor,gidx);

	_grid->GlobalCoorToRankIndex(rank,o_idx,i_idx,gcoor);

	assert(rank == _grid->ThisRank() );
	
	int l_idx=generator_idx(o_idx,i_idx);
	_generators[l_idx] = master_engine;
	if ( britney ) { 
	  Skip(_generators[l_idx],l_idx); // Skip to next RNG sequence
	} else { 	
	  Skip(_generators[l_idx],gidx); // Skip to next RNG sequence
	}
    });
#else 
    ////////////////////////////////////////////////////////////////
    // Machine and thread decomposition dependent seeding is efficient
    // and maximally parallel; but NOT reproducible from machine to machine. 
    // Not ideal, but fastest way to reseed all nodes.
    ////////////////////////////////////////////////////////////////
    {
      // Obtain one Reseed per processor
      int Nproc = _grid->ProcessorCount();
      std::vector<RngEngine> seeders(Nproc);
      int me= _grid->ThisRank();
      for(int p=0;p<Nproc;p++){
	seeders[p] = Reseed(master_engine);
      }
      master_engine = seeders[me];
    }

    {
      // Obtain one reseeded generator per thread      
      int Nthread = 32; // Hardwire a good level or parallelism
      std::vector<RngEngine> seeders(Nthread);
      for(int t=0;t<Nthread;t++){
	seeders[t] = Reseed(master_engine);
      }

      thread_for( t, Nthread, {
	// set up one per local site in threaded fashion
	std::vector<uint32_t> newseeds;
	std::uniform_int_distribution<uint32_t> uid;	
	for(int l=0;l<_grid->lSites();l++) {
	  if ( (l%Nthread)==t ) {
	    _generators[l] = Reseed(seeders[t],newseeds,uid);
	  }
	}
      });
    }
#endif
  }

  void Report(){
    std::cout << GridLogMessage << "Time spent in the fill() routine by GridParallelRNG: "<< _time_counter/1e3 << " ms" << std::endl;
  }


  ////////////////////////////////////////////////////////////////////////
  // Support for rigorous test of RNG's
  // Return uniform random uint32_t from requested site generator
  ////////////////////////////////////////////////////////////////////////
  uint32_t GlobalU01(int gsite){

    uint32_t the_number;
    // who
    int rank,o_idx,i_idx;
    Coordinate gcoor;
    _grid->GlobalIndexToGlobalCoor(gsite,gcoor);
    _grid->GlobalCoorToRankIndex(rank,o_idx,i_idx,gcoor);

    // draw
    int l_idx=generator_idx(o_idx,i_idx);
    if( rank == _grid->ThisRank() ){
      the_number = _uid[l_idx](_generators[l_idx]);
    }
      
    // share & return
    _grid->Broadcast(rank,(void *)&the_number,sizeof(the_number));
    return the_number;
  }

};

template <class vobj> inline void random(GridParallelRNG &rng,Lattice<vobj> &l)   { rng.fill(l,rng._uniform);  }
template <class vobj> inline void gaussian(GridParallelRNG &rng,Lattice<vobj> &l) { rng.fill(l,rng._gaussian); }
template <class vobj> inline void bernoulli(GridParallelRNG &rng,Lattice<vobj> &l){ rng.fill(l,rng._bernoulli);}

template <class sobj> inline void random(GridSerialRNG &rng,sobj &l)   { rng.fill(l,rng._uniform  ); }
template <class sobj> inline void gaussian(GridSerialRNG &rng,sobj &l) { rng.fill(l,rng._gaussian ); }
template <class sobj> inline void bernoulli(GridSerialRNG &rng,sobj &l){ rng.fill(l,rng._bernoulli); }

NAMESPACE_END(Grid);
#endif
