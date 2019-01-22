#ifndef Hadrons_MDistil_DistilVectors_hpp_
#define Hadrons_MDistil_DistilVectors_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Solver.hpp>
#include <Hadrons/EigenPack.hpp>
#include <Hadrons/A2AVectors.hpp>
#include <Hadrons/DilutedNoise.hpp>

BEGIN_HADRONS_NAMESPACE

template<typename LatticeObj>
class Perambulator : Serializable{
	  // TODO: The next line makes friends across all combinations
	  //     (not much of a problem given all public anyway ...)
	  //      FYI, the bug here was that I forgot that the friend is templated
  template<typename T> friend std::ostream & operator<<(std::ostream &os, const Perambulator<T>& p);
protected:
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS( Perambulator,
                                  std::string,             ID,          // Allows owner to specialise
                                  std::string,             Provenance,  // For info only
                                  std::vector<int>,        dimensions,
                                  std::vector<LatticeObj>, perambulator,
                                  // Following items are redundant, but useful
                                  int,     nd,          // Number of dimensions
                                  size_t,  NumElements); // Number of elements
protected:
  // Constructor common code
  inline void ConstructCommon(const int * Dimensions) {
    assert(nd > 0);
    dimensions.resize(nd);
    NumElements = 1;
    for(int i = 0 ; i < nd ; i++) {
      assert(Dimensions[i] > 0);
      NumElements *= (size_t) Dimensions[i];
      dimensions[i] = Dimensions[i];
    }
    //const LatticeObj perambulatorDefault;
    perambulator.resize(NumElements);//,perambulatorDefault);
  }
public:
  // Constructor with dimensions passed as std::vector<int>
  inline Perambulator(const std::vector<int> & Dimensions)
  : nd {(int) Dimensions.size()} {
    ConstructCommon( &Dimensions[0] ); }
  
  // Constructor with dimensions passed as std::vector<int>
  inline Perambulator(const std::vector<int> & Dimensions, const std::string sID)
  : nd {(int) Dimensions.size()}, ID(sID) {
    ConstructCommon( &Dimensions[0] ); }
  
  // Constructor with dimensions passed as std::vector<int>
  inline Perambulator(const std::vector<int> & Dimensions, const std::string sID, const std::string sProvenance)
  : nd {(int) Dimensions.size()}, ID(sID), Provenance(sProvenance) {
    ConstructCommon( &Dimensions[0] ); }
  
  // Constructor with dimensions passed as individual parameters
  // FYI: The caller is free to ignore the names and use the indices however they see fit
  inline Perambulator(int NumNoise, int NumEvec=1, int NumTime=1, int NumSpin=1, int I_k=1, int I_t=1, int I_s=1) {
    int Dimensions[]={NumNoise,NumEvec,NumTime,NumSpin,I_k,I_t,I_s};
    nd = sizeof(Dimensions)/sizeof(Dimensions[0]);
    while( nd > 1 && Dimensions[nd-1] == 1 )
      nd--;
    ConstructCommon( Dimensions );
  }
  
  inline LatticeObj & operator()(size_t count, const int * Coord) {
    assert( count == nd );
    assert( Coord );
    size_t idx = 0;
    // C memory order (???)
    for( int d = 0 ; d < nd ; d++ ) {
      assert( Coord[d] < dimensions[d] );
      idx *= (size_t) dimensions[d];
      idx += (size_t) Coord[d];
    }
    return perambulator[idx];
  }

  inline LatticeObj & operator()(const std::vector<int> Coord) {
    return operator()(Coord.size(), &Coord[0]);
  }
  
  inline LatticeObj & operator()(int idxNoise, int idxEvec=0, int idxTime=0, int idxSpin=0, int I_k=0, int I_t=0, int I_s=0) {
    int MyIndex[]={idxNoise,idxEvec,idxTime,idxSpin,I_k,I_t,I_s};
    int i = sizeof(MyIndex)/sizeof(MyIndex[0]);
    assert( i >= nd );
    while( i > nd )
      assert(MyIndex[--i] == 0);
    return operator()(i, MyIndex);
  }
};


/******************************************************************************
 *                         DistilVectors                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MDistil)

class DistilVectorsPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(DistilVectorsPar,
		                    std::string, noise,
		                    std::string, perambulator,
		                    std::string, eigenPack,
                                    bool, multiFile);
};

template <typename FImpl>
class TDistilVectors: public Module<DistilVectorsPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TDistilVectors(const std::string name);
    // destructor
    virtual ~TDistilVectors(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(DistilVectors, TDistilVectors<FIMPL>, MDistil);

/******************************************************************************
 *                 TDistilVectors implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TDistilVectors<FImpl>::TDistilVectors(const std::string name)
: Module<DistilVectorsPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TDistilVectors<FImpl>::getInput(void)
{
    std::vector<std::string> in;

    in.push_back(par().noise);
    in.push_back(par().perambulator);
    in.push_back(par().eigenPack);

    return in;
}

template <typename FImpl>
std::vector<std::string> TDistilVectors<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName() + "_rho", getName() + "_phi"};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDistilVectors<FImpl>::setup(void)
{
   auto &noise = envGet(std::vector<std::vector<std::vector<SpinVector>>>, par().noise);
  
   envCreate(std::vector<FermionField>, getName() + "_rho", 1, 
		                    noise.size(), envGetGrid(FermionField));
   envCreate(std::vector<FermionField>, getName() + "_phi", 1, 
                 	            noise.size(), envGetGrid(FermionField)); 


  GridCartesian * grid4d = env().getGrid();
  std::vector<int> latt_size   = GridDefaultLatt();
  std::vector<int> simd_layout = GridDefaultSimd(Nd, vComplex::Nsimd());
  std::vector<int> mpi_layout  = GridDefaultMpi();
  std::vector<int> simd_layout_3 = GridDefaultSimd(Nd-1, vComplex::Nsimd());
  latt_size[Nd-1] = 1;
  simd_layout_3.push_back( 1 );
  mpi_layout[Nd-1] = 1;
  GridCartesian * grid3d = new GridCartesian(latt_size,simd_layout_3,mpi_layout,*grid4d);


  envTmp(LatticeSpinColourVector, "tmp2",1,LatticeSpinColourVector(grid4d));
  envTmp(LatticeColourVector, "tmp_nospin",1,LatticeColourVector(grid4d));
  envTmp(LatticeSpinColourVector, "tmp3d",1,LatticeSpinColourVector(grid3d));
  envTmp(LatticeColourVector, "tmp3d_nospin",1,LatticeColourVector(grid3d));
  envTmp(LatticeSpinColourVector, "sink_tslice",1,LatticeSpinColourVector(grid3d));
  envTmp(LatticeSpinColourVector, "sink4d",1,LatticeSpinColourVector(grid4d));
  envTmp(LatticeColourVector, "evec3d",1,LatticeColourVector(grid3d));
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDistilVectors<FImpl>::execute(void)
{
   
    auto        &noise     = envGet(std::vector<std::vector<std::vector<SpinVector>>>, par().noise);
    //auto        &perambulator   = envGet(std::vector<SpinVector>, par().perambulator);
    auto        &perambulator   = envGet(Perambulator<SpinVector>, par().noise);
    auto        &epack   = envGet(Grid::Hadrons::EigenPack<LatticeColourVector>, par().eigenPack);
    auto        &rho       = envGet(std::vector<FermionField>, getName() + "_rho");
    auto        &phi       = envGet(std::vector<FermionField>, getName() + "_phi");


  envGetTmp(LatticeSpinColourVector, tmp2);
  envGetTmp(LatticeColourVector, tmp_nospin);
  envGetTmp(LatticeSpinColourVector, tmp3d);
  envGetTmp(LatticeColourVector, tmp3d_nospin);
  envGetTmp(LatticeSpinColourVector, sink_tslice);
  envGetTmp(LatticeSpinColourVector, sink4d);
  envGetTmp(LatticeColourVector, evec3d);

  GridCartesian * grid4d = env().getGrid();

  int Ntlocal = grid4d->LocalDimensions()[3];
  int Ntfirst = grid4d->LocalStarts()[3];

  int tsrc=0;
  int nnoise=1;
  int LI=6;
  int Ns=4;
  int Nt_inv=1;
  int Nt=64;
  int TI=64;
  int nvec=6;
  bool full_tdil=true;

  int vecindex;
  for (int inoise = 0; inoise < nnoise; inoise++) {
    for (int dk = 0; dk < LI; dk++) {
      for (int dt = 0; dt < Nt_inv; dt++) {
        if(full_tdil) dt=tsrc; //TODO: this works for now, as longs as tsrc=0, but will crash otherwise!
        for (int ds = 0; ds < Ns; ds++) {
          vecindex = inoise + nnoise * dk + nnoise * LI * ds + nnoise *LI * Ns*dt;
          rho[vecindex] = zero;
          tmp3d_nospin = zero;
          for (int it = dt; it < Nt; it += TI){
	    if( it >= Ntfirst && it < Ntfirst + Ntlocal ) {
              for (int ik = dk; ik < nvec; ik += LI){
                for (int is = ds; is < Ns; is += Ns){ //at the moment, full spin dilution is enforced
                  ExtractSliceLocal(evec3d,epack.evec[ik],0,it,3);
                  tmp3d_nospin = evec3d * noise[inoise][it][ik]()(is)(); //noises do not have to be a spin vector
                  tmp3d=zero;
                  pokeSpin(tmp3d,tmp3d_nospin,is);
                  tmp2=zero;
                  InsertSliceLocal(tmp3d,tmp2,0,it-Ntfirst,Grid::QCD::Tdir);
                  rho[vecindex] += tmp2;
                }
              }
            }
          }
        }
      }
    }
  }


  for (int inoise = 0; inoise < nnoise; inoise++) {
    for (int dk = 0; dk < LI; dk++) {
      for (int dt = 0; dt < Nt_inv; dt++) {
        for (int ds = 0; ds < Ns; ds++) {
          vecindex = inoise + nnoise * dk + nnoise * LI * ds + nnoise *LI * Ns*dt;
          phi[vecindex] = zero;
          for (int t = Ntfirst; t < Ntfirst + Ntlocal; t++) {
            sink_tslice=zero;
            for (int ivec = 0; ivec < nvec; ivec++) {
              ExtractSliceLocal(evec3d,epack.evec[ivec],0,t,3);
              sink_tslice += evec3d * perambulator(t, ivec, dk, inoise,dt,ds);
            }
            InsertSliceLocal(sink_tslice,phi[vecindex],0,t-Ntfirst,Grid::QCD::Tdir);
          }
        }
      }
    }
  }


}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MDistil_DistilVectors_hpp_
