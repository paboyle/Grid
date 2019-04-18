#ifndef Hadrons_MDistil_DistilSource_hpp_
#define Hadrons_MDistil_DistilSource_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Solver.hpp>
#include <Hadrons/EigenPack.hpp>
#include <Hadrons/A2AVectors.hpp>
#include <Hadrons/DilutedNoise.hpp>

// These are members of Distillation
#include <Hadrons/Modules/MDistil/Distil.hpp>

BEGIN_HADRONS_NAMESPACE


/******************************************************************************
 *                         DistilSource                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MDistil)

class DistilSourcePar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(DistilSourcePar,
		                    std::string, noise,
		                    std::string, eigenPack,
                                    bool, multiFile,
				    int, tsrc,
				    int, nnoise,
				    int, LI,
				    int, SI,
				    int, TI,
				    int, nvec,
				    int, Ns,
				    int, Nt,
				    int, Nt_inv);
};

template <typename FImpl>
class TDistilSource: public Module<DistilSourcePar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TDistilSource(const std::string name);
    // destructor
    virtual ~TDistilSource(void);
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
protected:
    // These variables are created in setup() and freed in Cleanup()
    GridCartesian * grid3d; // Owned by me, so I must delete it
    GridCartesian * grid4d; // Owned by environment (so I won't delete it)
protected:
    virtual void Cleanup(void);
};

MODULE_REGISTER_TMP(DistilSource, TDistilSource<FIMPL>, MDistil);

/******************************************************************************
 *                 TDistilSource implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TDistilSource<FImpl>::TDistilSource(const std::string name)
:  grid3d{nullptr}, grid4d{nullptr}, Module<DistilSourcePar>(name)
{}
// destructor
template <typename FImpl>
TDistilSource<FImpl>::~TDistilSource(void)
{
  Cleanup();
};

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TDistilSource<FImpl>::getInput(void)
{
    std::vector<std::string> in;

    in.push_back(par().noise);
    in.push_back(par().eigenPack);

    return in;
}

template <typename FImpl>
std::vector<std::string> TDistilSource<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDistilSource<FImpl>::setup(void)
{
    Cleanup();
   auto &noise = envGet(std::vector<Complex>, par().noise);

   int nnoise=par().nnoise;
   int LI=par().LI;
   int Ns=par().Ns;
   int SI=par().SI;
   int Nt_inv=par().Nt_inv;

   envCreate(std::vector<FermionField>, getName(), 1, 
		                    nnoise*LI*SI*Nt_inv, envGetGrid(FermionField));

  grid4d = env().getGrid();
  std::vector<int> latt_size   = GridDefaultLatt();
  std::vector<int> simd_layout = GridDefaultSimd(Nd, vComplex::Nsimd());
  std::vector<int> mpi_layout  = GridDefaultMpi();
  std::vector<int> simd_layout_3 = GridDefaultSimd(Nd-1, vComplex::Nsimd());
  latt_size[Nd-1] = 1;
  simd_layout_3.push_back( 1 );
  mpi_layout[Nd-1] = 1;
  grid3d = MakeLowerDimGrid(grid4d);


  envTmp(LatticeSpinColourVector, "tmp2",1,LatticeSpinColourVector(grid4d));
  envTmp(LatticeColourVector, "tmp_nospin",1,LatticeColourVector(grid4d));
  envTmp(LatticeSpinColourVector, "tmp3d",1,LatticeSpinColourVector(grid3d));
  envTmp(LatticeColourVector, "tmp3d_nospin",1,LatticeColourVector(grid3d));
  envTmp(LatticeSpinColourVector, "sink_tslice",1,LatticeSpinColourVector(grid3d));
  envTmp(LatticeColourVector, "evec3d",1,LatticeColourVector(grid3d));
}

// clean up any temporaries created by setup (that aren't stored in the environment)
template <typename FImpl>
void TDistilSource<FImpl>::Cleanup(void)
{
  if( grid3d != nullptr ) {
    delete grid3d;
    grid3d = nullptr;
  }
  grid4d = nullptr;
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDistilSource<FImpl>::execute(void)
{
   
    auto        &noise     = envGet(std::vector<Complex>, par().noise);
    auto        &epack   = envGet(Grid::Hadrons::EigenPack<LatticeColourVector>, par().eigenPack);
    auto        &rho       = envGet(std::vector<FermionField>, getName());

  envGetTmp(LatticeSpinColourVector, tmp2);
  envGetTmp(LatticeColourVector, tmp_nospin);
  envGetTmp(LatticeSpinColourVector, tmp3d);
  envGetTmp(LatticeColourVector, tmp3d_nospin);
  envGetTmp(LatticeSpinColourVector, sink_tslice);
  envGetTmp(LatticeColourVector, evec3d);


  int Ntlocal = grid4d->LocalDimensions()[3];
  int Ntfirst = grid4d->LocalStarts()[3];

  int tsrc=par().tsrc;
  int nnoise=par().nnoise;
  int LI=par().LI;
  int Ns=par().Ns;
  int Nt_inv=par().Nt_inv; // TODO: No input, but define through Nt, TI
  int Nt=par().Nt;
  int TI=par().TI;
  int nvec=par().nvec;
  int SI=par().SI;
  
  bool full_tdil=(TI==Nt);

  int vecindex;
  int t_inv;
  for (int inoise = 0; inoise < nnoise; inoise++) {
    for (int dk = 0; dk < LI; dk++) {
      for (int dt = 0; dt < Nt_inv; dt++) {
        for (int ds = 0; ds < SI; ds++) {
          vecindex = inoise + nnoise * dk + nnoise * LI * ds + nnoise *LI * SI*dt;
          rho[vecindex] = zero;
          tmp3d_nospin = zero;
          for (int it = dt; it < Nt; it += TI){
            if (full_tdil) t_inv = tsrc; else t_inv = it;
            if( t_inv >= Ntfirst && t_inv < Ntfirst + Ntlocal ) {
              for (int ik = dk; ik < nvec; ik += LI){
                for (int is = ds; is < Ns; is += SI){ 
                  ExtractSliceLocal(evec3d,epack.evec[ik],0,t_inv,3);
                  tmp3d_nospin = evec3d * noise[inoise + nnoise*(t_inv + Nt*(ik+nvec*is))];
                  tmp3d=zero;
                  pokeSpin(tmp3d,tmp3d_nospin,is);
                  tmp2=zero;
                  InsertSliceLocal(tmp3d,tmp2,0,t_inv-Ntfirst,Grid::QCD::Tdir);
                  rho[vecindex] += tmp2;
                }
              }
            }
          }
        }
      }
    }
  }


  // TEST TO SEE WHETHER THIS MIGHT BE THE MEMORY LEAK
  Cleanup();

}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MDistil_DistilSource_hpp_
