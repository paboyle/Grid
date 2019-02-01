#ifndef Hadrons_MDistil_DistilVectors_hpp_
#define Hadrons_MDistil_DistilVectors_hpp_

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
   //auto &noise = envGet(std::vector<std::vector<std::vector<SpinVector>>>, par().noise);
   auto &noise = envGet(std::vector<Complex>, par().noise);

   int nnoise=par().nnoise;
   int LI=par().LI;
   int Ns=par().Ns;
   int Nt_inv=par().Nt_inv;

   envCreate(std::vector<FermionField>, getName() + "_rho", 1, 
		                    nnoise*LI*Ns*Nt_inv, envGetGrid(FermionField));
   envCreate(std::vector<FermionField>, getName() + "_phi", 1, 
                 	            nnoise*LI*Ns*Nt_inv, envGetGrid(FermionField)); 


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
  envTmp(LatticeColourVector, "evec3d",1,LatticeColourVector(grid3d));
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDistilVectors<FImpl>::execute(void)
{
   
    //auto        &noise     = envGet(std::vector<std::vector<std::vector<SpinVector>>>, par().noise);
    auto        &noise     = envGet(std::vector<Complex>, par().noise);
    auto        &perambulator   = envGet(Perambulator<SpinVector COMMA 6>, par().perambulator);
    auto        &epack   = envGet(Grid::Hadrons::EigenPack<LatticeColourVector>, par().eigenPack);
    auto        &rho       = envGet(std::vector<FermionField>, getName() + "_rho");
    auto        &phi       = envGet(std::vector<FermionField>, getName() + "_phi");


  envGetTmp(LatticeSpinColourVector, tmp2);
  envGetTmp(LatticeColourVector, tmp_nospin);
  envGetTmp(LatticeSpinColourVector, tmp3d);
  envGetTmp(LatticeColourVector, tmp3d_nospin);
  envGetTmp(LatticeSpinColourVector, sink_tslice);
  envGetTmp(LatticeColourVector, evec3d);

  GridCartesian * grid4d = env().getGrid();

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
  
  bool full_tdil=(TI==Nt);

  int vecindex;
  int t_inv;
  for (int inoise = 0; inoise < nnoise; inoise++) {
    for (int dk = 0; dk < LI; dk++) {
      for (int dt = 0; dt < Nt_inv; dt++) {
        for (int ds = 0; ds < Ns; ds++) {
          vecindex = inoise + nnoise * dk + nnoise * LI * ds + nnoise *LI * Ns*dt;
          rho[vecindex] = zero;
          tmp3d_nospin = zero;
          for (int it = dt; it < Nt; it += TI){
            if (full_tdil) t_inv = tsrc; else t_inv = it;
            if( t_inv >= Ntfirst && t_inv < Ntfirst + Ntlocal ) {
              for (int ik = dk; ik < nvec; ik += LI){
                for (int is = ds; is < Ns; is += Ns){ //at the moment, full spin dilution is enforced
                  ExtractSliceLocal(evec3d,epack.evec[ik],0,t_inv,3);
                  tmp3d_nospin = evec3d * noise[inoise + nnoise*(t_inv + Nt*(ik+nvec*is))];
                  //tmp3d_nospin = evec3d * noise[inoise][it][ik]()(is)(); //noises do not have to be a spin vector
                  tmp3d=zero;
                  pokeSpin(tmp3d,tmp3d_nospin,is);
                  tmp2=zero;
                  InsertSliceLocal(tmp3d,tmp2,0,t_inv-Ntfirst,Grid::QCD::Tdir);
                  rho[vecindex] += tmp2;
                }
              }
            }
          }
          SpinColourVector scv0;
          std::vector<int> siteFirst(grid4d->Nd(),0);
          peekSite(scv0, rho[vecindex], siteFirst);
          auto & cplx0 = scv0()(0)(0);
          std::cout <<  "site0[rho(vecindex = " << vecindex << ")] = " << cplx0 << std::endl;
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
          SpinColourVector scv0;
          std::vector<int> siteFirst(grid4d->Nd(),0);
          peekSite(scv0, phi[vecindex], siteFirst);
          auto & cplx0 = scv0()(0)(0);
          std::cout <<  "site0[phi(vecindex = " << vecindex << ")] = " << cplx0 << std::endl;
        }
      }
    }
  }

  std::cout << "size rho" << rho.size() << std::endl;
  std::cout << "size phi" << phi.size() << std::endl;

}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MDistil_DistilVectors_hpp_
