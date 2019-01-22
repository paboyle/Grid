#ifndef Hadrons_MDistil_PerambLight_hpp_
#define Hadrons_MDistil_PerambLight_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Solver.hpp>
#include <Hadrons/EigenPack.hpp>
#include <Hadrons/A2AVectors.hpp>
#include <Hadrons/DilutedNoise.hpp>
#include <Hadrons/DistilVectors.hpp>

BEGIN_HADRONS_NAMESPACE


/******************************************************************************
 *                             PerambLight                                    *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MDistil)

class PerambLightPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(PerambLightPar,
		                    std::string, noise,
		                    std::string, eigenPack,
                                    bool, multiFile);
};

template <typename FImpl>
class TPerambLight: public Module<PerambLightPar>
{
public:
    // constructor
    TPerambLight(const std::string name);
    // destructor
    virtual ~TPerambLight(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(PerambLight, TPerambLight<FIMPL>, MDistil);

/******************************************************************************
 *                 TPerambLight implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TPerambLight<FImpl>::TPerambLight(const std::string name)
: Module<PerambLightPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TPerambLight<FImpl>::getInput(void)
{
    std::vector<std::string> in;

    in.push_back(par().noise);
    in.push_back(par().eigenPack);
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TPerambLight<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName() + "_perambulator_light"};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TPerambLight<FImpl>::setup(void)
{

   auto &noise = envGet(std::vector<std::vector<std::vector<SpinVector>>>, par().noise);
  
   int nvec = 6;
   int Nt=64;

   envCreate(Perambulator<SpinVector>, getName() + "_perambulator_light", 1, 
		                    noise.size() *nvec*Nt);

}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TPerambLight<FImpl>::execute(void)
{
/*
    for (int inoise = 0; inoise < nnoise; inoise++) {
      for (int dk = 0; dk < LI; dk++) {
        for (int dt = 0; dt < Nt_inv; dt++) {
          if(full_tdil) dt=tsrc; //this works for now, as longs as tsrc=0, but will crash otherwise!
          for (int ds = 0; ds < Ns; ds++) {
            std::cout <<  "LapH source vector from noise " << inoise << " and dilution component (d_k,d_t,d_alpha) : (" << dk << ","<< dt << "," << ds << ")" << std::endl;
            dist_source = zero;
            tmp3d_nospin = zero;
            evec3d = zero;
            for (int it = dt; it < Nt; it += TI){
              if( it >= Ntfirst && it < Ntfirst + Ntlocal ) {
                for (int ik = dk; ik < nvec; ik += LI){
                  for (int is = ds; is < Ns; is += Ns){ //at the moment, full spin dilution is enforced
                    std::cout <<  "LapH source vector from noise " << it << " and dilution component (d_k,d_t,d_alpha) : (" << ik << ","<< is << ")" << std::endl;
                    ExtractSliceLocal(evec3d,eig4d.evec[ik],0,it,3);
                    tmp3d_nospin = evec3d * noises[inoise][it][ik]()(is)(); //noises do not have to be a spin vector
                    tmp3d=zero;
                    pokeSpin(tmp3d,tmp3d_nospin,is);
                    tmp2=zero;
#ifdef USE_LOCAL_SLICES
                    InsertSliceLocal(tmp3d,tmp2,0,it-Ntfirst,Grid::QCD::Tdir);
#else
                    InsertSlice(tmp3d,tmp2,it,Grid::QCD::Tdir);
#endif
                    dist_source += tmp2;
                  }
                }
              }
            }
            std::cout <<  "Inversion for noise " << inoise << " and dilution component (d_k,d_t,d_alpha) : (" << dk << ","<< dt << "," << ds << ")" << std::endl;
            result=zero;
            LatticeFermion src5(FGrid);
            LatticeFermion sol5(FGrid);
            Dop.ImportPhysicalFermionSource(dist_source,src5);
            SchurSolver(Dop,src5,sol5);
            Dop.ExportPhysicalFermionSolution(sol5,result); //These are the meson sinks
            if (compute_current_sink)
              current_sink[inoise+nnoise*(dk+LI*(dt+Nt_inv*ds))] = result;
            std::cout <<  "Contraction of perambulator from noise " << inoise << " and dilution component (d_k,d_t,d_alpha) : (" << dk << ","<< dt << "," << ds << ")" << std::endl;
            for (int is = 0; is < Ns; is++) {
              result_nospin = peekSpin(result,is);
              for (int t = Ntfirst; t < Ntfirst + Ntlocal; t++) {
#ifdef USE_LOCAL_SLICES
                ExtractSliceLocal(result_3d,result_nospin,0,t-Ntfirst,Grid::QCD::Tdir);
#else
                ExtractSlice(result_3d,result_nospin,t,3);
#endif
                for (int ivec = 0; ivec < nvec; ivec++) {
                  ExtractSliceLocal(evec3d,eig4d.evec[ivec],0,t,3);
                  pokeSpin(perambulator(t, ivec, dk, inoise,dt,ds),innerProduct(evec3d, result_3d),is);
                }
          }
        }
      }
    }
    std::cout <<  "perambulator done" << std::endl;
    perambulator.SliceShare( grid3d, grid4d );

    // THIS IS WHERE WE WANT TO SAVE THE PERAMBULATORS TO DISK
    perambulator.WriteTemporary(std::string(pszPerambPack));
*/
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MDistil_PerambLight_hpp_
