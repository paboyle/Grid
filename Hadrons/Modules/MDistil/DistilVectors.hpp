#ifndef Hadrons_MDistil_DistilVectors_hpp_
#define Hadrons_MDistil_DistilVectors_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *      Class to generate Distillation $\varrho$ and $\varphi$ vectors        *
 ******************************************************************************/
/*
 template <typename FImpl>
 class DistilVectorsFromPerambulator
 {
 public:
   FERM_TYPE_ALIASES(FImpl,);
   SOLVER_TYPE_ALIASES(FImpl,);
 public:
   DistilVectorsFromPerambulator(FMat &action);
   virtual ~DistilVectorsFromPerambulator(void) = default;
   void makeRho(FermionField &rhoOut,
   const FermionField &evec, const Complex &noises);
   void makePhi(FermionField &phiOut,
   const FermionField &evec, const SpinVector &Perambulators);
 private:
   GridBase                                 *fGrid_, *frbGrid_, *gGrid_;
   bool                                     is5d_;
   FermionField                             src_o_, sol_e_, sol_o_, tmp_, tmp5_;
 };

*/


/******************************************************************************
 *               A2AVectorsSchurDiagTwo template implementation               *
 ******************************************************************************/
/*

template <typename FImpl>
A2AVectorsSchurDiagTwo<FImpl>::A2AVectorsSchurDiagTwo(FMat &action)
: action_(action)
, fGrid_(action_.FermionGrid())
, frbGrid_(action_.FermionRedBlackGrid())
, gGrid_(action_.GaugeGrid())
, src_o_(frbGrid_)
, sol_e_(frbGrid_)
, sol_o_(frbGrid_)
, tmp_(frbGrid_)
, tmp5_(fGrid_)
, op_(action_)
{}

template <typename FImpl>
void A2AVectorsSchurDiagTwo<FImpl>::makeRho(FermionField &rhoOut, const FermionField &evec, const Complex &noises)
{

  LatticeSpinColourVector tmp2(grid4d);
  LatticeSpinColourVector tmp3d(grid3d);
  LatticeColourVector tmp3d_nospin(grid3d);
  LatticeColourVector tmp_nospin(grid4d);
  LatticeColourVector evec3d(grid3d);


  for (int inoise = 0; inoise < nnoise; inoise++) {
    for (int dk = 0; dk < LI; dk++) {
      for (int dt = 0; dt < Nt_inv; dt++) {
        if(full_tdil) dt=tsrc; //TODO: this works for now, as longs as tsrc=0, but will crash otherwise!
        for (int ds = 0; ds < Ns; ds++) {
          vecindex = inoise + nnoise * dk + nnoise * LI * ds + nnoise *LI * Ns*dt;
          sources_tsrc[vecindex] = zero;
          tmp3d_nospin = zero;
          for (int it = dt; it < Nt; it += TI){
            if( it >= Ntfirst && it < Ntfirst + Ntlocal ) {
              for (int ik = dk; ik < nvec; ik += LI){
                for (int is = ds; is < Ns; is += Ns){ //at the moment, full spin dilution is enforced
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
                  sources_tsrc[vecindex] += tmp2;
                }
              }
            }
          }
        }
      }
    }
  }


}

*/



/******************************************************************************
 *                         DistilVectors                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MDistil)

class DistilVectorsPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(DistilVectorsPar,
                                    unsigned int, i);
};

template <typename FImpl>
class TDistilVectors: public Module<DistilVectorsPar>
{
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
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TDistilVectors<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDistilVectors<FImpl>::setup(void)
{
    
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDistilVectors<FImpl>::execute(void)
{
    
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MDistil_DistilVectors_hpp_
