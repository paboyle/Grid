#include "Grid.h"


using namespace std;
using namespace Grid;
using namespace Grid::QCD;

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  const int Ls = 8;

  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);

  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  GridSerialRNG    sRNG;
  GridParallelRNG  pRNG(UGrid);
  sRNG.SeedRandomDevice();
  pRNG.SeedRandomDevice();

  LatticeLorentzColourMatrix     U(UGrid);

  SU3::HotConfiguration(pRNG, U);

  // simplify template declaration? Strip the lorentz from the second template
  WilsonGaugeActionR Waction(5.6);

  Real mass=0.04;
  Real pv  =1.0;
  RealD M5=1.5;

  typedef typename GparityDomainWallFermionR::FermionField FermionField;

  const int nu = 3;
  std::vector<int> twists(Nd,0);
  twists[nu] = 1;
  GparityDomainWallFermionR::ImplParams params;
  params.twists = twists;

  GparityDomainWallFermionR DenOp(U,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,params);
  GparityDomainWallFermionR NumOp(U,*FGrid,*FrbGrid,*UGrid,*UrbGrid,pv,M5,params);
  
  ConjugateGradient<FermionField>  CG(1.0e-8,10000);

  TwoFlavourEvenOddRatioPseudoFermionAction<GparityWilsonImplR> Nf2(NumOp, DenOp,CG,CG);
  
  //Collect actions
  ActionLevel<LatticeGaugeField> Level1;
  Level1.push_back(&Nf2);
  Level1.push_back(&Waction);
  ActionSet<LatticeGaugeField> FullSet;
  FullSet.push_back(Level1);

  // Create integrator
  typedef MinimumNorm2<LatticeGaugeField>  IntegratorType;// change here to change the algorithm
  IntegratorParameters MDpar(20);
  IntegratorType MDynamics(UGrid,MDpar, FullSet);

  // Create HMC
  NerscHmcCheckpointer<LatticeGaugeField> Checkpoint(std::string("ckpoint_lat"),std::string("ckpoint_rng"),1);
  HMCparameters HMCpar;
  HybridMonteCarlo<LatticeGaugeField,IntegratorType>  HMC(HMCpar, MDynamics,sRNG,pRNG,U);
  HMC.AddObservable(&Checkpoint);

  // Run it
  HMC.evolve();


}
