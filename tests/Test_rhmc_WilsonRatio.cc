#include "Grid.h"


using namespace std;
using namespace Grid;
using namespace Grid::QCD;

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  std::vector<int> latt_size   = GridDefaultLatt();
  std::vector<int> simd_layout = GridDefaultSimd(4,vComplex::Nsimd());
  std::vector<int> mpi_layout  = GridDefaultMpi();
  
  GridCartesian            Fine(latt_size,simd_layout,mpi_layout);
  GridRedBlackCartesian  RBFine(latt_size,simd_layout,mpi_layout);


  std::vector<int> seeds({6,7,8,80});
  GridParallelRNG  pRNG(&Fine);
  pRNG.SeedFixedIntegers(seeds);

  std::vector<int> seedsS({1,2,3,4});
  GridSerialRNG    sRNG;
  sRNG.SeedFixedIntegers(seedsS);

  LatticeLorentzColourMatrix     U(&Fine);

  SU3::HotConfiguration(pRNG, U);

  // simplify template declaration? Strip the lorentz from the second template
  WilsonGaugeActionR Waction(5.6);

  Real mass=-0.77;
  Real pv  =0.0;
  WilsonFermionR DenOp(U,Fine,RBFine,mass);
  WilsonFermionR NumOp(U,Fine,RBFine,pv);
  
  // erange,maxiter,resid,npoly
  OneFlavourRationalParams Params(1.0e-2,64.0,1000,1.0e-6,6);
  OneFlavourRatioRationalPseudoFermionAction<WilsonImplR> WilsonNf1a(NumOp,DenOp,Params);
  OneFlavourRatioRationalPseudoFermionAction<WilsonImplR> WilsonNf1b(NumOp,DenOp,Params);
  
  //Collect actions
  ActionLevel<LatticeGaugeField> Level1;
  Level1.push_back(&WilsonNf1a);
  Level1.push_back(&WilsonNf1b);
  Level1.push_back(&Waction);
  ActionSet<LatticeGaugeField> FullSet;
  FullSet.push_back(Level1);


  // Create integrator
  typedef MinimumNorm2<LatticeGaugeField>  IntegratorAlgorithm;// change here to change the algorithm
  //  typedef LeapFrog  IntegratorAlgorithm;// change here to change the algorithm
  IntegratorParameters MDpar(20);
  IntegratorAlgorithm MDynamics(&Fine,MDpar, FullSet);

  // Create HMC
  HMCparameters HMCpar;
  HybridMonteCarlo<LatticeGaugeField,IntegratorAlgorithm>  HMC(HMCpar, MDynamics, sRNG, pRNG);

  HMC.evolve(U);

}
