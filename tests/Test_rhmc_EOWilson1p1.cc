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
  GridParallelRNG  pRNG(&Fine);
  pRNG.SeedRandomDevice();
  LatticeLorentzColourMatrix     U(&Fine);

  SU3::HotConfiguration(pRNG, U);

  // simplify template declaration? Strip the lorentz from the second template
  WilsonGaugeAction<LatticeLorentzColourMatrix, LatticeColourMatrix> Waction(5.6);

  Real mass=-0.77;
  WilsonFermionR FermOp(U,Fine,RBFine,mass);

  // 1+1 flavour
  OneFlavourRationalParams Params(1.0e-4,64.0,1000,1.0e-6);
  OneFlavourEvenOddRationalPseudoFermionAction<WilsonImplR> WilsonNf1a(FermOp,Params);
  OneFlavourEvenOddRationalPseudoFermionAction<WilsonImplR> WilsonNf1b(FermOp,Params);
  
  //Collect actions
  ActionLevel Level1;
  Level1.push_back(&WilsonNf1a);
  Level1.push_back(&WilsonNf1b);
  Level1.push_back(&Waction);

  ActionSet FullSet;
  FullSet.push_back(Level1);

  // Create integrator
  typedef MinimumNorm2  IntegratorAlgorithm;// change here to change the algorithm

  IntegratorParameters MDpar(12,20,1.0);
  Integrator<IntegratorAlgorithm> MDynamics(&Fine,MDpar, FullSet);

  // Create HMC
  HMCparameters HMCpar;
  HybridMonteCarlo<IntegratorAlgorithm>  HMC(HMCpar, MDynamics);
  HMC.evolve(U);

}
