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

  latt_size.resize(4);

  latt_size[0] = 8;
  latt_size[1] = 8;
  latt_size[2] = 8;
  latt_size[3] = 8;
  double volume = latt_size[0]*latt_size[1]*latt_size[2]*latt_size[3];
  
  GridCartesian           Fine(latt_size,simd_layout,mpi_layout);


  std::vector<int> seeds({6,7,8,80});
  GridParallelRNG  pRNG(&Fine);
  pRNG.SeedFixedIntegers(seeds);

  std::vector<int> seedsS({1,2,3,4});
  GridSerialRNG    sRNG;
  sRNG.SeedFixedIntegers(seedsS);

  LatticeGaugeField U(&Fine);

  SU3::HotConfiguration(pRNG, U);

  // simplify template declaration? Strip the lorentz from the second template
  WilsonGaugeActionR Waction(6.0);

  //Collect actions
  ActionLevel<LatticeGaugeField> Level1;
  Level1.push_back(&Waction);
  ActionSet<LatticeGaugeField> FullSet;
  FullSet.push_back(Level1);

  // Create integrator
  typedef MinimumNorm2<LatticeGaugeField>  IntegratorAlgorithm;// change here to modify the algorithm
  IntegratorParameters MDpar(20);
  IntegratorAlgorithm  MDynamics(&Fine,MDpar, FullSet);

  // Create HMC
  HMCparameters HMCpar;
  HybridMonteCarlo<LatticeGaugeField,IntegratorAlgorithm>  HMC(HMCpar, MDynamics, sRNG, pRNG);

  HMC.evolve(U);

}
