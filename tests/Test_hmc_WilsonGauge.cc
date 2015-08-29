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
  GridParallelRNG  pRNG(&Fine);
  std::vector<int> seeds({1,2,3,4,5,6,7,8});
  pRNG.SeedFixedIntegers(seeds);
  LatticeGaugeField U(&Fine);

  SU3::HotConfiguration(pRNG, U);
 

  // simplify template declaration? Strip the lorentz from the second template
  WilsonGaugeAction<LatticeGaugeField, LatticeColourMatrix> Waction(6.0);

  //Collect actions
  ActionLevel Level1;
  Level1.push_back(&Waction);
  ActionSet FullSet;
  FullSet.push_back(Level1);

  // Create integrator
  typedef MinimumNorm2  IntegratorAlgorithm;// change here to modify the algorithm
  IntegratorParameters MDpar(12,20,1.0);
  Integrator<IntegratorAlgorithm> MDynamics(&Fine,MDpar, FullSet);

  // Create HMC
  HMCparameters HMCpar;
  HybridMonteCarlo<IntegratorAlgorithm>  HMC(HMCpar, MDynamics);

  HMC.evolve(U);

}
