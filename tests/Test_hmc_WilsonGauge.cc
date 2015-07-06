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
  pRNG.SeedRandomDevice();
  LatticeLorentzColourMatrix     U(&Fine);

  SU3::HotConfiguration(pRNG, U);
 

  // simplify template?
  WilsonGaugeAction<LatticeLorentzColourMatrix, LatticeColourMatrix> Waction(6.0);

  //Collect actions
  ActionLevel Level1;
  Level1.push_back(&Waction);
  ActionSet FullSet;
  FullSet.push_back(Level1);

  // Create integrator
  IntegratorParameters MDpar(12,30,1.0);
  std::vector<int> rel ={1};
  Integrator<LeapFrog> MDleapfrog(MDpar, FullSet,rel);

  // Create HMC
  HMCparameters HMCpar;
  HybridMonteCarlo<LeapFrog>  HMC(HMCpar, MDleapfrog, &Fine);

  HMC.evolve(U);

}
