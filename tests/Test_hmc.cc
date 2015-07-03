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
  GridParallelRNG       FineRNG(&Fine);
  GridSerialRNG       SerialRNG;
  FineRNG.SeedRandomDevice();
  
  WilsonGaugeAction<LatticeLorentzColourMatrix, LatticeColourMatrix> Waction(6.0);
  //Collect actions
  ActionLevel Level;
  Level.push_back(&Waction);
  

  //Integrator<IntegratorLeapFrog(12,10,1.0);
}
