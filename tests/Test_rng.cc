#include <Grid.h>

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  std::vector<int> latt_size   = GridDefaultLatt();
  std::vector<int> simd_layout = GridDefaultSimd(4,vComplex::Nsimd());
  std::vector<int> mpi_layout  = GridDefaultMpi();
     
  GridCartesian     Grid(latt_size,simd_layout,mpi_layout);

  std::vector<int> seeds({1,2,3,4});

  GridSerialRNG             sRNG;   sRNG.SeedRandomDevice();
  GridSerialRNG            fsRNG;  fsRNG.SeedFixedIntegers(seeds);

  GridParallelRNG           pRNG(&Grid);   pRNG.SeedRandomDevice();
  GridParallelRNG          fpRNG(&Grid);  fpRNG.SeedFixedIntegers(seeds);

  SpinMatrix rnd  ; 
  random(sRNG,rnd);
  std::cout<<GridLogMessage<<"Random Spin Matrix (random_device)\n"<< rnd<<std::endl;

  random(fsRNG,rnd);
  std::cout<<GridLogMessage<<"Random Spin Matrix (fixed seed)\n"<< rnd<<std::endl;

  SpinVector rv; 
  random(sRNG,rv);
  std::cout<<GridLogMessage<<"Random Spin Vector (random device)\n"<< rv<<std::endl;

  random(fsRNG,rv);
  std::cout<<GridLogMessage<<"Random Spin Vector (fixed seed)\n"<< rv<<std::endl;

  gaussian(fsRNG,rv);
  std::cout<<GridLogMessage<<"Gaussian Spin Vector (fixed seed)\n"<< rv<<std::endl;

  LatticeColourVector lcv(&Grid);
  random(pRNG,lcv);
  std::cout<<GridLogMessage<<"Random Lattice Colour Vector (random device)\n"<< lcv<<std::endl;

  random(fpRNG,lcv);
  std::cout<<GridLogMessage<<"Random Lattice Colour Vector (fixed seed)\n"<< lcv<<std::endl;

  Grid_finalize();
}
