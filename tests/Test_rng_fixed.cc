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

  GridSerialRNG            fsRNG;  fsRNG.SeedFixedIntegers(seeds);
  GridParallelRNG          fpRNG(&Grid);  fpRNG.SeedFixedIntegers(seeds);

  vComplexF tmp; random(fsRNG,tmp);
  std::cout<<"Random vComplexF (fixed seed)\n"<< tmp<<std::endl;
  std::cout<<"conjugate(tmp)\n"<< conjugate(tmp)<<std::endl;
  std::cout<<"conjugate(tmp)*tmp\n"<< conjugate(tmp)*tmp<<std::endl;
  std::cout<<"innerProduct"<< innerProduct(tmp,tmp)<<std::endl;
  std::cout<<"Reduce(innerProduct)"<< Reduce(innerProduct(tmp,tmp))<<std::endl;

  SpinMatrix rnd  ; 
  random(fsRNG,rnd);
  std::cout<<"Random Spin Matrix (fixed seed)\n"<< rnd<<std::endl;

  SpinVector rv; 
  random(fsRNG,rv);
  std::cout<<"Random Spin Vector (fixed seed)\n"<< rv<<std::endl;

  gaussian(fsRNG,rv);
  std::cout<<"Gaussian Spin Vector (fixed seed)\n"<< rv<<std::endl;

  LatticeColourVector lcv(&Grid);

  LatticeFermion src(&Grid); random(fpRNG,src);
  std::cout << "src norm : " << norm2(src)<<std::endl;
  std::cout << "src " << src<<std::endl;

  random(fpRNG,lcv);
  std::cout<<"Random Lattice Colour Vector (fixed seed)\n"<< lcv<<std::endl;


  Grid_finalize();
}
