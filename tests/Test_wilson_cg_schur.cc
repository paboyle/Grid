#include <Grid.h>

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

template<class d>
struct scal {
  d internal;
};

  Gamma::GammaMatrix Gmu [] = {
    Gamma::GammaX,
    Gamma::GammaY,
    Gamma::GammaZ,
    Gamma::GammaT
  };

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);


  std::vector<int> latt_size   = GridDefaultLatt();
  std::vector<int> simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
  std::vector<int> mpi_layout  = GridDefaultMpi();
  GridCartesian               Grid(latt_size,simd_layout,mpi_layout);
  GridRedBlackCartesian     RBGrid(latt_size,simd_layout,mpi_layout);

  std::vector<int> seeds({1,2,3,4});
  GridParallelRNG          pRNG(&Grid);  pRNG.SeedFixedIntegers(seeds);

  LatticeGaugeField Umu(&Grid); random(pRNG,Umu);

  LatticeFermion    src(&Grid); random(pRNG,src);
  LatticeFermion result(&Grid); result=zero;
  LatticeFermion resid(&Grid); 

  RealD mass=0.5;
  WilsonFermionR Dw(Umu,Grid,RBGrid,mass);

  ConjugateGradient<LatticeFermion> CG(1.0e-8,10000);
  SchurRedBlackDiagMooeeSolve<LatticeFermion> SchurSolver(CG);

  SchurSolver(Dw,src,result);
  
  Grid_finalize();
}
