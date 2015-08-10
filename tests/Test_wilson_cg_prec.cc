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

  LatticeFermion src(&Grid); random(pRNG,src);
  RealD nrm = norm2(src);
  LatticeFermion result(&Grid); result=zero;
  LatticeGaugeField Umu(&Grid); random(pRNG,Umu);

  std::vector<LatticeColourMatrix> U(4,&Grid);

  for(int mu=0;mu<Nd;mu++){
    U[mu] = PeekIndex<LorentzIndex>(Umu,mu);
  }
  
  RealD mass=0.5;
  WilsonFermionR Dw(Umu,Grid,RBGrid,mass);

  //  HermitianOperator<WilsonFermion,LatticeFermion> HermOp(Dw);
  //  ConjugateGradient<LatticeFermion> CG(1.0e-8,10000);
  //  CG(HermOp,src,result);
  
  LatticeFermion    src_o(&RBGrid);
  LatticeFermion result_o(&RBGrid);
  pickCheckerboard(Odd,src_o,src);
  result_o=zero;

  SchurDiagMooeeOperator<WilsonFermionR,LatticeFermion> HermOpEO(Dw);
  ConjugateGradient<LatticeFermion> CG(1.0e-8,10000);
  CG(HermOpEO,src_o,result_o);

  Grid_finalize();
}
