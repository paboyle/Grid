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

  const int Ls=4;
  const int L=4;
  std::vector<int> latt_2f(Nd,L);
  std::vector<int> latt_1f(Nd,L); latt_1f[0] = 2L;

  std::vector<int> simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
  std::vector<int> mpi_layout  = GridDefaultMpi(); //node layout

  GridCartesian         * UGrid_1f   = SpaceTimeGrid::makeFourDimGrid(latt_1f, simd_layout, mpi_layout);
  GridRedBlackCartesian * UrbGrid_1f = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid_1f);
  GridCartesian         * FGrid_1f   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid_1f);
  GridRedBlackCartesian * FrbGrid_1f = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid_1f);

  std::vector<int> seeds4({1,2,3,4});
  std::vector<int> seeds5({5,6,7,8});
  GridParallelRNG          RNG5_1f(FGrid_1f);  RNG5_1f.SeedFixedIntegers(seeds5);
  GridParallelRNG          RNG4_1f(UGrid_1f);  RNG4_1f.SeedFixedIntegers(seeds4);

  LatticeFermion    src_1f(FGrid_1f); random(RNG5_1f,src_1f);
  LatticeFermion result_1f(FGrid_1f); result_1f=zero;
  LatticeGaugeField Umu_1f(UGrid_1f); random(RNG4_1f,Umu_1f);

  //Coordinate grid for reference
  LatticeInteger xcoor_1f(UGrid_1f);
  LatticeCoordinate(xcoor_1f,0);

  //Copy-conjugate the gauge field
  //First C-shift the lattice by Lx/2
  {
    LatticeGaugeField Umu_shift = conjugate( Cshift(Umu_1f,0,L) );
    Umu_1f = where( xcoor_1f >= Integer(L), Umu_shift, Umu_1f );
  }

  //Make the gauge field antiperiodic in x-direction
  Umu_1f = where(xcoor_1f == Integer(L-1), -Umu_1f, Umu_1f);
    
  RealD mass=0.1;
  RealD M5=1.8;
  DomainWallFermionR Ddwf(Umu_1f,*FGrid_1f,*FrbGrid_1f,*UGrid_1f,*UrbGrid_1f,mass,M5);

  LatticeFermion    src_o_1f(FrbGrid_1f);
  LatticeFermion result_o_1f(FrbGrid_1f);
  pickCheckerboard(Odd,src_o_1f,src_1f);
  result_o_1f=zero;

  SchurDiagMooeeOperator<DomainWallFermionR,LatticeFermion> HermOpEO(Ddwf);
  ConjugateGradient<LatticeFermion> CG(1.0e-8,10000);
  CG(HermOpEO,src_o_1f,result_o_1f);

  Grid_finalize();
}
