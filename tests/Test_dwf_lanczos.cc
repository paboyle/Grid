#include <Grid.h>

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  const int Ls=8;

  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  std::vector<int> seeds4({1,2,3,4});
  std::vector<int> seeds5({5,6,7,8});
  GridParallelRNG          RNG5(FGrid);  RNG5.SeedFixedIntegers(seeds5);
  GridParallelRNG          RNG4(UGrid);  RNG4.SeedFixedIntegers(seeds4);

  LatticeFermion    src(FGrid); gaussian(RNG5,src);
  LatticeGaugeField Umu(UGrid); 
  SU3::HotConfiguration(RNG4, Umu);

  std::vector<LatticeColourMatrix> U(4,UGrid);
  for(int mu=0;mu<Nd;mu++){
    U[mu] = PeekIndex<LorentzIndex>(Umu,mu);
  }
  
  RealD mass=0.1;
  RealD M5=1.8;
  DomainWallFermionR Ddwf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);

  MdagMLinearOperator<DomainWallFermionR,LatticeFermion> HermOp(Ddwf);

  const int Nk = 30;
  const int Np = 10;
  const int Nm = Nk+Np;
  const int MaxIt= 10000;
  RealD resid = 1.0e-8;

  std::vector<double> Coeffs(1,1.0);
  Polynomial<LatticeFermion> PolyX(Coeffs);
  ImplicitlyRestartedLanczos<LatticeFermion> IRL(HermOp,PolyX,Nk,Nm,resid,MaxIt);

  
  std::vector<RealD>          eval(Nm);
  std::vector<LatticeFermion> evec(Nm,FGrid);
  for(int i=0;i<Nm;i++){
    std::cout << i<<" / "<< Nm<< " grid pointer "<<evec[i]._grid<<std::endl;
  };

  int Nconv;
  IRL.calc(eval,evec,
	   src,
	   Nconv);


  Grid_finalize();
}
