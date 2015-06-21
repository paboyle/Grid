#include <Grid.h>
#include <algorithms/iterative/PrecGeneralisedConjugateResidual.h>
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



  const int Ls=8;

  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplexF::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);


  std::vector<int> seeds4({1,2,3,4});
  std::vector<int> seeds5({5,6,7,8});
  GridParallelRNG          RNG5(FGrid);  RNG5.SeedFixedIntegers(seeds5);
  GridParallelRNG          RNG4(UGrid);  RNG4.SeedFixedIntegers(seeds4);

  LatticeFermion    src(FGrid); random(RNG5,src);
  LatticeFermion result(FGrid); result=zero;
  LatticeGaugeField Umu(UGrid); 

  SU3::HotConfiguration(RNG4,Umu);

  TrivialPrecon<LatticeFermion> simple;

  PrecGeneralisedConjugateResidual<LatticeFermion> PGCR(1.0e-6,10000,simple,4,160);

  ConjugateResidual<LatticeFermion> CR(1.0e-6,10000);

  ConjugateGradient<LatticeFermion> CG(1.0e-6,10000);
  
  RealD mass=0.5;
  RealD M5=1.8;
  DomainWallFermion Ddwf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);

  std::cout<<"*********************************************************"<<std::endl;
  std::cout<<"* Solving with MdagM VPGCR "<<std::endl;
  std::cout<<"*********************************************************"<<std::endl;
  MdagMLinearOperator<DomainWallFermion,LatticeFermion> HermOp(Ddwf);
  result=zero;
  PGCR(HermOp,src,result);

  std::cout<<"*********************************************************"<<std::endl;
  std::cout<<"* Solving with g5-VPGCR "<<std::endl;
  std::cout<<"*********************************************************"<<std::endl;
  Gamma5R5HermitianLinearOperator<DomainWallFermion,LatticeFermion> g5HermOp(Ddwf);
  result=zero;
  PGCR(g5HermOp,src,result);

  std::cout<<"*********************************************************"<<std::endl;
  std::cout<<"* Solving with MdagM-CR "<<std::endl;
  std::cout<<"*********************************************************"<<std::endl;
  result=zero;
  CR(HermOp,src,result);

  std::cout<<"*********************************************************"<<std::endl;
  std::cout<<"* Solving with g5-CR "<<std::endl;
  std::cout<<"*********************************************************"<<std::endl;
  result=zero;
  CR(g5HermOp,src,result);

  std::cout<<"*********************************************************"<<std::endl;
  std::cout<<"* Solving with MdagM-CG "<<std::endl;
  std::cout<<"*********************************************************"<<std::endl;
  result=zero;
  CG(HermOp,src,result);

  Grid_finalize();
}
