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

  const int Ls=9;

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
  LatticeFermion    ref(FGrid); ref=zero;
  LatticeFermion    tmp(FGrid);
  LatticeFermion    err(FGrid);
  LatticeGaugeField Umu(UGrid); random(RNG4,Umu);

  std::vector<LatticeColourMatrix> U(4,UGrid);
  for(int mu=0;mu<Nd;mu++){
    U[mu] = PeekIndex<LorentzIndex>(Umu,mu);
  }
  
  RealD mass=0.1;
  RealD M5=1.8;

  {
    OverlapWilsonContFracTanhFermion Dcf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,1.0);
    HermitianLinearOperator<OverlapWilsonContFracTanhFermion,LatticeFermion> HermIndefOp(Dcf);

    HermIndefOp.Op(src,ref);
    HermIndefOp.OpDiag(src,result);
    
    for(int d=0;d<4;d++){
      HermIndefOp.OpDir(src,tmp,d,+1); result=result+tmp; 
      std::cout<<"dir "<<d<<" tmp "<<norm2(tmp)<<std::endl;
      HermIndefOp.OpDir(src,tmp,d,-1); result=result+tmp;
      std::cout<<"dir "<<d<<" tmp "<<norm2(tmp)<<std::endl;
    }
    err = result-ref;
    std::cout<<"Error "<<norm2(err)<<std::endl;
  }

  {
    OverlapWilsonPartialFractionTanhFermion Dpf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,1.0);
    HermitianLinearOperator<OverlapWilsonPartialFractionTanhFermion,LatticeFermion> HermIndefOp(Dpf);
    
    HermIndefOp.Op(src,ref);
    HermIndefOp.OpDiag(src,result);
    
    for(int d=0;d<4;d++){
      HermIndefOp.OpDir(src,tmp,d,+1); result=result+tmp; 
      std::cout<<"dir "<<d<<" tmp "<<norm2(tmp)<<std::endl;
      HermIndefOp.OpDir(src,tmp,d,-1); result=result+tmp;
      std::cout<<"dir "<<d<<" tmp "<<norm2(tmp)<<std::endl;
    }

    err = result-ref;
    std::cout<<"Error "<<norm2(err)<<std::endl;
  }


  Grid_finalize();
}
