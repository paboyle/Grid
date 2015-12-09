#include <Grid.h>

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

RealD InverseApproximation(RealD x){
  return 1.0/x;
}
RealD SqrtApproximation(RealD x){
  return std::sqrt(x);
}
RealD StepFunction(RealD x){
  if ( x<0.1 )  return 1.0;
  else return 0.0;
}


int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  GridCartesian *grid = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), 
						       GridDefaultSimd(Nd,vComplex::Nsimd()),
						       GridDefaultMpi());

  double     lo=0.1;
  double     hi=64.0;

  Chebyshev<LatticeFermion> ChebyInv(lo,hi,2000,InverseApproximation);


  {
    std::ofstream of("chebyinv");
    ChebyInv.csv(of);
  }

  ChebyInv.JacksonSmooth();

  {
    std::ofstream of("chebyinvjack");
    ChebyInv.csv(of);
  }


  Chebyshev<LatticeFermion> ChebyStep(lo,hi,200,StepFunction);
  {
    std::ofstream of("chebystep");
    ChebyStep.csv(of);
  }


  ChebyStep.JacksonSmooth();

  {
    std::ofstream of("chebystepjack");
    ChebyStep.csv(of);
  }

  lo=-8;
  hi=8;
  Chebyshev<LatticeFermion> ChebyIndefInv(lo,hi,40,InverseApproximation);
  {
    std::ofstream of("chebyindefinv");
    ChebyIndefInv.csv(of);
  }

  lo=0;
  hi=64;
  Chebyshev<LatticeFermion> ChebyNE(lo,hi,40,InverseApproximation);
  {
    std::ofstream of("chebyNE");
    ChebyNE.csv(of);
  }

  Grid_finalize();
}
