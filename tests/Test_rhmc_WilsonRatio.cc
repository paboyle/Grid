#include "Grid.h"

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

namespace Grid { 
  namespace QCD { 


class HmcRunner : public NerscHmcRunner {
public:

  void BuildTheAction (int argc, char **argv)

  {
    typedef WilsonImplR ImplPolicy;
    typedef WilsonFermionR FermionAction;
    typedef typename FermionAction::FermionField FermionField;

    UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
    UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  
    FGrid   = UGrid;
    FrbGrid = UrbGrid;

    // temporarily need a gauge field
    LatticeGaugeField  U(UGrid);

    // Gauge action
    WilsonGaugeActionR Waction(5.6);

    RealD mass=-0.77;
    RealD pv  =0.0;
    FermionAction DenOp(U,*FGrid,*FrbGrid,mass);
    FermionAction NumOp(U,*FGrid,*FrbGrid,pv);
  

  // erange,maxiter,resid,npoly
    OneFlavourRationalParams Params(1.0e-2,64.0,1000,1.0e-6,6);
    OneFlavourRatioRationalPseudoFermionAction<WilsonImplR> WilsonNf1a(NumOp,DenOp,Params);
    OneFlavourRatioRationalPseudoFermionAction<WilsonImplR> WilsonNf1b(NumOp,DenOp,Params);
  
    //Collect actions
    ActionLevel<LatticeGaugeField> Level1;
    Level1.push_back(&WilsonNf1a);
    Level1.push_back(&WilsonNf1b);
    Level1.push_back(&Waction);

    TheAction.push_back(Level1);

    Run(argc,argv);
  };

};

}}

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  int threads = GridThread::GetThreads();
  std::cout<<GridLogMessage << "Grid is setup to use "<<threads<<" threads"<<std::endl;

  HmcRunner TheHMC;
  
  TheHMC.BuildTheAction(argc,argv);

}

