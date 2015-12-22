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

    Real mass=-0.77;
    FermionAction FermOp(U,*FGrid,*FrbGrid,mass);

    // 1+1 flavour
    OneFlavourRationalParams Params(1.0e-4,64.0,1000,1.0e-6);
    OneFlavourRationalPseudoFermionAction<WilsonImplR> WilsonNf1a(FermOp,Params);
    OneFlavourRationalPseudoFermionAction<WilsonImplR> WilsonNf1b(FermOp,Params);

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

