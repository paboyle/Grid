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
  
    ConjugateGradient<FermionField>  CG(1.0e-8,10000);

    TwoFlavourPseudoFermionAction<ImplPolicy> Nf2(FermOp,CG,CG);
  
    //Collect actions
    ActionLevel<LatticeGaugeField> Level1(1);
    Level1.push_back(&Nf2);

    ActionLevel<LatticeGaugeField> Level2(4);
    Level2.push_back(&Waction);

    TheAction.push_back(Level1);
    TheAction.push_back(Level2);

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


