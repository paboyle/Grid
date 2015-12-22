#include "Grid.h"
using namespace std;
using namespace Grid;
using namespace Grid::QCD;

namespace Grid { 
  namespace QCD { 

class NerscHmcRunner {
public:

  enum StartType_t { ColdStart, HotStart, TepidStart, CheckpointStart };

  ActionSet<LatticeGaugeField> TheAction;

  GridCartesian         * UGrid   ;
  GridCartesian         * FGrid   ;
  GridRedBlackCartesian * UrbGrid ;
  GridRedBlackCartesian * FrbGrid ;

  void BuildTheAction (int argc, char **argv)
  {
    const int Ls = 8;

    UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
    UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  
    FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
    FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

    // temporarily need a gauge field
    LatticeGaugeField  U(UGrid);

    // Gauge action
    WilsonGaugeActionR Waction(5.6);

    Real mass=0.04;
    Real pv  =1.0;
    RealD M5=1.5;
    DomainWallFermionR DenOp(U,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);
    DomainWallFermionR NumOp(U,*FGrid,*FrbGrid,*UGrid,*UrbGrid,pv,M5);
  
    ConjugateGradient<LatticeFermion>  CG(1.0e-8,10000);
    TwoFlavourEvenOddRatioPseudoFermionAction<WilsonImplR> Nf2(NumOp, DenOp,CG,CG);
  
    //Collect actions
    ActionLevel<LatticeGaugeField> Level1;
    Level1.push_back(&Nf2);
    Level1.push_back(&Waction);
    TheAction.push_back(Level1);

    Run(argc,argv);
  };
  
  void Run (int argc, char  **argv){

    StartType_t StartType = HotStart;

    std::string arg;

    if( GridCmdOptionExists(argv,argv+argc,"--StartType") ){
      arg = GridCmdOptionPayload(argv,argv+argc,"--StartType");
      if ( arg == "HotStart" ) { StartType = HotStart; }
      else if ( arg == "ColdStart" ) { StartType = ColdStart; }
      else if ( arg == "TepidStart" ) { StartType = TepidStart; }
      else if ( arg == "CheckpointStart" ) { StartType = CheckpointStart; }
      else assert(0);
    }

    int StartTraj = 0;
    if( GridCmdOptionExists(argv,argv+argc,"--StartTrajectory") ){
      arg= GridCmdOptionPayload(argv,argv+argc,"--StartTrajectory");
      std::vector<int> ivec(0);
      GridCmdOptionIntVector(arg,ivec);
      StartTraj = ivec[0];
    }    

    int NumTraj = 1;
    if( GridCmdOptionExists(argv,argv+argc,"--Trajectories") ){
      arg= GridCmdOptionPayload(argv,argv+argc,"--Trajectories");
      std::vector<int> ivec(0);
      GridCmdOptionIntVector(arg,ivec);
      NumTraj = ivec[0];
    }

    // Create integrator
    typedef MinimumNorm2<LatticeGaugeField>  IntegratorType;// change here to change the algorithm
    IntegratorParameters MDpar(20);
    IntegratorType MDynamics(UGrid,MDpar, TheAction);

    // Checkpoint strategy
    NerscHmcCheckpointer<LatticeGaugeField> Checkpoint(std::string("ckpoint_lat"),std::string("ckpoint_rng"),1);
    PlaquetteLogger<LatticeGaugeField> PlaqLog(std::string("plaq"));

    HMCparameters HMCpar;
    HMCpar.StartTrajectory = StartTraj;
    HMCpar.Trajectories    = NumTraj;
    
    GridSerialRNG    sRNG;
    GridParallelRNG  pRNG(UGrid);
    LatticeGaugeField  U(UGrid);

    std::vector<int> SerSeed({1,2,3,4,5});
    std::vector<int> ParSeed({6,7,8,9,10});

    if ( StartType == HotStart ) {
      // Hot start
      HMCpar.NoMetropolisUntil =0;
      HMCpar.MetropolisTest = true;
      sRNG.SeedFixedIntegers(SerSeed);
      pRNG.SeedFixedIntegers(ParSeed);
      SU3::HotConfiguration(pRNG, U);
    } else if ( StartType == ColdStart ) { 
      // Cold start
      HMCpar.NoMetropolisUntil =0;
      HMCpar.MetropolisTest = true;
      sRNG.SeedFixedIntegers(SerSeed);
      pRNG.SeedFixedIntegers(ParSeed);
      SU3::ColdConfiguration(pRNG, U);
    } else if ( StartType == TepidStart ) {       
      // Tepid start
      HMCpar.NoMetropolisUntil =0;
      HMCpar.MetropolisTest = true;
      sRNG.SeedFixedIntegers(SerSeed);
      pRNG.SeedFixedIntegers(ParSeed);
      SU3::TepidConfiguration(pRNG, U);
    } else if ( StartType == CheckpointStart ) { 
      HMCpar.NoMetropolisUntil =0;
      HMCpar.MetropolisTest = true;
      // CheckpointRestart
      Checkpoint.CheckpointRestore(StartTraj, U, sRNG, pRNG);
    }

    HybridMonteCarlo<LatticeGaugeField,IntegratorType>  HMC(HMCpar, MDynamics,sRNG,pRNG,U);
    HMC.AddObservable(&Checkpoint);
    HMC.AddObservable(&PlaqLog);
    
    // Run it
    HMC.evolve();
    
  }
  
};

}}

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  int threads = GridThread::GetThreads();
  std::cout<<GridLogMessage << "Grid is setup to use "<<threads<<" threads"<<std::endl;

  NerscHmcRunner TheHMC;
  
  TheHMC.BuildTheAction(argc,argv);

}
