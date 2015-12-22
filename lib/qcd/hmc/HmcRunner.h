#ifndef HMC_RUNNER
#define HMC_RUNNER

namespace Grid{
  namespace QCD{


class NerscHmcRunner {
public:

  enum StartType_t { ColdStart, HotStart, TepidStart, CheckpointStart };

  ActionSet<LatticeGaugeField> TheAction;

  GridCartesian         * UGrid   ;
  GridCartesian         * FGrid   ;
  GridRedBlackCartesian * UrbGrid ;
  GridRedBlackCartesian * FrbGrid ;

  virtual void BuildTheAction (int argc, char **argv) = 0;

  
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
#endif
