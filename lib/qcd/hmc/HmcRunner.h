    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/hmc/HmcRunner.h

    Copyright (C) 2015

Author: paboyle <paboyle@ph.ed.ac.uk>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

    See the full license in the file "LICENSE" in the top level distribution directory
    *************************************************************************************/
    /*  END LEGAL */
#ifndef HMC_RUNNER
#define HMC_RUNNER

namespace Grid{
  namespace QCD{


template<class Gimpl>
class NerscHmcRunnerTemplate {
public:

  INHERIT_GIMPL_TYPES(Gimpl);

  enum StartType_t { ColdStart, HotStart, TepidStart, CheckpointStart };

  ActionSet<GaugeField> TheAction;

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
    typedef MinimumNorm2<GaugeField>  IntegratorType;// change here to change the algorithm
    IntegratorParameters MDpar(20);
    IntegratorType MDynamics(UGrid,MDpar, TheAction);

    // Checkpoint strategy
    NerscHmcCheckpointer<Gimpl> Checkpoint(std::string("ckpoint_lat"),std::string("ckpoint_rng"),1);
    PlaquetteLogger<Gimpl>      PlaqLog(std::string("plaq"));

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
      HMCpar.NoMetropolisUntil =10;
      HMCpar.MetropolisTest = true;
      sRNG.SeedFixedIntegers(SerSeed);
      pRNG.SeedFixedIntegers(ParSeed);
      SU3::HotConfiguration(pRNG, U);
    } else if ( StartType == ColdStart ) { 
      // Cold start
      HMCpar.NoMetropolisUntil =10;
      HMCpar.MetropolisTest = true;
      sRNG.SeedFixedIntegers(SerSeed);
      pRNG.SeedFixedIntegers(ParSeed);
      SU3::ColdConfiguration(pRNG, U);
    } else if ( StartType == TepidStart ) {       
      // Tepid start
      HMCpar.NoMetropolisUntil =10;
      HMCpar.MetropolisTest = true;
      sRNG.SeedFixedIntegers(SerSeed);
      pRNG.SeedFixedIntegers(ParSeed);
      SU3::TepidConfiguration(pRNG, U);
    } else if ( StartType == CheckpointStart ) { 
      HMCpar.NoMetropolisUntil =10;
      HMCpar.MetropolisTest = true;
      // CheckpointRestart
      Checkpoint.CheckpointRestore(StartTraj, U, sRNG, pRNG);
    }

    HybridMonteCarlo<GaugeField,IntegratorType>  HMC(HMCpar, MDynamics,sRNG,pRNG,U);
    HMC.AddObservable(&Checkpoint);
    HMC.AddObservable(&PlaqLog);
    
    // Run it
    HMC.evolve();
    
  }
  
};

 typedef NerscHmcRunnerTemplate<PeriodicGimplR> NerscHmcRunner;
 typedef NerscHmcRunnerTemplate<PeriodicGimplF> NerscHmcRunnerF;
 typedef NerscHmcRunnerTemplate<PeriodicGimplD> NerscHmcRunnerD;

 typedef NerscHmcRunnerTemplate<PeriodicGimplR> PeriodicNerscHmcRunner;
 typedef NerscHmcRunnerTemplate<PeriodicGimplF> PeriodicNerscHmcRunnerF;
 typedef NerscHmcRunnerTemplate<PeriodicGimplD> PeriodicNerscHmcRunnerD;

 typedef NerscHmcRunnerTemplate<ConjugateGimplR> ConjugateNerscHmcRunner;
 typedef NerscHmcRunnerTemplate<ConjugateGimplF> ConjugateNerscHmcRunnerF;
 typedef NerscHmcRunnerTemplate<ConjugateGimplD> ConjugateNerscHmcRunnerD;

}}
#endif
