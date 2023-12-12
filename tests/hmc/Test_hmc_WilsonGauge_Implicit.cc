/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./tests/Test_hmc_WilsonFermionGauge.cc

Copyright (C) 2015

Author: Peter Boyle <pabobyle@ph.ed.ac.uk>
Author: neo <cossu@post.kek.jp>
Author: Guido Cossu <guido.cossu@ed.ac.uk>

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

See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
/*  END LEGAL */
#include <Grid/Grid.h>

#define USE_OBC
#define DO_IMPLICIT


int main(int argc, char **argv) 
{
  using namespace Grid;

  Grid_init(&argc, &argv);
  GridLogLayout();
  
  std::string arg;
 
  HMCparameters HMCparams;
#if 1
  {
    XmlReader  HMCrd("HMCparameters.xml");
    read(HMCrd,"HMCparameters",HMCparams);
  }
#else
//IntegratorParameters MD;
  std::vector<int> steps(0);
  if( GridCmdOptionExists(argv,argv+argc,"--MDsteps") ){
    arg= GridCmdOptionPayload(argv,argv+argc,"--MDsteps");
    GridCmdOptionIntVector(arg,steps);
    assert(steps.size()==1);
  }
  MD.trajL   = 0.001*std::sqrt(2.);
  MD.MDsteps = 1;
  if (steps.size()>0) MD.MDsteps = steps[0];
  if( GridCmdOptionExists(argv,argv+argc,"--trajL") ){
    arg= GridCmdOptionPayload(argv,argv+argc,"--trajL");
    std::vector<int> traj(0);
    GridCmdOptionIntVector(arg,traj);
    assert(traj.size()==1);
    MD.trajL *= double(traj[0]);
  }
  MD.RMHMCTol=1e-8;
  MD.RMHMCCGTol=1e-8;
  std::cout << "RMHMCTol= "<<  MD.RMHMCTol<<" RMHMCCGTol= "<<MD.RMHMCCGTol<<std::endl;

  HMCparameters HMCparams;
  HMCparams.StartTrajectory  = 0;
  HMCparams.Trajectories     = 1;
  HMCparams.NoMetropolisUntil=  100;
  // "[HotStart, ColdStart, TepidStart, CheckpointStart]\n";
  HMCparams.StartingType     =std::string("ColdStart");
  HMCparams.Kappa=0.01; //checking against trivial. Pathetic.
  HMCparams.MD = MD;
#endif



   // Typedefs to simplify notation
#ifdef DO_IMPLICIT
  typedef GenericHMCRunner<ImplicitMinimumNorm2> HMCWrapper;  // Uses the default minimum norm
//  typedef GenericHMCRunner<ImplicitCampostrini> HMCWrapper;  // 4th order
  HMCparams.MD.name    = std::string("ImplicitMinimumNorm2");
#else
  typedef GenericHMCRunner<MinimumNorm2> HMCWrapper;  // Uses the default minimum norm
  HMCparams.MD.name    = std::string("MinimumNorm2");
#endif



  // Possibile to create the module by hand 
  // hardcoding parameters or using a Reader


  // Checkpointer definition
  CheckpointerParameters CPparams;  
  CPparams.config_prefix = "ckpoint_lat";
  CPparams.rng_prefix = "ckpoint_rng";
  CPparams.saveInterval = 1;
  CPparams.format = "IEEE64BIG";
  
  HMCWrapper TheHMC(HMCparams);
  // Grid from the command line
  TheHMC.Resources.AddFourDimGrid("gauge");
  TheHMC.Resources.LoadNerscCheckpointer(CPparams);

  RNGModuleParameters RNGpar;
  RNGpar.serial_seeds = "1 2 3 4 5";
  RNGpar.parallel_seeds = "6 7 8 9 10";
  TheHMC.Resources.SetRNGSeeds(RNGpar);

  // Construct observables
  // here there is too much indirection 
  typedef PlaquetteMod<HMCWrapper::ImplPolicy> PlaqObs;
  typedef TopologicalChargeMod<HMCWrapper::ImplPolicy> QObs;
  TheHMC.Resources.AddObservable<PlaqObs>();
  TopologyObsParameters TopParams;
  TopParams.interval = 1;
  TopParams.do_smearing = true;
//  TopParams.Smearing.steps = 1600;
//  TopParams.Smearing.step_size = 0.01;
  TopParams.Smearing.init_step_size = 0.01;
  TopParams.Smearing.meas_interval = 10;
  TopParams.Smearing.maxTau = 16.0; 
//  TheHMC.Resources.AddObservable<QObs>(TopParams);
  //////////////////////////////////////////////

  /////////////////////////////////////////////////////////////
  // Collect actions, here use more encapsulation
  // need wrappers of the fermionic classes 
  // that have a complex construction
  // standard

  RealD beta = 6.4;
  std::cout << "Wilson Gauge beta= " <<beta <<std::endl;
#ifndef USE_OBC
  WilsonGaugeActionR Waction(beta);
#else
  std::vector<Complex> boundaryG = {1,1,1,0};
  WilsonGaugeActionR::ImplParams ParamsG(boundaryG);
  WilsonGaugeActionR Waction(beta,ParamsG);
  std::cout << "boundaryG = " <<boundaryG  <<std::endl;
#endif

  
  ActionLevel<HMCWrapper::Field> Level1(1);
  Level1.push_back(&Waction);
  TheHMC.TheAction.push_back(Level1);

  TheHMC.ReadCommandLine(argc, argv); // these can be parameters from file
  std::cout << "trajL= " <<TheHMC.Parameters.MD.trajL <<" steps= "<<TheHMC.Parameters.MD.MDsteps << " integrator= "<<TheHMC.Parameters.MD.name<<std::endl;

  NoSmearing<HMCWrapper::ImplPolicy> S;
#ifndef DO_IMPLICIT
  TrivialMetric<HMCWrapper::ImplPolicy::Field> Mtr;
#else
// g_x3_2
    LaplacianRatParams gpar(2),mpar(2);
    gpar.offset = 1.;
    gpar.a0[0] = 500.;
    gpar.a1[0] = 0.;
    gpar.b0[0] = 0.25;
    gpar.b1[0] = 1.;
    gpar.a0[1] = -500.;
    gpar.a1[1] = 0.;
    gpar.b0[1] = 0.36;
    gpar.b1[1] = 1.2;
    gpar.b2=1.;

    mpar.offset = 1.;
    mpar.a0[0] =  -0.850891906532;
    mpar.a1[0] = -1.54707654538;
    mpar. b0[0] = 2.85557166137;
    mpar. b1[0] = 5.74194794773;
    mpar.a0[1] = -13.5120056831218384729709214298;
    mpar.a1[1] = 1.54707654538396877086370295729;
    mpar.b0[1] = 19.2921090880640520026645390317;
    mpar.b1[1] = -3.54194794773029020262811172870;
    mpar.b2=1.;
    for(int i=0;i<2;i++){
       gpar.a1[i] *=16.;
       gpar.b1[i] *=16.;
       mpar.a1[i] *=16.;
       mpar.b1[i] *=16.;
    }
    gpar.b2 *= 16.*16.;
    mpar.b2 *= 16.*16.;

    ConjugateGradient<LatticeGaugeField> CG(1.0e-8,10000);
    LaplacianParams LapPar(0.0001, 1.0, 10000, 1e-8, 12, 64);

    std::cout << GridLogMessage << "LaplacianRat " << std::endl;

    gpar.tolerance=HMCparams.MD.RMHMCCGTol;
    mpar.tolerance=HMCparams.MD.RMHMCCGTol;

    std::cout << GridLogMessage << "gpar offset= " << gpar.offset <<std::endl;
    std::cout << GridLogMessage << " a0= " << gpar.a0 <<std::endl;
    std::cout << GridLogMessage << " a1= " << gpar.a1 <<std::endl;
    std::cout << GridLogMessage << " b0= " << gpar.b0 <<std::endl;
    std::cout << GridLogMessage << " b1= " << gpar.b1 <<std::endl;
    std::cout << GridLogMessage << " b2= " << gpar.b2 <<std::endl ;;

    std::cout << GridLogMessage << "mpar offset= " << mpar.offset <<std::endl;
    std::cout << GridLogMessage << " a0= " << mpar.a0 <<std::endl;
    std::cout << GridLogMessage << " a1= " << mpar.a1 <<std::endl;
    std::cout << GridLogMessage << " b0= " << mpar.b0 <<std::endl;
    std::cout << GridLogMessage << " b1= " << mpar.b1 <<std::endl;
    std::cout << GridLogMessage << " b2= " << mpar.b2 <<std::endl;
//  Assumes PeriodicGimplR or D at the moment
    Coordinate latt  = GridDefaultLatt();
    Coordinate mpi   = GridDefaultMpi();
    auto UGrid = TheHMC.Resources.GetCartesian("gauge");
    Coordinate simdF = GridDefaultSimd(Nd,vComplexF::Nsimd());
    auto UGrid_f   = SpaceTimeGrid::makeFourDimGrid(latt,simdF,mpi);
    std::cout << GridLogMessage << " UGrid= " << UGrid <<std::endl;
    std::cout << GridLogMessage << " UGrid_f= " << UGrid_f <<std::endl;

    LaplacianAdjointRat<HMCWrapper::ImplPolicy, PeriodicGimplF> Mtr(UGrid, UGrid_f,CG, gpar, mpar);
#endif
 
  {
    XmlWriter HMCwr("HMCparameters.xml.out");
    write(HMCwr,"HMCparameters",TheHMC.Parameters);
  }

  TheHMC.Run(S,Mtr);  // no smearing

  Grid_finalize();

} // main
