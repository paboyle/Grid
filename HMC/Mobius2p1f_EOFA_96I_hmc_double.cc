/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./tests/Test_hmc_EODWFRatio.cc

Copyright (C) 2015-2016

Author: Peter Boyle <pabobyle@ph.ed.ac.uk>
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



int main(int argc, char **argv) {
  using namespace Grid;

  std::cout << " Grid Initialise "<<std::endl;
  
  Grid_init(&argc, &argv);

  CartesianCommunicator::BarrierWorld();
  std::cout << GridLogMessage << " Clock skew check" <<std::endl;
  
  int threads = GridThread::GetThreads();

   // Typedefs to simplify notation
  typedef WilsonImplD FermionImplPolicy;
  typedef MobiusFermionD FermionAction;
  typedef MobiusEOFAFermionD FermionEOFAAction;
  typedef typename FermionAction::FermionField FermionField;

  typedef WilsonImplF FermionImplPolicyF;
  typedef MobiusFermionF FermionActionF;
  typedef MobiusEOFAFermionF FermionEOFAActionF;
  typedef typename FermionActionF::FermionField FermionFieldF;

  typedef Grid::XmlReader       Serialiser;

  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  IntegratorParameters MD;
  //  typedef GenericHMCRunner<LeapFrog> HMCWrapper;
  //  MD.name    = std::string("Leap Frog");
  typedef GenericHMCRunner<ForceGradient> HMCWrapper;
  MD.name    = std::string("Force Gradient");
  //  typedef GenericHMCRunner<MinimumNorm2> HMCWrapper;
  //  MD.name    = std::string("MinimumNorm2");
  // TrajL = 2
  // 4/2 => 0.6 dH
  // 3/3 => 0.8 dH .. depth 3, slower
  //MD.MDsteps =  4;
  MD.MDsteps =  8;
  MD.trajL   = 0.5;

  HMCparameters HMCparams;
  HMCparams.StartTrajectory  = 1077;
  HMCparams.Trajectories     = 20;
  HMCparams.NoMetropolisUntil=  0;
  // "[HotStart, ColdStart, TepidStart, CheckpointStart]\n";
  HMCparams.StartingType     =std::string("ColdStart");
  //  HMCparams.StartingType     =std::string("CheckpointStart");
  HMCparams.MD = MD;
  HMCWrapper TheHMC(HMCparams);

  // Grid from the command line arguments --grid and --mpi
  TheHMC.Resources.AddFourDimGrid("gauge"); // use default simd lanes decomposition

  CheckpointerParameters CPparams;
  CPparams.config_prefix = "ckpoint_HMC_lat";
  CPparams.rng_prefix    = "ckpoint_HMC_rng";
  CPparams.saveInterval  = 1;
  CPparams.format        = "IEEE64BIG";
  TheHMC.Resources.LoadNerscCheckpointer(CPparams);
  std::cout << "loaded NERSC checpointer"<<std::endl;
  RNGModuleParameters RNGpar;
  RNGpar.serial_seeds = "1 2 3 4 5 6 7 8 9 10";
  RNGpar.parallel_seeds = "6 7 8 9 10";
  TheHMC.Resources.SetRNGSeeds(RNGpar);

  // Construct observables
  // here there is too much indirection
  typedef PlaquetteMod<HMCWrapper::ImplPolicy> PlaqObs;
  TheHMC.Resources.AddObservable<PlaqObs>();
  //////////////////////////////////////////////

  const int Ls      = 12;
  RealD M5  = 1.8;
  RealD b   = 1.5;
  RealD c   = 0.5;
  RealD beta         = 2.13;
  //  Real light_mass   = 5.4e-4;
  Real light_mass     = 7.8e-4;
  //  Real light_mass     = 7.8e-3;
  Real strange_mass = 0.0362;
  Real pv_mass      = 1.0;
  std::vector<Real> hasenbusch({ 0.005, 0.0145, 0.045, 0.108, 0.25, 0.35 , 0.51, 0.6, 0.8 }); // Updated
  //std::vector<Real> hasenbusch({ 0.0145, 0.045, 0.108, 0.25, 0.35 , 0.51, 0.6, 0.8 }); // Updated

  auto GridPtr   = TheHMC.Resources.GetCartesian();
  auto GridRBPtr = TheHMC.Resources.GetRBCartesian();

  typedef SchurDiagMooeeOperator<FermionAction ,FermionField > LinearOperatorD;
  typedef SchurDiagMooeeOperator<FermionEOFAAction ,FermionField > LinearOperatorEOFAD;

  ////////////////////////////////////////////////////////////////
  // Domain decomposed
  ////////////////////////////////////////////////////////////////
  Coordinate latt4  = GridPtr->GlobalDimensions();
  Coordinate mpi    = GridPtr->ProcessorGrid();
  Coordinate shm;

  GlobalSharedMemory::GetShmDims(mpi,shm);

  //////////////////////////
  // Fermion Grids
  //////////////////////////
  auto FGrid     = SpaceTimeGrid::makeFiveDimGrid(Ls,GridPtr);
  auto FrbGrid   = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,GridPtr);

  IwasakiGaugeActionR GaugeAction(beta);

  // temporarily need a gauge field
  LatticeGaugeFieldD  U(GridPtr); U=Zero();

  std::cout << GridLogMessage << " Running the HMC "<< std::endl;
  TheHMC.ReadCommandLine(argc,argv);  // params on CML or from param file
  TheHMC.initializeGaugeFieldAndRNGs(U);
  std::cout << "loaded NERSC gauge field"<<std::endl;

  // These lines are unecessary if BC are all periodic
  std::vector<Complex> boundary = {1,1,1,-1};
  FermionAction::ImplParams Params(boundary);

  //  double StoppingCondition = 1e-14;
  //  double MDStoppingCondition = 1e-9;
  double StoppingCondition = 1e-14;
  double MDStoppingCondition = 1e-9;
  double MDStoppingConditionLoose = 1e-9;
  double MDStoppingConditionStrange = 1e-9;
  double MaxCGIterations = 50000;
  ConjugateGradient<FermionField>  CG(StoppingCondition,MaxCGIterations);
  ConjugateGradient<FermionField>  MDCG(MDStoppingCondition,MaxCGIterations);

  ////////////////////////////////////
  // Collect actions
  ////////////////////////////////////
  ActionLevel<HMCWrapper::Field> Level1(1);
  ActionLevel<HMCWrapper::Field> Level2(2);
  ActionLevel<HMCWrapper::Field> Level3(4);

  ////////////////////////////////////
  // Strange action
  ////////////////////////////////////
  FermionAction StrangeOp (U,*FGrid,*FrbGrid,*GridPtr,*GridRBPtr,strange_mass,M5,b,c, Params);
  FermionAction StrangePauliVillarsOp(U,*FGrid,*FrbGrid,*GridPtr,*GridRBPtr,pv_mass,  M5,b,c, Params);

  // Probably dominates the force - back to EOFA.
  OneFlavourRationalParams SFRp;
  SFRp.lo       = 0.8;
  SFRp.hi       = 30.0;
  SFRp.MaxIter  = 10000;
  SFRp.tolerance= 1.0e-12;
  SFRp.mdtolerance= 1.0e-9;
  SFRp.degree   = 10;
  SFRp.precision= 50;
  
  MobiusEOFAFermionD Strange_Op_L (U , *FGrid , *FrbGrid , *GridPtr , *GridRBPtr , strange_mass, strange_mass, pv_mass, 0.0, -1, M5, b, c);
  MobiusEOFAFermionD Strange_Op_R (U , *FGrid , *FrbGrid , *GridPtr , *GridRBPtr , pv_mass, strange_mass,      pv_mass, -1.0, 1, M5, b, c);
  ConjugateGradient<FermionField>      ActionCG(StoppingCondition,MaxCGIterations);
  ConjugateGradient<FermionField>  DerivativeCG(MDStoppingCondition,MaxCGIterations);
  LinearOperatorEOFAD Strange_LinOp_L (Strange_Op_L);
  LinearOperatorEOFAD Strange_LinOp_R (Strange_Op_R);

  ExactOneFlavourRatioPseudoFermionAction<FermionImplPolicy> 
    EOFA(Strange_Op_L, Strange_Op_R, 
	 ActionCG, 
	 ActionCG, ActionCG,
	 DerivativeCG, DerivativeCG,
	 SFRp, true);
  Level2.push_back(&EOFA);

  ////////////////////////////////////
  // up down action
  ////////////////////////////////////
  std::vector<Real> light_den;
  std::vector<Real> light_num;

  int n_hasenbusch = hasenbusch.size();
  light_den.push_back(light_mass); 
  for(int h=0;h<n_hasenbusch;h++){
    light_den.push_back(hasenbusch[h]);
  }

  for(int h=0;h<n_hasenbusch;h++){
    light_num.push_back(hasenbusch[h]);
  }
  light_num.push_back(pv_mass);

  std::vector<FermionAction *> Numerators;
  std::vector<FermionAction *> Denominators;
  std::vector<TwoFlavourEvenOddRatioPseudoFermionAction<FermionImplPolicy> *> Quotients;
  
  std::vector<OneFlavourEvenOddRatioRationalPseudoFermionAction<FermionImplPolicy> *> Bdys;

  typedef SchurDiagMooeeOperator<FermionAction ,FermionField > LinearOperatorD;
  std::vector<LinearOperatorD *> LinOpD;
  
  for(int h=0;h<n_hasenbusch+1;h++){
    std::cout << GridLogMessage
	      << " 2f quotient Action ";
    std::cout << "det D("<<light_den[h]<<")";
    std::cout << "/ det D("<<light_num[h]<<")";
    std::cout << std::endl;

    FermionAction::ImplParams ParamsNum(boundary);
    FermionAction::ImplParams ParamsDen(boundary);
    
    Numerators.push_back  (new FermionAction(U,*FGrid,*FrbGrid,*GridPtr,*GridRBPtr,light_num[h],M5,b,c, ParamsNum));
    Denominators.push_back(new FermionAction(U,*FGrid,*FrbGrid,*GridPtr,*GridRBPtr,light_den[h],M5,b,c, ParamsDen));

    LinOpD.push_back(new LinearOperatorD(*Denominators[h]));

    double conv  = MDStoppingCondition;
    if (h<3) conv= MDStoppingConditionLoose; // Relax on first two hasenbusch factors
    
    Quotients.push_back (new TwoFlavourEvenOddRatioPseudoFermionAction<FermionImplPolicy>(*Numerators[h],*Denominators[h],MDCG,CG,CG));
  }
  int nquo=Quotients.size();
  for(int h=0;h<nquo;h++){
    Level1.push_back(Quotients[h]);
  }

  /////////////////////////////////////////////////////////////
  // Gauge action
  /////////////////////////////////////////////////////////////
  Level3.push_back(&GaugeAction);
  TheHMC.TheAction.push_back(Level1);
  TheHMC.TheAction.push_back(Level2);
  TheHMC.TheAction.push_back(Level3);
  std::cout << GridLogMessage << " Action complete "<< std::endl;
  /////////////////////////////////////////////////////////////

  TheHMC.Run();  // no smearing

  Grid_finalize();
} // main



