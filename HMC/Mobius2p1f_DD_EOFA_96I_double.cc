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

  Grid_init(&argc, &argv);

  CartesianCommunicator::BarrierWorld();
  std::cout << GridLogMessage << " Clock skew check" <<std::endl;
  
  int threads = GridThread::GetThreads();

   // Typedefs to simplify notation
  typedef WilsonImplD FermionImplPolicy;
  typedef MobiusFermionD FermionAction;
  typedef MobiusEOFAFermionD FermionEOFAAction;
  typedef typename FermionAction::FermionField FermionField;

  typedef Grid::XmlReader       Serialiser;

  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  IntegratorParameters MD;
  //  typedef GenericHMCRunner<LeapFrog> HMCWrapper;
  //  MD.name    = std::string("Leap Frog");
  typedef GenericHMCRunner<ForceGradient> HMCWrapper;
  MD.name    = std::string("Force Gradient");
  //typedef GenericHMCRunner<MinimumNorm2> HMCWrapper;
  // MD.name    = std::string("MinimumNorm2");
  // TrajL = 2
  // 4/2 => 0.6 dH
  // 3/3 => 0.8 dH .. depth 3, slower
  //MD.MDsteps =  4;
  MD.MDsteps =  3;
  MD.trajL   = 0.5;

  HMCparameters HMCparams;
  HMCparams.StartTrajectory  = 1077;
  HMCparams.Trajectories     = 1;
  HMCparams.NoMetropolisUntil=  0;
  // "[HotStart, ColdStart, TepidStart, CheckpointStart]\n";
  //  HMCparams.StartingType     =std::string("ColdStart");
  HMCparams.StartingType     =std::string("CheckpointStart");
  HMCparams.MD = MD;
  HMCWrapper TheHMC(HMCparams);

  // Grid from the command line arguments --grid and --mpi
  TheHMC.Resources.AddFourDimGrid("gauge"); // use default simd lanes decomposition

  CheckpointerParameters CPparams;
  CPparams.config_prefix = "ckpoint_DDHMC_lat";
  CPparams.rng_prefix    = "ckpoint_DDHMC_rng";
  CPparams.saveInterval  = 1;
  CPparams.format        = "IEEE64BIG";
  TheHMC.Resources.LoadNerscCheckpointer(CPparams);
  std::cout << "loaded NERSC checpointer"<<std::endl;
  RNGModuleParameters RNGpar;
  RNGpar.serial_seeds = "1 2 3 4 5";
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
  Real beta         = 2.13;
  //  Real light_mass   = 5.4e-4;
  Real light_mass     = 7.8e-4;
  Real light_mass_dir = 0.01;
  Real strange_mass = 0.0362;
  Real pv_mass      = 1.0;
  std::vector<Real> hasenbusch({ 0.01, 0.045, 0.108, 0.25, 0.51 , pv_mass });
  //  std::vector<Real> hasenbusch({ light_mass, 0.01, 0.045, 0.108, 0.25, 0.51 , pv_mass });
  //  std::vector<Real> hasenbusch({ light_mass, 0.005, 0.0145, 0.045, 0.108, 0.25, 0.51 , pv_mass }); // Updated
  //  std::vector<Real> hasenbusch({ light_mass, 0.0145, 0.045, 0.108, 0.25, 0.51 , 0.75 , pv_mass });

  int SP_iters=9000;
  
  RationalActionParams OFRp; // Up/down
  OFRp.lo       = 6.0e-5;
  OFRp.hi       = 90.0;
  OFRp.inv_pow  = 2;
  OFRp.MaxIter  = SP_iters; // get most shifts by 2000, stop sharing space
  OFRp.action_tolerance= 1.0e-8;
  OFRp.action_degree   = 18;
  OFRp.md_tolerance= 1.0e-7;
  OFRp.md_degree   = 14;
  //  OFRp.degree   = 20; converges
  //  OFRp.degree   = 16;
  OFRp.precision= 80;
  OFRp.BoundsCheckFreq=0;
  std::vector<RealD> ActionTolByPole({
      //      1.0e-8,1.0e-8,1.0e-8,1.0e-8,
      3.0e-7,1.0e-7,1.0e-8,1.0e-8,
      1.0e-8,1.0e-8,1.0e-8,1.0e-8,
      1.0e-8,1.0e-8,1.0e-8,1.0e-8,
      1.0e-8,1.0e-8,1.0e-8,1.0e-8,
      1.0e-8,1.0e-8
    });
  std::vector<RealD> MDTolByPole({
      //      1.6e-5,5.0e-6,1.0e-6,3.0e-7, // soften convergence more more
      //      1.0e-6,3.0e-7,1.0e-7,1.0e-7,
      1.0e-5,1.0e-6,1.0e-7,1.0e-7, // soften convergence
      1.0e-8,1.0e-8,1.0e-8,1.0e-8,
      1.0e-8,1.0e-8,1.0e-8,1.0e-8,
      1.0e-8,1.0e-8
    });

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
  
  Coordinate CommDim(Nd);
  for(int d=0;d<Nd;d++) CommDim[d]= (mpi[d]/shm[d])>1 ? 1 : 0;

  Coordinate NonDirichlet(Nd+1,0);
  Coordinate Dirichlet(Nd+1,0);
  Dirichlet[1] = CommDim[0]*latt4[0]/mpi[0] * shm[0];
  Dirichlet[2] = CommDim[1]*latt4[1]/mpi[1] * shm[1];
  Dirichlet[3] = CommDim[2]*latt4[2]/mpi[2] * shm[2];
  Dirichlet[4] = CommDim[3]*latt4[3]/mpi[3] * shm[3];
  //Dirichlet[1] = 0;
  //Dirichlet[2] = 0;
  //Dirichlet[3] = 0;

  // 
  Coordinate Block4(Nd);
  Block4[0] = Dirichlet[1];
  Block4[1] = Dirichlet[2];
  Block4[2] = Dirichlet[3];
  Block4[3] = Dirichlet[4];

  int Width=4;
  TheHMC.Resources.SetMomentumFilter(new DDHMCFilter<WilsonImplD::Field>(Block4,Width));

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
  FermionAction::ImplParams ParamsDir(boundary);

  Params.dirichlet=NonDirichlet;
  ParamsDir.dirichlet=Dirichlet;
  ParamsDir.partialDirichlet=0;
  std::cout << GridLogMessage<< "Partial Dirichlet depth is "<<dwf_compressor_depth<<std::endl;

  //  double StoppingCondition = 1e-14;
  //  double MDStoppingCondition = 1e-9;
  double StoppingCondition = 1e-8;
  double MDStoppingCondition = 1e-8;
  double MDStoppingConditionLoose = 1e-8;
  double MDStoppingConditionStrange = 1e-8;
  double MaxCGIterations = 300000;
  ConjugateGradient<FermionField>  CG(StoppingCondition,MaxCGIterations);
  ConjugateGradient<FermionField>  MDCG(MDStoppingCondition,MaxCGIterations);

  ////////////////////////////////////
  // Collect actions
  ////////////////////////////////////
  ActionLevel<HMCWrapper::Field> Level1(1);
  ActionLevel<HMCWrapper::Field> Level2(3);
  ActionLevel<HMCWrapper::Field> Level3(15);

  ////////////////////////////////////
  // Strange action
  ////////////////////////////////////
  FermionAction StrangeOp (U,*FGrid,*FrbGrid,*GridPtr,*GridRBPtr,strange_mass,M5,b,c, Params);
  FermionAction StrangePauliVillarsOp(U,*FGrid,*FrbGrid,*GridPtr,*GridRBPtr,pv_mass,  M5,b,c, Params);

  // Probably dominates the force - back to EOFA.
  OneFlavourRationalParams SFRp;
  SFRp.lo       = 0.1;
  SFRp.hi       = 25.0;
  SFRp.MaxIter  = 10000;
  SFRp.tolerance= 1.0e-8;
  SFRp.mdtolerance= 2.0e-6;
  SFRp.degree   = 12;
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
  std::vector<int> dirichlet_den;
  std::vector<int> dirichlet_num;

  int n_hasenbusch = hasenbusch.size();
  light_den.push_back(light_mass);  dirichlet_den.push_back(0);
  for(int h=0;h<n_hasenbusch;h++){
    light_den.push_back(hasenbusch[h]); dirichlet_den.push_back(1);
  }

  for(int h=0;h<n_hasenbusch;h++){
    light_num.push_back(hasenbusch[h]); dirichlet_num.push_back(1);
  }
  light_num.push_back(pv_mass);  dirichlet_num.push_back(0);

  std::vector<FermionAction *> Numerators;
  std::vector<FermionAction *> Denominators;
  std::vector<TwoFlavourEvenOddRatioPseudoFermionAction<FermionImplPolicy> *> Quotients;
  
  std::vector<GeneralEvenOddRatioRationalPseudoFermionAction<FermionImplPolicy> *> Bdys;

  typedef SchurDiagMooeeOperator<FermionAction ,FermionField > LinearOperatorD;
  std::vector<LinearOperatorD *> LinOpD;
  
  for(int h=0;h<n_hasenbusch+1;h++){
    std::cout << GridLogMessage
	      << " 2f quotient Action ";
    std::cout << "det D("<<light_den[h]<<")";
    if ( dirichlet_den[h] ) std::cout << "^dirichlet    ";
    std::cout << "/ det D("<<light_num[h]<<")";
    if ( dirichlet_num[h] ) std::cout << "^dirichlet    ";
    std::cout << std::endl;

    FermionAction::ImplParams ParamsNum(boundary);
    FermionAction::ImplParams ParamsDen(boundary);
    
    if ( dirichlet_num[h]==1) ParamsNum.dirichlet = Dirichlet;
    else                      ParamsNum.dirichlet = NonDirichlet;

    if ( dirichlet_den[h]==1) ParamsDen.dirichlet = Dirichlet;
    else                      ParamsDen.dirichlet = NonDirichlet;

    if ( dirichlet_num[h]==1) ParamsNum.partialDirichlet = 1;
    else                      ParamsNum.partialDirichlet = 0;

    if ( dirichlet_den[h]==1) ParamsDen.partialDirichlet = 1;
    else                      ParamsDen.partialDirichlet = 0;
    
    Numerators.push_back  (new FermionAction(U,*FGrid,*FrbGrid,*GridPtr,*GridRBPtr,light_num[h],M5,b,c, ParamsNum));
    Denominators.push_back(new FermionAction(U,*FGrid,*FrbGrid,*GridPtr,*GridRBPtr,light_den[h],M5,b,c, ParamsDen));

    LinOpD.push_back(new LinearOperatorD(*Denominators[h]));

    double conv  = MDStoppingCondition;
    if (h<3) conv= MDStoppingConditionLoose; // Relax on first two hasenbusch factors
    
    if(h!=0) {
      Quotients.push_back (new TwoFlavourEvenOddRatioPseudoFermionAction<FermionImplPolicy>(*Numerators[h],*Denominators[h],MDCG,CG));
    } else {
      Bdys.push_back( new GeneralEvenOddRatioRationalPseudoFermionAction<FermionImplPolicy>(*Numerators[h],*Denominators[h],OFRp));
      Bdys.push_back( new GeneralEvenOddRatioRationalPseudoFermionAction<FermionImplPolicy>(*Numerators[h],*Denominators[h],OFRp));
    }
  }
  for(int h=0;h<Bdys.size();h++){
    Bdys[h]->SetTolerances(ActionTolByPole,MDTolByPole);
  }
  int nquo=Quotients.size();
  Level1.push_back(Bdys[0]);
  Level1.push_back(Bdys[1]);
  Level2.push_back(Quotients[0]);
  for(int h=1;h<nquo-1;h++){
    Level2.push_back(Quotients[h]);
  }
  Level2.push_back(Quotients[nquo-1]);

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



