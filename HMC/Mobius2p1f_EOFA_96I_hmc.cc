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

NAMESPACE_BEGIN(Grid);

template<class FermionOperatorD, class FermionOperatorF, class SchurOperatorD, class  SchurOperatorF> 
  class MixedPrecisionConjugateGradientOperatorFunction : public OperatorFunction<typename FermionOperatorD::FermionField> {
  public:
    typedef typename FermionOperatorD::FermionField FieldD;
    typedef typename FermionOperatorF::FermionField FieldF;

    using OperatorFunction<FieldD>::operator();

    RealD   Tolerance;
    RealD   InnerTolerance; //Initial tolerance for inner CG. Defaults to Tolerance but can be changed
    Integer MaxInnerIterations;
    Integer MaxOuterIterations;
    GridBase* SinglePrecGrid4; //Grid for single-precision fields
    GridBase* SinglePrecGrid5; //Grid for single-precision fields
    RealD OuterLoopNormMult; //Stop the outer loop and move to a final double prec solve when the residual is OuterLoopNormMult * Tolerance

    FermionOperatorF &FermOpF;
    FermionOperatorD &FermOpD;;
    SchurOperatorF &LinOpF;
    SchurOperatorD &LinOpD;

    Integer TotalInnerIterations; //Number of inner CG iterations
    Integer TotalOuterIterations; //Number of restarts
    Integer TotalFinalStepIterations; //Number of CG iterations in final patch-up step

    MixedPrecisionConjugateGradientOperatorFunction(RealD tol, 
						    Integer maxinnerit, 
						    Integer maxouterit, 
						    GridBase* _sp_grid4, 
						    GridBase* _sp_grid5, 
						    FermionOperatorF &_FermOpF,
						    FermionOperatorD &_FermOpD,
						    SchurOperatorF   &_LinOpF,
						    SchurOperatorD   &_LinOpD): 
      LinOpF(_LinOpF),
      LinOpD(_LinOpD),
      FermOpF(_FermOpF),
      FermOpD(_FermOpD),
      Tolerance(tol), 
      InnerTolerance(tol), 
      MaxInnerIterations(maxinnerit), 
      MaxOuterIterations(maxouterit), 
      SinglePrecGrid4(_sp_grid4),
      SinglePrecGrid5(_sp_grid5),
      OuterLoopNormMult(100.) 
    { 
      /* Debugging instances of objects; references are stored
      std::cout << GridLogMessage << " Mixed precision CG wrapper LinOpF " <<std::hex<< &LinOpF<<std::dec <<std::endl;
      std::cout << GridLogMessage << " Mixed precision CG wrapper LinOpD " <<std::hex<< &LinOpD<<std::dec <<std::endl;
      std::cout << GridLogMessage << " Mixed precision CG wrapper FermOpF " <<std::hex<< &FermOpF<<std::dec <<std::endl;
      std::cout << GridLogMessage << " Mixed precision CG wrapper FermOpD " <<std::hex<< &FermOpD<<std::dec <<std::endl;
      */
    };

    void operator()(LinearOperatorBase<FieldD> &LinOpU, const FieldD &src, FieldD &psi) {

      std::cout << GridLogMessage << " Mixed precision CG wrapper operator() "<<std::endl;

      SchurOperatorD * SchurOpU = static_cast<SchurOperatorD *>(&LinOpU);
      
      //      std::cout << GridLogMessage << " Mixed precision CG wrapper operator() FermOpU " <<std::hex<< &(SchurOpU->_Mat)<<std::dec <<std::endl;
      //      std::cout << GridLogMessage << " Mixed precision CG wrapper operator() FermOpD " <<std::hex<< &(LinOpD._Mat) <<std::dec <<std::endl;
      // Assumption made in code to extract gauge field
      // We could avoid storing LinopD reference alltogether ?
      assert(&(SchurOpU->_Mat)==&(LinOpD._Mat));

      ////////////////////////////////////////////////////////////////////////////////////
      // Must snarf a single precision copy of the gauge field in Linop_d argument
      ////////////////////////////////////////////////////////////////////////////////////
      typedef typename FermionOperatorF::GaugeField GaugeFieldF;
      typedef typename FermionOperatorF::GaugeLinkField GaugeLinkFieldF;
      typedef typename FermionOperatorD::GaugeField GaugeFieldD;
      typedef typename FermionOperatorD::GaugeLinkField GaugeLinkFieldD;

      GridBase * GridPtrF = SinglePrecGrid4;
      GridBase * GridPtrD = FermOpD.Umu.Grid();
      GaugeFieldF     U_f  (GridPtrF);
      GaugeLinkFieldF Umu_f(GridPtrF);
      //      std::cout << " Dim gauge field "<<GridPtrF->Nd()<<std::endl; // 4d
      //      std::cout << " Dim gauge field "<<GridPtrD->Nd()<<std::endl; // 4d

      ////////////////////////////////////////////////////////////////////////////////////
      // Moving this to a Clone method of fermion operator would allow to duplicate the 
      // physics parameters and decrease gauge field copies
      ////////////////////////////////////////////////////////////////////////////////////
      GaugeLinkFieldD Umu_d(GridPtrD);
      for(int mu=0;mu<Nd*2;mu++){ 
	Umu_d = PeekIndex<LorentzIndex>(FermOpD.Umu, mu);
	precisionChange(Umu_f,Umu_d);
	PokeIndex<LorentzIndex>(FermOpF.Umu, Umu_f, mu);
      }
      pickCheckerboard(Even,FermOpF.UmuEven,FermOpF.Umu);
      pickCheckerboard(Odd ,FermOpF.UmuOdd ,FermOpF.Umu);

      ////////////////////////////////////////////////////////////////////////////////////
      // Make a mixed precision conjugate gradient
      ////////////////////////////////////////////////////////////////////////////////////
#if 1
      RealD delta=1.e-4;
      std::cout << GridLogMessage << "Calling reliable update Conjugate Gradient" <<std::endl;
      ConjugateGradientReliableUpdate<FieldD,FieldF> MPCG(Tolerance,MaxInnerIterations*MaxOuterIterations,delta,SinglePrecGrid5,LinOpF,LinOpD);
#else      
      std::cout << GridLogMessage << "Calling mixed precision Conjugate Gradient" <<std::endl;
      MixedPrecisionConjugateGradient<FieldD,FieldF> MPCG(Tolerance,MaxInnerIterations,MaxOuterIterations,SinglePrecGrid5,LinOpF,LinOpD);
#endif
      MPCG(src,psi);
    }
  };

NAMESPACE_END(Grid);


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

  typedef WilsonImplF FermionImplPolicyF;
  typedef MobiusFermionF FermionActionF;
  typedef MobiusEOFAFermionF FermionEOFAActionF;
  typedef typename FermionActionF::FermionField FermionFieldF;

  typedef WilsonImplD2 FermionImplPolicyD2;
  typedef MobiusFermionD2 FermionActionD2;
  typedef MobiusEOFAFermionD2 FermionEOFAActionD2;
  typedef typename FermionActionD2::FermionField FermionFieldD2;

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
  MD.MDsteps =  12;
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
  CPparams.config_prefix = "ckpoint_HMC_lat";
  CPparams.rng_prefix    = "ckpoint_HMC_rng";
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
  Real strange_mass = 0.0362;
  Real pv_mass      = 1.0;
  //  std::vector<Real> hasenbusch({ 0.01, 0.045, 0.108, 0.25, 0.51 , pv_mass });
  //  std::vector<Real> hasenbusch({ light_mass, 0.01, 0.045, 0.108, 0.25, 0.51 , pv_mass });
  std::vector<Real> hasenbusch({ 0.005, 0.0145, 0.045, 0.108, 0.25, 0.51 , pv_mass }); // Updated
  //  std::vector<Real> hasenbusch({ light_mass, 0.0145, 0.045, 0.108, 0.25, 0.51 , 0.75 , pv_mass });

  auto GridPtr   = TheHMC.Resources.GetCartesian();
  auto GridRBPtr = TheHMC.Resources.GetRBCartesian();

  typedef SchurDiagMooeeOperator<FermionActionF,FermionFieldF> LinearOperatorF;
  typedef SchurDiagMooeeOperator<FermionAction ,FermionField > LinearOperatorD;
  typedef SchurDiagMooeeOperator<FermionEOFAActionF,FermionFieldF> LinearOperatorEOFAF;
  typedef SchurDiagMooeeOperator<FermionEOFAAction ,FermionField > LinearOperatorEOFAD;
  typedef MixedPrecisionConjugateGradientOperatorFunction<MobiusFermionD,MobiusFermionF,LinearOperatorD,LinearOperatorF> MxPCG;
  typedef MixedPrecisionConjugateGradientOperatorFunction<MobiusEOFAFermionD,MobiusEOFAFermionF,LinearOperatorEOFAD,LinearOperatorEOFAF> MxPCG_EOFA;

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

  //////////////////////////
  // Fermion Grids
  //////////////////////////
  auto FGrid     = SpaceTimeGrid::makeFiveDimGrid(Ls,GridPtr);
  auto FrbGrid   = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,GridPtr);

  Coordinate simdF = GridDefaultSimd(Nd,vComplexF::Nsimd());
  auto GridPtrF   = SpaceTimeGrid::makeFourDimGrid(latt4,simdF,mpi);
  auto GridRBPtrF = SpaceTimeGrid::makeFourDimRedBlackGrid(GridPtrF);
  auto FGridF     = SpaceTimeGrid::makeFiveDimGrid(Ls,GridPtrF);
  auto FrbGridF   = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,GridPtrF);

  IwasakiGaugeActionR GaugeAction(beta);

  // temporarily need a gauge field
  LatticeGaugeFieldD  U(GridPtr); U=Zero();
  LatticeGaugeFieldF  UF(GridPtrF); UF=Zero();
  LatticeGaugeFieldD2 UD2(GridPtrF); UD2=Zero();

  std::cout << GridLogMessage << " Running the HMC "<< std::endl;
  TheHMC.ReadCommandLine(argc,argv);  // params on CML or from param file
  TheHMC.initializeGaugeFieldAndRNGs(U);
  std::cout << "loaded NERSC gauge field"<<std::endl;

  // These lines are unecessary if BC are all periodic
  std::vector<Complex> boundary = {1,1,1,-1};
  FermionAction::ImplParams Params(boundary);
  FermionActionF::ImplParams ParamsF(boundary);
  Params.dirichlet=NonDirichlet;
  ParamsF.dirichlet=NonDirichlet;

  //  double StoppingCondition = 1e-14;
  //  double MDStoppingCondition = 1e-9;
  double StoppingCondition = 1e-8;
  double MDStoppingCondition = 1e-7;
  double MDStoppingConditionLoose = 1e-7;
  double MDStoppingConditionStrange = 1e-7;
  double MaxCGIterations = 300000;
  ConjugateGradient<FermionField>  CG(StoppingCondition,MaxCGIterations);
  ConjugateGradient<FermionField>  MDCG(MDStoppingCondition,MaxCGIterations);

  ////////////////////////////////////
  // Collect actions
  ////////////////////////////////////
  //  ActionLevel<HMCWrapper::Field> Level1(1);
  ActionLevel<HMCWrapper::Field> Level2(1);
  ActionLevel<HMCWrapper::Field> Level3(15);

  ////////////////////////////////////
  // Strange action
  ////////////////////////////////////
  FermionAction StrangeOp (U,*FGrid,*FrbGrid,*GridPtr,*GridRBPtr,strange_mass,M5,b,c, Params);
  FermionAction StrangePauliVillarsOp(U,*FGrid,*FrbGrid,*GridPtr,*GridRBPtr,pv_mass,  M5,b,c, Params);

  // Probably dominates the force - back to EOFA.
  OneFlavourRationalParams SFRp;
  SFRp.lo       = 0.25;
  SFRp.hi       = 25.0;
  SFRp.MaxIter  = 10000;
  SFRp.tolerance= 1.0e-5;
  SFRp.mdtolerance= 2.0e-4;
  SFRp.degree   = 8;
  SFRp.precision= 50;
  
  MobiusEOFAFermionD Strange_Op_L (U , *FGrid , *FrbGrid , *GridPtr , *GridRBPtr , strange_mass, strange_mass, pv_mass, 0.0, -1, M5, b, c);
  MobiusEOFAFermionF Strange_Op_LF(UF, *FGridF, *FrbGridF, *GridPtrF, *GridRBPtrF, strange_mass, strange_mass, pv_mass, 0.0, -1, M5, b, c);
  MobiusEOFAFermionD Strange_Op_R (U , *FGrid , *FrbGrid , *GridPtr , *GridRBPtr , pv_mass, strange_mass,      pv_mass, -1.0, 1, M5, b, c);
  MobiusEOFAFermionF Strange_Op_RF(UF, *FGridF, *FrbGridF, *GridPtrF, *GridRBPtrF, pv_mass, strange_mass,      pv_mass, -1.0, 1, M5, b, c);
  ConjugateGradient<FermionField>      ActionCG(StoppingCondition,MaxCGIterations);
  ConjugateGradient<FermionField>  DerivativeCG(MDStoppingCondition,MaxCGIterations);
  LinearOperatorEOFAD Strange_LinOp_L (Strange_Op_L);
  LinearOperatorEOFAD Strange_LinOp_R (Strange_Op_R);
  LinearOperatorEOFAF Strange_LinOp_LF(Strange_Op_LF);
  LinearOperatorEOFAF Strange_LinOp_RF(Strange_Op_RF);

  const int MX_inner = 1000;
  MxPCG_EOFA ActionCGL(StoppingCondition,
		       MX_inner,
		       MaxCGIterations,
		       GridPtrF,
		       FrbGridF,
		       Strange_Op_LF,Strange_Op_L,
		       Strange_LinOp_LF,Strange_LinOp_L);

  MxPCG_EOFA DerivativeCGL(MDStoppingConditionStrange,
			   MX_inner,
			   MaxCGIterations,
			   GridPtrF,
			   FrbGridF,
			   Strange_Op_LF,Strange_Op_L,
			   Strange_LinOp_LF,Strange_LinOp_L);
  
  MxPCG_EOFA ActionCGR(StoppingCondition,
		       MX_inner,
		       MaxCGIterations,
		       GridPtrF,
		       FrbGridF,
		       Strange_Op_RF,Strange_Op_R,
		       Strange_LinOp_RF,Strange_LinOp_R);
  
  MxPCG_EOFA DerivativeCGR(MDStoppingConditionStrange,
			   MX_inner,
			   MaxCGIterations,
			   GridPtrF,
			   FrbGridF,
			   Strange_Op_RF,Strange_Op_R,
			   Strange_LinOp_RF,Strange_LinOp_R);

  ExactOneFlavourRatioPseudoFermionAction<FermionImplPolicy> 
    EOFA(Strange_Op_L, Strange_Op_R, 
	 ActionCG, 
	 ActionCGL, ActionCGR,
	 DerivativeCGL, DerivativeCGR,
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
    light_den.push_back(hasenbusch[h]); dirichlet_den.push_back(0);
  }

  for(int h=0;h<n_hasenbusch;h++){
    light_num.push_back(hasenbusch[h]); dirichlet_num.push_back(0);
  }
  light_num.push_back(pv_mass);  dirichlet_num.push_back(0);

  std::vector<FermionAction *> Numerators;
  std::vector<FermionAction *> Denominators;
  std::vector<FermionActionF *> NumeratorsF;
  std::vector<FermionActionF *> DenominatorsF;
  std::vector<FermionActionD2 *> NumeratorsD2;
  std::vector<FermionActionD2 *> DenominatorsD2;
  std::vector<TwoFlavourEvenOddRatioPseudoFermionAction<FermionImplPolicy> *> Quotients;
  std::vector<MxPCG *> ActionMPCG;
  std::vector<MxPCG *> MPCG;
  
#define MIXED_PRECISION
#ifdef MIXED_PRECISION
  std::vector<OneFlavourEvenOddRatioRationalMixedPrecPseudoFermionAction<FermionImplPolicy,FermionImplPolicyF,FermionImplPolicyD2> *> Bdys;
#else
  std::vector<OneFlavourEvenOddRatioRationalPseudoFermionAction<FermionImplPolicy> *> Bdys;
#endif

  typedef SchurDiagMooeeOperator<FermionActionF,FermionFieldF> LinearOperatorF;
  typedef SchurDiagMooeeOperator<FermionAction ,FermionField > LinearOperatorD;
  std::vector<LinearOperatorD *> LinOpD;
  std::vector<LinearOperatorF *> LinOpF; 
  
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
    FermionActionF::ImplParams ParamsDenF(boundary);
    FermionActionF::ImplParams ParamsNumF(boundary);
    
    ParamsNum.dirichlet = NonDirichlet;
    ParamsDen.dirichlet = NonDirichlet;

    ParamsNum.partialDirichlet = 0;
    ParamsDen.partialDirichlet = 0;
    
    Numerators.push_back  (new FermionAction(U,*FGrid,*FrbGrid,*GridPtr,*GridRBPtr,light_num[h],M5,b,c, ParamsNum));
    Denominators.push_back(new FermionAction(U,*FGrid,*FrbGrid,*GridPtr,*GridRBPtr,light_den[h],M5,b,c, ParamsDen));

    ParamsDenF.dirichlet = ParamsDen.dirichlet;
    ParamsDenF.partialDirichlet = ParamsDen.partialDirichlet;
    DenominatorsF.push_back(new FermionActionF(UF,*FGridF,*FrbGridF,*GridPtrF,*GridRBPtrF,light_den[h],M5,b,c, ParamsDenF));

    ParamsNumF.dirichlet = ParamsNum.dirichlet;
    ParamsNumF.partialDirichlet = ParamsNum.partialDirichlet;
    NumeratorsF.push_back  (new FermionActionF(UF,*FGridF,*FrbGridF,*GridPtrF,*GridRBPtrF,light_num[h],M5,b,c, ParamsNumF));

    LinOpD.push_back(new LinearOperatorD(*Denominators[h]));
    LinOpF.push_back(new LinearOperatorF(*DenominatorsF[h]));

    double conv  = MDStoppingCondition;
    if (h<3) conv= MDStoppingConditionLoose; // Relax on first two hasenbusch factors
    const int MX_inner = 5000;
    MPCG.push_back(new MxPCG(conv,
			     MX_inner,
			     MaxCGIterations,
			     GridPtrF,
			     FrbGridF,
			     *DenominatorsF[h],*Denominators[h],
			     *LinOpF[h], *LinOpD[h]) );

    ActionMPCG.push_back(new MxPCG(StoppingCondition,
				   MX_inner,
				   MaxCGIterations,
				   GridPtrF,
				   FrbGridF,
				   *DenominatorsF[h],*Denominators[h],
				   *LinOpF[h], *LinOpD[h]) );

    
    Quotients.push_back (new TwoFlavourEvenOddRatioPseudoFermionAction<FermionImplPolicy>(*Numerators[h],*Denominators[h],*MPCG[h],*ActionMPCG[h],CG));
  }
  int nquo=Quotients.size();
  for(int h=0;h<nquo;h++){
    Level2.push_back(Quotients[h]);
  }

  /////////////////////////////////////////////////////////////
  // Gauge action
  /////////////////////////////////////////////////////////////
  Level3.push_back(&GaugeAction);
  //  TheHMC.TheAction.push_back(Level1);
  TheHMC.TheAction.push_back(Level2);
  TheHMC.TheAction.push_back(Level3);
  std::cout << GridLogMessage << " Action complete "<< std::endl;
  /////////////////////////////////////////////////////////////

  TheHMC.Run();  // no smearing

  Grid_finalize();
} // main



