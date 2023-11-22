/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: 

Copyright (C) 2015-2016

Author: Peter Boyle <pabobyle@ph.ed.ac.uk>
Author: Guido Cossu
Author: David Murphy
Author: Chulwoo Jung <chulwoo@bnl.gov>

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

#ifdef GRID_DEFAULT_PRECISION_DOUBLE
#define MIXED_PRECISION
#endif
// second level EOFA
#undef EOFA_H
#define USE_OBC
#undef DO_IMPLICIT

NAMESPACE_BEGIN(Grid);

  /*
   * Need a plan for gauge field update for mixed precision in HMC                      (2x speed up)
   *    -- Store the single prec action operator.
   *    -- Clone the gauge field from the operator function argument.
   *    -- Build the mixed precision operator dynamically from the passed operator and single prec clone.
   */

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
      MixedPrecisionConjugateGradient<FieldD,FieldF> MPCG(Tolerance,MaxInnerIterations,MaxOuterIterations,SinglePrecGrid5,LinOpF,LinOpD);
      std::cout << GridLogMessage << "Calling mixed precision Conjugate Gradient" <<std::endl;
      MPCG(src,psi);
    }
  };

NAMESPACE_END(Grid);


int main(int argc, char **argv) {
  using namespace Grid;

  Grid_init(&argc, &argv);
  int threads = GridThread::GetThreads();
  // here make a routine to print all the relevant information on the run
  std::cout << GridLogMessage << "Grid is setup to use " << threads << " threads" << std::endl;

   // Typedefs to simplify notation
  typedef WilsonImplR FermionImplPolicy;
  typedef MobiusFermionD FermionAction;
  typedef MobiusFermionF FermionActionF;
  typedef MobiusEOFAFermionD FermionEOFAAction;
  typedef MobiusEOFAFermionF FermionEOFAActionF;
  typedef typename FermionAction::FermionField FermionField;
  typedef typename FermionActionF::FermionField FermionFieldF;

  typedef Grid::XmlReader       Serialiser;
  
  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  HMCparameters HMCparams;
#if 1
  {
    XmlReader  HMCrd("HMCparameters.xml");
    read(HMCrd,"HMCparameters",HMCparams);
  }
#else
  {
//    HMCparameters HMCparams;
  //  "[HotStart, ColdStart, TepidStart, CheckpointStart]\n";
  //  HMCparams.StartingType     =std::string("ColdStart");
    HMCparams.StartingType     =std::string("CheckpointStart");
    HMCparams.StartTrajectory  =7;
    HMCparams.SW  =4;
    HMCparams.Trajectories     =1000;
    HMCparams.NoMetropolisUntil=0;
    HMCparams.MD.name          =std::string("Force Gradient");
    HMCparams.MD.MDsteps       = 10;
    HMCparams.MD.trajL         = 1.0;
  }
#endif

#ifdef DO_IMPLICIT
//    typedef GenericHMCRunner<ImplicitLeapFrog> HMCWrapper; 
  typedef GenericHMCRunner<ImplicitMinimumNorm2> HMCWrapper; 
  HMCparams.MD.name          =std::string("ImplicitMinimumNorm2");
#else
//  typedef GenericHMCRunner<LeapFrog> HMCWrapper; 
//  typedef GenericHMCRunner<ForceGradient> HMCWrapper; 
  typedef GenericHMCRunner<MinimumNorm2> HMCWrapper; 
  HMCparams.MD.name          =std::string("MinimumNorm2");
#endif

  std::cout << GridLogMessage<< HMCparams <<std::endl;
  HMCWrapper TheHMC(HMCparams);
  TheHMC.ReadCommandLine(argc, argv);
  { 
    XmlWriter HMCwr("HMCparameters.xml.out");
    write(HMCwr,"HMCparameters",TheHMC.Parameters);
  }

  // Grid from the command line arguments --grid and --mpi
  TheHMC.Resources.AddFourDimGrid("gauge"); // use default simd lanes decomposition
  
  CheckpointerParameters CPparams;
  CPparams.config_prefix = "ckpoint_lat";
  CPparams.rng_prefix    = "ckpoint_rng";
  CPparams.saveInterval  = 1;
  CPparams.format        = "IEEE64BIG";
  TheHMC.Resources.LoadNerscCheckpointer(CPparams);

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
  Real beta         = 5.983;
  std::cout << GridLogMessage << " beta  "<< beta << std::endl;
  Real light_mass   = 0.00049;
  Real strange_mass = 0.0158;
  Real charm_mass = 0.191;
  Real pv_mass    = 1.0;
  RealD M5  = 1.4;
  RealD b   = 2.0; 
  RealD c   = 1.0;

  // Copied from paper
//  std::vector<Real> hasenbusch({ 0.045 }); // Paper values from F1 incorrect run
  std::vector<Real> hasenbusch({ 0.0038, 0.0145, 0.045, 0.108 , 0.25, 0.51 }); // Paper values from F1 incorrect run
  std::vector<Real> hasenbusch2({ 0.4 }); // Paper values from F1 incorrect run

//  RealD eofa_mass=0.05 ;

  ///////////////////////////////////////////////////////////////////////////////////////////////
  //Bad choices with large dH. Equalising force L2 norm was not wise.
  ///////////////////////////////////////////////////////////////////////////////////////////////
  //std::vector<Real> hasenbusch({ 0.03, 0.2, 0.3, 0.5, 0.8 }); 

  auto GridPtr   = TheHMC.Resources.GetCartesian();
  auto GridRBPtr = TheHMC.Resources.GetRBCartesian();
  auto FGrid     = SpaceTimeGrid::makeFiveDimGrid(Ls,GridPtr);
  auto FrbGrid   = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,GridPtr);

  Coordinate latt  = GridDefaultLatt();
  Coordinate mpi   = GridDefaultMpi();
  Coordinate simdF = GridDefaultSimd(Nd,vComplexF::Nsimd());
  Coordinate simdD = GridDefaultSimd(Nd,vComplexD::Nsimd());
  auto GridPtrF   = SpaceTimeGrid::makeFourDimGrid(latt,simdF,mpi);
//auto UGrid_f    = SpaceTimeGrid::makeFourDimGrid(latt,simdF,mpi);
  auto GridRBPtrF = SpaceTimeGrid::makeFourDimRedBlackGrid(GridPtrF);
  auto FGridF     = SpaceTimeGrid::makeFiveDimGrid(Ls,GridPtrF);
  auto FrbGridF   = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,GridPtrF);


#ifndef USE_OBC
//  IwasakiGaugeActionR GaugeAction(beta);
  WilsonGaugeActionR GaugeAction(beta);
#else
  std::vector<Complex> boundaryG = {1,1,1,0};
  WilsonGaugeActionR::ImplParams ParamsG(boundaryG);
  WilsonGaugeActionR GaugeAction(beta,ParamsG);
#endif

  // temporarily need a gauge field
  LatticeGaugeField U(GridPtr);
  LatticeGaugeFieldF UF(GridPtrF);

  // These lines are unecessary if BC are all periodic
#ifndef USE_OBC
  std::vector<Complex> boundary = {1,1,1,-1};
#else
  std::vector<Complex> boundary = {1,1,1,0};
#endif
  FermionAction::ImplParams Params(boundary);
  FermionActionF::ImplParams ParamsF(boundary);
  
  double ActionStoppingCondition     = 1e-8;
  double DerivativeStoppingCondition = 1e-8;
  double MaxCGIterations =  100000;

  ////////////////////////////////////
  // Collect actions
  ////////////////////////////////////
  ActionLevel<HMCWrapper::Field> Level1(1);
  ActionLevel<HMCWrapper::Field> Level2(HMCparams.SW);

  ////////////////////////////////////
  // Strange action
  ////////////////////////////////////
  typedef SchurDiagMooeeOperator<FermionActionF,FermionFieldF> LinearOperatorF;
  typedef SchurDiagMooeeOperator<FermionAction ,FermionField > LinearOperatorD;
  typedef SchurDiagMooeeOperator<FermionEOFAActionF,FermionFieldF> LinearOperatorEOFAF;
  typedef SchurDiagMooeeOperator<FermionEOFAAction ,FermionField > LinearOperatorEOFAD;

  typedef MixedPrecisionConjugateGradientOperatorFunction<MobiusFermionD,MobiusFermionF,LinearOperatorD,LinearOperatorF> MxPCG;
  typedef MixedPrecisionConjugateGradientOperatorFunction<MobiusEOFAFermionD,MobiusEOFAFermionF,LinearOperatorEOFAD,LinearOperatorEOFAF> MxPCG_EOFA;

  // DJM: setup for EOFA ratio (Mobius)
  OneFlavourRationalParams OFRp;
  OFRp.lo       = 0.99; // How do I know this on F1?
  OFRp.hi       = 20;
  OFRp.MaxIter  = 100000;
  OFRp.tolerance= 1.0e-12;
  OFRp.degree   = 12;
  OFRp.precision= 50;

  
  MobiusEOFAFermionD Strange_Op_L (U , *FGrid , *FrbGrid , *GridPtr , *GridRBPtr , strange_mass, strange_mass, charm_mass, 0.0, -1, M5, b, c);
  MobiusEOFAFermionF Strange_Op_LF(UF, *FGridF, *FrbGridF, *GridPtrF, *GridRBPtrF, strange_mass, strange_mass, charm_mass, 0.0, -1, M5, b, c);
  MobiusEOFAFermionD Strange_Op_R (U , *FGrid , *FrbGrid , *GridPtr , *GridRBPtr , charm_mass, strange_mass,      charm_mass, -1.0, 1, M5, b, c);
  MobiusEOFAFermionF Strange_Op_RF(UF, *FGridF, *FrbGridF, *GridPtrF, *GridRBPtrF, charm_mass, strange_mass,      charm_mass, -1.0, 1, M5, b, c);
  
#ifdef EOFA_H
  MobiusEOFAFermionD Strange2_Op_L (U , *FGrid , *FrbGrid , *GridPtr , *GridRBPtr , eofa_mass, eofa_mass, charm_mass , 0.0, -1, M5, b, c);
  MobiusEOFAFermionF Strange2_Op_LF(UF, *FGridF, *FrbGridF, *GridPtrF, *GridRBPtrF, eofa_mass, eofa_mass, charm_mass , 0.0, -1, M5, b, c);
  MobiusEOFAFermionD Strange2_Op_R (U , *FGrid , *FrbGrid , *GridPtr , *GridRBPtr , charm_mass , eofa_mass,      charm_mass , -1.0, 1, M5, b, c);
  MobiusEOFAFermionF Strange2_Op_RF(UF, *FGridF, *FrbGridF, *GridPtrF, *GridRBPtrF, charm_mass , eofa_mass,      charm_mass , -1.0, 1, M5, b, c);
#endif

  ConjugateGradient<FermionField>      ActionCG(ActionStoppingCondition,MaxCGIterations);
  ConjugateGradient<FermionField>  DerivativeCG(DerivativeStoppingCondition,MaxCGIterations);
#ifdef MIXED_PRECISION
  const int MX_inner = 5000;

  // Mixed precision EOFA
  LinearOperatorEOFAD Strange_LinOp_L (Strange_Op_L);
  LinearOperatorEOFAD Strange_LinOp_R (Strange_Op_R);
  LinearOperatorEOFAF Strange_LinOp_LF(Strange_Op_LF);
  LinearOperatorEOFAF Strange_LinOp_RF(Strange_Op_RF);

#ifdef EOFA_H
  // Mixed precision EOFA
  LinearOperatorEOFAD Strange2_LinOp_L (Strange2_Op_L);
  LinearOperatorEOFAD Strange2_LinOp_R (Strange2_Op_R);
  LinearOperatorEOFAF Strange2_LinOp_LF(Strange2_Op_LF);
  LinearOperatorEOFAF Strange2_LinOp_RF(Strange2_Op_RF);
#endif

  MxPCG_EOFA ActionCGL(ActionStoppingCondition,
		       MX_inner,
		       MaxCGIterations,
		       GridPtrF,
		       FrbGridF,
		       Strange_Op_LF,Strange_Op_L,
		       Strange_LinOp_LF,Strange_LinOp_L);

#ifdef EOFA_H
  MxPCG_EOFA ActionCGL2(ActionStoppingCondition,
		       MX_inner,
		       MaxCGIterations,
		       GridPtrF,
		       FrbGridF,
		       Strange2_Op_LF,Strange2_Op_L,
		       Strange2_LinOp_LF,Strange2_LinOp_L);
#endif

  MxPCG_EOFA DerivativeCGL(DerivativeStoppingCondition,
			   MX_inner,
			   MaxCGIterations,
			   GridPtrF,
			   FrbGridF,
			   Strange_Op_LF,Strange_Op_L,
			   Strange_LinOp_LF,Strange_LinOp_L);

#ifdef EOFA_H
  MxPCG_EOFA DerivativeCGL2(DerivativeStoppingCondition,
			   MX_inner,
			   MaxCGIterations,
			   GridPtrF,
			   FrbGridF,
			   Strange2_Op_LF,Strange2_Op_L,
			   Strange2_LinOp_LF,Strange2_LinOp_L);
#endif
  
  MxPCG_EOFA ActionCGR(ActionStoppingCondition,
		       MX_inner,
		       MaxCGIterations,
		       GridPtrF,
		       FrbGridF,
		       Strange_Op_RF,Strange_Op_R,
		       Strange_LinOp_RF,Strange_LinOp_R);
  
#ifdef EOFA_H
  MxPCG_EOFA ActionCGR2(ActionStoppingCondition,
		       MX_inner,
		       MaxCGIterations,
		       GridPtrF,
		       FrbGridF,
		       Strange2_Op_RF,Strange2_Op_R,
		       Strange2_LinOp_RF,Strange2_LinOp_R);
#endif
  
  MxPCG_EOFA DerivativeCGR(DerivativeStoppingCondition,
			   MX_inner,
			   MaxCGIterations,
			   GridPtrF,
			   FrbGridF,
			   Strange_Op_RF,Strange_Op_R,
			   Strange_LinOp_RF,Strange_LinOp_R);
  
#ifdef EOFA_H
  MxPCG_EOFA DerivativeCGR2(DerivativeStoppingCondition,
			   MX_inner,
			   MaxCGIterations,
			   GridPtrF,
			   FrbGridF,
			   Strange2_Op_RF,Strange2_Op_R,
			   Strange2_LinOp_RF,Strange2_LinOp_R);
#endif
  
  ExactOneFlavourRatioPseudoFermionAction<FermionImplPolicy> 
    EOFA(Strange_Op_L, Strange_Op_R, 
	 ActionCG, 
	 ActionCGL, ActionCGR,
	 DerivativeCGL, DerivativeCGR,
	 OFRp, true);
  
#ifdef EOFA_H
  ExactOneFlavourRatioPseudoFermionAction<FermionImplPolicy> 
    EOFA2(Strange2_Op_L, Strange2_Op_R, 
	 ActionCG, 
	 ActionCGL2, ActionCGR2,
	 DerivativeCGL2, DerivativeCGR2,
	 OFRp, true);
#endif

  Level1.push_back(&EOFA);
#ifdef EOFA_H
  Level1.push_back(&EOFA2);
#endif

#else
  ExactOneFlavourRatioPseudoFermionAction<FermionImplPolicy> 
    EOFA(Strange_Op_L, Strange_Op_R, 
	 ActionCG, 
	 ActionCG, ActionCG,
	 ActionCG, ActionCG,
	 //         DerivativeCG, DerivativeCG,
	 OFRp, true);
  Level1.push_back(&EOFA);
#endif

  ////////////////////////////////////
  // up down action
  ////////////////////////////////////
  std::vector<Real> light_den;
  std::vector<Real> light_num;

  int n_hasenbusch = hasenbusch.size();
  light_den.push_back(light_mass);
  for(int h=0;h<n_hasenbusch;h++){
    light_den.push_back(hasenbusch[h]);
    light_num.push_back(hasenbusch[h]);
  }
  light_num.push_back(pv_mass);

  int n_hasenbusch2 = hasenbusch2.size();
  light_den.push_back(charm_mass);
  for(int h=0;h<n_hasenbusch2;h++){
    light_den.push_back(hasenbusch2[h]);
    light_num.push_back(hasenbusch2[h]);
  }
  light_num.push_back(pv_mass);


  //////////////////////////////////////////////////////////////
  // Forced to replicate the MxPCG and DenominatorsF etc.. because
  // there is no convenient way to "Clone" physics params from double op
  // into single op for any operator pair.
  // Same issue prevents using MxPCG in the Heatbath step
  //////////////////////////////////////////////////////////////
  std::vector<FermionAction *> Numerators;
  std::vector<FermionAction *> Denominators;
  std::vector<TwoFlavourEvenOddRatioPseudoFermionAction<FermionImplPolicy> *> Quotients;
  std::vector<MxPCG *> ActionMPCG;
  std::vector<MxPCG *> MPCG;
  std::vector<FermionActionF *> DenominatorsF;
  std::vector<LinearOperatorD *> LinOpD;
  std::vector<LinearOperatorF *> LinOpF; 

  for(int h=0;h<light_den.size();h++){

    std::cout << GridLogMessage << " 2f quotient Action  "<< light_num[h] << " / " << light_den[h]<< std::endl;

    Numerators.push_back  (new FermionAction(U,*FGrid,*FrbGrid,*GridPtr,*GridRBPtr,light_num[h],M5,b,c, Params));
    Denominators.push_back(new FermionAction(U,*FGrid,*FrbGrid,*GridPtr,*GridRBPtr,light_den[h],M5,b,c, Params));

#ifdef MIXED_PRECISION
    ////////////////////////////////////////////////////////////////////////////
    // Mixed precision CG for 2f force
    ////////////////////////////////////////////////////////////////////////////
    double DerivativeStoppingConditionLoose = 1e-8;

    DenominatorsF.push_back(new FermionActionF(UF,*FGridF,*FrbGridF,*GridPtrF,*GridRBPtrF,light_den[h],M5,b,c, ParamsF));
    LinOpD.push_back(new LinearOperatorD(*Denominators[h]));
    LinOpF.push_back(new LinearOperatorF(*DenominatorsF[h]));

    double conv  = DerivativeStoppingCondition;
    if (h<3) conv= DerivativeStoppingConditionLoose; // Relax on first two hasenbusch factors
    MPCG.push_back(new MxPCG(conv,
			     MX_inner,
			     MaxCGIterations,
			     GridPtrF,
			     FrbGridF,
			     *DenominatorsF[h],*Denominators[h],
			     *LinOpF[h], *LinOpD[h]) );

    ActionMPCG.push_back(new MxPCG(ActionStoppingCondition,
				   MX_inner,
				   MaxCGIterations,
				   GridPtrF,
				   FrbGridF,
				   *DenominatorsF[h],*Denominators[h],
				   *LinOpF[h], *LinOpD[h]) );

    // Heatbath not mixed yet. As inverts numerators not so important as raised mass.
    Quotients.push_back (new TwoFlavourEvenOddRatioPseudoFermionAction<FermionImplPolicy>(*Numerators[h],*Denominators[h],*MPCG[h],*ActionMPCG[h],ActionCG));
#else
    ////////////////////////////////////////////////////////////////////////////
    // Standard CG for 2f force
    ////////////////////////////////////////////////////////////////////////////
    Quotients.push_back   (new TwoFlavourEvenOddRatioPseudoFermionAction<FermionImplPolicy>(*Numerators[h],*Denominators[h],DerivativeCG,ActionCG));
#endif

  }

  for(int h=0;h<n_hasenbusch+1;h++){
    Level1.push_back(Quotients[h]);
  }

  /////////////////////////////////////////////////////////////
  // Gauge action
  /////////////////////////////////////////////////////////////
  Level2.push_back(&GaugeAction);
  TheHMC.TheAction.push_back(Level1);
  TheHMC.TheAction.push_back(Level2);
  std::cout << GridLogMessage << " Action complete "<< std::endl;

  /////////////////////////////////////////////////////////////
  // HMC parameters are serialisable

  NoSmearing<HMCWrapper::ImplPolicy> S;
#ifndef DO_IMPLICIT
  TrivialMetric<HMCWrapper::ImplPolicy::Field> Mtr;
#else
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
    auto UGrid = TheHMC.Resources.GetCartesian("gauge");
//  auto UGrid_f   = SpaceTimeGrid::makeFourDimGrid(latt,simdF,mpi);
//  auto GridPtrF   = SpaceTimeGrid::makeFourDimGrid(latt,simdF,mpi);
//    std::cout << GridLogMessage << " UGrid= " << UGrid <<std::endl;
//    std::cout << GridLogMessage << " UGrid_f= " << UGrid_f <<std::endl;

    LaplacianAdjointRat<HMCWrapper::ImplPolicy, PeriodicGimplF> Mtr(UGrid, GridPtrF ,CG, gpar, mpar);
#endif

  std::cout << GridLogMessage << " Running the HMC "<< std::endl;
  TheHMC.Run(S,Mtr);  // no smearing

  Grid_finalize();
} // main



