/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: 

Copyright (C) 2015-2016

Author: Peter Boyle <pabobyle@ph.ed.ac.uk>

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
#include <Grid/qcd/action/momentum/DirichletFilter.h>
#include <Grid/qcd/action/momentum/DDHMCfilter.h>
#include <Grid/qcd/action/fermion/DirichletFermionOperator.h>
#include <Grid/qcd/action/fermion/SchurFactoredFermionOperator.h>

#include <Grid/qcd/action/pseudofermion/DomainDecomposedBoundaryTwoFlavourPseudoFermion.h>
#include <Grid/qcd/action/pseudofermion/DomainDecomposedBoundaryTwoFlavourRatioPseudoFermion.h>
#include <Grid/qcd/action/pseudofermion/DomainDecomposedBoundaryTwoFlavourBosonPseudoFermion.h>



NAMESPACE_BEGIN(Grid);

#define MIXED_PRECISION
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
    {     };

    void operator()(LinearOperatorBase<FieldD> &LinOpU, const FieldD &src, FieldD &psi)
    {

      SchurOperatorD * SchurOpU = static_cast<SchurOperatorD *>(&LinOpU);
      
      // Assumption made in code to extract gauge field
      // We could avoid storing LinopD reference alltogether ?
      assert(&(SchurOpU->_Mat)==&(LinOpD._Mat));

      ////////////////////////////////////////////////////////////////////////////////////
      // Must snarf a single precision copy of the gauge field in Linop_d argument
      ////////////////////////////////////////////////////////////////////////////////////
      typedef typename FermionOperatorF::GaugeField GaugeFieldF;
      typedef typename FermionOperatorD::GaugeField GaugeFieldD;
      typedef typename FermionOperatorF::GaugeLinkField GaugeLinkFieldF;
      typedef typename FermionOperatorD::GaugeLinkField GaugeLinkFieldD;

      GridBase * GridPtrF = SinglePrecGrid4;
      GridBase * GridPtrD = FermOpD.GaugeGrid();

      ////////////////////////////////////////////////////////////////////////////////////
      // Moving this to a Clone method of fermion operator would allow to duplicate the 
      // physics parameters and decrease gauge field copies
      ////////////////////////////////////////////////////////////////////////////////////
      auto &Umu_d = FermOpD.GetDoubledGaugeField();
      auto &Umu_f = FermOpF.GetDoubledGaugeField();
      auto &Umu_fe= FermOpF.GetDoubledGaugeFieldE();
      auto &Umu_fo= FermOpF.GetDoubledGaugeFieldO();
      precisionChange(Umu_f,Umu_d);
      pickCheckerboard(Even,Umu_fe,Umu_f);
      pickCheckerboard(Odd ,Umu_fo,Umu_f);

      ////////////////////////////////////////////////////////////////////////////////////
      // Could test to make sure that LinOpF and LinOpD agree to single prec?
      ////////////////////////////////////////////////////////////////////////////////////
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
  typedef MobiusFermionR FermionAction;
  typedef MobiusFermionF FermionActionF;
  typedef DirichletFermionOperator<WilsonImplR> DirichletFermion;
  typedef DirichletFermionOperator<WilsonImplF> DirichletFermionF;

  typedef MobiusEOFAFermionR FermionEOFAAction;
  typedef typename FermionAction::FermionField FermionField;
  typedef typename FermionActionF::FermionField FermionFieldF;

  typedef Grid::XmlReader       Serialiser;
  
  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  IntegratorParameters MD;
  //  typedef GenericHMCRunner<LeapFrog> HMCWrapper; 
  //  MD.name    = std::string("Leap Frog");
  //  typedef GenericHMCRunner<ForceGradient> HMCWrapper; 
  //  MD.name    = std::string("Force Gradient");
  typedef GenericHMCRunner<MinimumNorm2> HMCWrapper; 
  MD.name    = std::string("MinimumNorm2");
  MD.MDsteps = 12;
  MD.trajL   = 1.0;
  
  HMCparameters HMCparams;
  HMCparams.StartTrajectory  = 26;
  HMCparams.Trajectories     = 1000;
  HMCparams.NoMetropolisUntil=  10;
  //  "[HotStart, ColdStart, TepidStart, CheckpointStart]\n";
  //  HMCparams.StartingType     =std::string("ColdStart");
  HMCparams.StartingType     =std::string("CheckpointStart");
  HMCparams.MD = MD;
  HMCWrapper TheHMC(HMCparams);

  // Grid from the command line arguments --grid and --mpi
  TheHMC.Resources.AddFourDimGrid("gauge"); // use default simd lanes decomposition
  
  CheckpointerParameters CPparams;
  CPparams.config_prefix = "ckpoint_EOFA4D_lat";
  CPparams.rng_prefix    = "ckpoint_EOFA4D_rng";
  CPparams.saveInterval  = 1;
  CPparams.format        = "IEEE64BIG";
  TheHMC.Resources.LoadNerscCheckpointer(CPparams);

  RNGModuleParameters RNGpar;
  RNGpar.serial_seeds = "1 2 3 4 5";
  RNGpar.parallel_seeds = "6 7 8 9 10";
  TheHMC.Resources.SetRNGSeeds(RNGpar);

  // Momentum Dirichlet
  Coordinate Block({16,16,16,16});

  //  TheHMC.Resources.SetMomentumFilter(new DirichletFilter<WilsonImplR::Field>(Block));
  TheHMC.Resources.SetMomentumFilter(new DDHMCFilter<WilsonImplR::Field>(Block,1));
  // Construct observables
  // here there is too much indirection 
  typedef PlaquetteMod<HMCWrapper::ImplPolicy> PlaqObs;
  TheHMC.Resources.AddObservable<PlaqObs>();
  //////////////////////////////////////////////

  const int Ls      = 16;
  Real beta         = 2.13;
  Real light_mass   = 0.01;
  Real strange_mass = 0.04;
  Real pv_mass      = 1.0;
  RealD M5  = 1.8;
  RealD b   = 1.0; 
  RealD c   = 0.0;

  std::vector<Real> hasenbusch({ 0.04, 0.4, 0.7 });
  
  auto GridPtr   = TheHMC.Resources.GetCartesian();
  auto GridRBPtr = TheHMC.Resources.GetRBCartesian();
  auto FGrid     = SpaceTimeGrid::makeFiveDimGrid(Ls,GridPtr);
  auto FrbGrid   = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,GridPtr);

  Coordinate latt  = GridDefaultLatt();
  Coordinate mpi   = GridDefaultMpi();
  Coordinate simdF = GridDefaultSimd(Nd,vComplexF::Nsimd());
  Coordinate simdD = GridDefaultSimd(Nd,vComplexD::Nsimd());
  auto GridPtrF   = SpaceTimeGrid::makeFourDimGrid(latt,simdF,mpi);
  auto GridRBPtrF = SpaceTimeGrid::makeFourDimRedBlackGrid(GridPtrF);
  auto FGridF     = SpaceTimeGrid::makeFiveDimGrid(Ls,GridPtrF);
  auto FrbGridF   = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,GridPtrF);

  IwasakiGaugeActionR GaugeAction(beta);

  // temporarily need a gauge field
  LatticeGaugeField U(GridPtr);
  LatticeGaugeFieldF UF(GridPtrF);

  // These lines are unecessary if BC are all periodic
  std::vector<Complex> boundary = {1,1,1,-1};
  FermionAction::ImplParams Params(boundary);
  FermionAction::ImplParams DirichletParams(boundary);
  DirichletParams.locally_periodic=true;

  double ActionStoppingCondition     = 1e-10;
  double DerivativeStoppingCondition = 1e-8;
  double MaxCGIterations = 30000;

  ////////////////////////////////////
  // Collect actions
  ////////////////////////////////////
  ActionLevel<HMCWrapper::Field> Level1(1);
  ActionLevel<HMCWrapper::Field> Level2(8);

  ConjugateGradient<FermionField>      ActionCG(ActionStoppingCondition,MaxCGIterations);
  ConjugateGradient<FermionField>  DerivativeCG(DerivativeStoppingCondition,MaxCGIterations);

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

  //////////////////////////////////////////////////////////////
  // Forced to replicate the MxPCG and DenominatorsF etc.. because
  // there is no convenient way to "Clone" physics params from double op
  // into single op for any operator pair.
  // Same issue prevents using MxPCG in the Heatbath step
  //////////////////////////////////////////////////////////////
  typedef SchurDiagMooeeOperator<DirichletFermionF,FermionFieldF> DirichletLinearOperatorF;
  typedef SchurDiagMooeeOperator<DirichletFermion ,FermionField > DirichletLinearOperatorD;
  typedef SchurDiagMooeeDagOperator<DirichletFermionF,FermionFieldF> DirichletLinearOperatorDagF;
  typedef SchurDiagMooeeDagOperator<DirichletFermion ,FermionField > DirichletLinearOperatorDagD;
  typedef SchurDiagMooeeOperator<FermionActionF,FermionFieldF>    PeriodicLinearOperatorF;
  typedef SchurDiagMooeeOperator<FermionAction ,FermionField >    PeriodicLinearOperatorD;
  typedef SchurDiagMooeeDagOperator<FermionActionF,FermionFieldF> PeriodicLinearOperatorDagF;
  typedef SchurDiagMooeeDagOperator<FermionAction ,FermionField > PeriodicLinearOperatorDagD;

  typedef MixedPrecisionConjugateGradientOperatorFunction<DirichletFermion,
							  DirichletFermionF,
							  DirichletLinearOperatorD,
							  DirichletLinearOperatorF> DirichletMxPCG;
  typedef MixedPrecisionConjugateGradientOperatorFunction<DirichletFermion,
							  DirichletFermionF,
							  DirichletLinearOperatorD,
							  DirichletLinearOperatorF> DirichletMxDagPCG;
  typedef MixedPrecisionConjugateGradientOperatorFunction<FermionAction,
							  FermionActionF,
							  PeriodicLinearOperatorD,
							  PeriodicLinearOperatorF> PeriodicMxPCG;
  typedef MixedPrecisionConjugateGradientOperatorFunction<FermionAction,
							  FermionActionF,
							  PeriodicLinearOperatorD,
							  PeriodicLinearOperatorF> PeriodicMxDagPCG;

  //  std::vector<FermionActionF *> DenominatorsF;
  /////////////////////////////////////////////////
  // These are consumed/owned by the Dirichlet wrappers
  /////////////////////////////////////////////////
  std::vector<FermionAction *> DNumerators;
  std::vector<FermionActionF *> DNumeratorsF;
  std::vector<FermionAction *> DDenominators;
  std::vector<FermionActionF *> DDenominatorsF;
  /////////////////////////////////////////////////
  // Dirichlet wrappers
  /////////////////////////////////////////////////
  std::vector<DirichletFermion *> DirichletNumerators;
  std::vector<DirichletFermionF *> DirichletNumeratorsF;
  std::vector<DirichletFermion *> DirichletDenominators;
  std::vector<DirichletFermionF *> DirichletDenominatorsF;
  
  std::vector<TwoFlavourEvenOddRatioPseudoFermionAction<FermionImplPolicy> *> Quotients;
  
  std::vector<DomainDecomposedBoundaryTwoFlavourRatioPseudoFermion<DomainWallFermionR::Impl_t> *> BoundaryQuotients;

  std::vector<DirichletMxPCG *> ActionMPCG;
  std::vector<DirichletMxPCG *> MPCG;
  std::vector<DirichletLinearOperatorD *> LinOpD;
  std::vector<DirichletLinearOperatorF *> LinOpF; 

  const int MX_inner = 1000;

  for(int h=0;h<n_hasenbusch+1;h++){

    std::cout << GridLogMessage << " 2f quotient Action  "<< light_num[h] << " / " << light_den[h]<< std::endl;

    DNumerators.push_back  (new FermionAction(U,*FGrid,*FrbGrid,*GridPtr,*GridRBPtr,light_num[h],M5,b,c, DirichletParams));
    DNumeratorsF.push_back (new FermionActionF(UF,*FGridF,*FrbGridF,*GridPtrF,*GridRBPtrF,light_num[h],M5,b,c, DirichletParams));
    DDenominators.push_back(new FermionAction(U,*FGrid,*FrbGrid,*GridPtr,*GridRBPtr,light_den[h],M5,b,c, DirichletParams));
    DDenominatorsF.push_back(new FermionActionF(UF,*FGridF,*FrbGridF,*GridPtrF,*GridRBPtrF,light_den[h],M5,b,c, DirichletParams));

    DirichletNumerators.push_back  (new  DirichletFermion (*DNumerators[h],Block));
    DirichletNumeratorsF.push_back (new  DirichletFermionF(*DNumeratorsF[h],Block));
    DirichletDenominators.push_back(new  DirichletFermion (*DDenominators[h],Block));
    DirichletDenominatorsF.push_back(new DirichletFermionF(*DDenominatorsF[h],Block));

    // Dirichlet Schur even odd MpsDagMpc operators on local domains
    LinOpD.push_back(new DirichletLinearOperatorD(*DirichletDenominators[h]));
    LinOpF.push_back(new DirichletLinearOperatorF(*DirichletDenominatorsF[h]));

    // Derivative
    MPCG.push_back(new DirichletMxPCG(DerivativeStoppingCondition,
				      MX_inner,
				      MaxCGIterations,
				      GridPtrF,
				      FrbGridF,
				      *DirichletDenominatorsF[h],*DirichletDenominators[h],
				      *LinOpF[h], *LinOpD[h]) );

    // Action
    ActionMPCG.push_back(new DirichletMxPCG(ActionStoppingCondition,
					    MX_inner,
					    MaxCGIterations,
					    GridPtrF,
					    FrbGridF,
					    *DirichletDenominatorsF[h],*DirichletDenominators[h],
					    *LinOpF[h], *LinOpD[h]) );
    
    ////////////////////////////////////////////////////////////////////////////
    // Standard CG for 2f force
    ////////////////////////////////////////////////////////////////////////////
    //    Quotients.push_back (new TwoFlavourEvenOddRatioPseudoFermionAction<FermionImplPolicy>(*Numerators[h],*Denominators[h],*MPCG[h],*ActionMPCG[h],ActionCG));
    Quotients.push_back   (new
			   TwoFlavourEvenOddRatioPseudoFermionAction<FermionImplPolicy>
			   (*DirichletNumerators[h],
			    *DirichletDenominators[h],
			    *ActionMPCG[h],
			    *ActionMPCG[h],
			    ActionCG));

    Level1.push_back(Quotients[h]);
  }

  /////////////////////////////////////////////////////////////
  // Boundary action
  /////////////////////////////////////////////////////////////
  int l_idx = 0;
  int pv_idx = n_hasenbusch;

  std::cout << GridLogMessage<<" Boundary action masses " <<light_num[l_idx]<<" / "<<light_den[pv_idx]<<std::endl;
  /*
  typedef SchurDiagMooeeOperator<DirichletFermionF,FermionFieldF> DirichletLinearOperatorF;
  typedef SchurDiagMooeeOperator<DirichletFermion ,FermionField > DirichletLinearOperatorD;
  typedef SchurDiagMooeeDagOperator<DirichletFermionF,FermionFieldF> DirichletLinearOperatorDagF;
  typedef SchurDiagMooeeDagOperator<DirichletFermion ,FermionField > DirichletLinearOperatorDagD;
  typedef SchurDiagMooeeOperator<FermionActionF,FermionFieldF>    PeriodicLinearOperatorF;
  typedef SchurDiagMooeeOperator<FermionAction ,FermionField >    PeriodicLinearOperatorD;
  typedef SchurDiagMooeeDagOperator<FermionActionF,FermionFieldF> PeriodicLinearOperatorDagF;
  typedef SchurDiagMooeeDagOperator<FermionAction ,FermionField > PeriodicOperatorDagD;
  */
  
  typedef MixedPrecisionConjugateGradientOperatorFunction<DirichletFermion,
							  DirichletFermionF,
							  DirichletLinearOperatorD,
							  DirichletLinearOperatorF> DirichletMxPCG;

  typedef MixedPrecisionConjugateGradientOperatorFunction<DirichletFermion,
							  DirichletFermionF,
							  DirichletLinearOperatorD,
							  DirichletLinearOperatorF> DirichletMxDagPCG;

  typedef MixedPrecisionConjugateGradientOperatorFunction<FermionAction,
							  FermionActionF,
							  PeriodicLinearOperatorD,
							  PeriodicLinearOperatorF> PeriodicMxPCG;

  typedef MixedPrecisionConjugateGradientOperatorFunction<FermionAction,
							  FermionActionF,
							  PeriodicLinearOperatorD,
							  PeriodicLinearOperatorF> PeriodicMxDagPCG;

  FermionAction  PeriNumerator   (U,*FGrid,*FrbGrid,*GridPtr,*GridRBPtr,light_num[pv_idx],M5,b,c, Params);
  FermionAction  PeriDenominator (U,*FGrid,*FrbGrid,*GridPtr,*GridRBPtr,light_den[l_idx] ,M5,b,c, Params);
  FermionActionF PeriNumeratorF  (UF,*FGridF,*FrbGridF,*GridPtrF,*GridRBPtrF,light_num[pv_idx],M5,b,c, Params);
  FermionActionF PeriDenominatorF(UF,*FGridF,*FrbGridF,*GridPtrF,*GridRBPtrF,light_den[l_idx] ,M5,b,c, Params);

  DirichletFermion  & DirichletNumerator   = *DirichletNumerators[pv_idx];
  DirichletFermionF & DirichletNumeratorF  = *DirichletNumeratorsF[pv_idx];
  DirichletFermion  & DirichletDenominator = *DirichletDenominators[l_idx];
  DirichletFermionF & DirichletDenominatorF= *DirichletDenominatorsF[l_idx];

  DirichletLinearOperatorD    DirichletNumeratorLinOpD(DirichletNumerator);
  DirichletLinearOperatorF    DirichletNumeratorLinOpF(DirichletNumeratorF);
  DirichletLinearOperatorDagD DirichletNumeratorLinOpDagD(DirichletNumerator);
  DirichletLinearOperatorDagF DirichletNumeratorLinOpDagF(DirichletNumeratorF);
  
  DirichletLinearOperatorD    DirichletDenominatorLinOpD(DirichletDenominator);
  DirichletLinearOperatorF    DirichletDenominatorLinOpF(DirichletDenominatorF);
  DirichletLinearOperatorDagD DirichletDenominatorLinOpDagD(DirichletDenominator);
  DirichletLinearOperatorDagF DirichletDenominatorLinOpDagF(DirichletDenominatorF);

  PeriodicLinearOperatorF     PeriodicNumeratorLinOpF(PeriNumeratorF);
  PeriodicLinearOperatorD     PeriodicNumeratorLinOpD(PeriNumerator);
  PeriodicLinearOperatorDagF  PeriodicNumeratorLinOpDagF(PeriNumeratorF);
  PeriodicLinearOperatorDagD  PeriodicNumeratorLinOpDagD(PeriNumerator);

  PeriodicLinearOperatorF     PeriodicDenominatorLinOpF(PeriDenominatorF);
  PeriodicLinearOperatorD     PeriodicDenominatorLinOpD(PeriDenominator);
  PeriodicLinearOperatorDagF  PeriodicDenominatorLinOpDagF(PeriDenominatorF);
  PeriodicLinearOperatorDagD  PeriodicDenominatorLinOpDagD(PeriDenominator);

  ConjugateGradient<FermionField>  OmegaSolverCG   (ActionStoppingCondition,MaxCGIterations);
  ConjugateGradient<FermionField>  OmegaDagSolverCG(ActionStoppingCondition,MaxCGIterations);
  ConjugateGradient<FermionField>  DSolverCG       (ActionStoppingCondition,MaxCGIterations);
  ConjugateGradient<FermionField>  DdagSolverCG    (ActionStoppingCondition,MaxCGIterations);

  DirichletMxPCG NumeratorOmegaSolver(ActionStoppingCondition,
				      MX_inner,
				      MaxCGIterations,
				      GridPtrF,
				      FrbGridF,
				      DirichletNumeratorF,DirichletNumerator,
				      DirichletNumeratorLinOpF, DirichletNumeratorLinOpD);

  DirichletMxDagPCG NumeratorOmegaDagSolver(ActionStoppingCondition,
					    MX_inner,
					    MaxCGIterations,
					    GridPtrF,
					    FrbGridF,
					    DirichletNumeratorF,DirichletNumerator,
					    DirichletNumeratorLinOpDagF,DirichletNumeratorLinOpDagD);

  PeriodicMxPCG NumeratorDSolver(ActionStoppingCondition,
				 MX_inner,
				 MaxCGIterations,
				 GridPtrF,
				 FrbGridF,
				 PeriNumeratorF,PeriNumerator,
				 PeriodicNumeratorLinOpF, PeriodicNumeratorLinOpD);

  PeriodicMxDagPCG NumeratorDdagSolver(ActionStoppingCondition,
				       MX_inner,
				       MaxCGIterations,
				       GridPtrF,
				       FrbGridF,
				       PeriNumeratorF,PeriNumerator,
				       PeriodicNumeratorLinOpDagF, PeriodicNumeratorLinOpDagD);
  
  SchurFactoredFermionOperator<DomainWallFermionR::Impl_t> BoundaryNumerator(PeriNumerator,
									     DirichletNumerator,
									     NumeratorOmegaSolver,NumeratorOmegaDagSolver,
									     NumeratorDSolver,NumeratorDdagSolver,
									     // ActionCG,ActionCG,
									     // ActionCG,ActionCG,
									     Block);
  
  DirichletMxPCG DenominatorOmegaSolver(ActionStoppingCondition,
				      MX_inner,
				      MaxCGIterations,
				      GridPtrF,
				      FrbGridF,
				      DirichletDenominatorF,DirichletDenominator,
				      DirichletDenominatorLinOpF, DirichletDenominatorLinOpD);

  DirichletMxDagPCG DenominatorOmegaDagSolver(ActionStoppingCondition,
					    MX_inner,
					    MaxCGIterations,
					    GridPtrF,
					    FrbGridF,
					    DirichletDenominatorF,DirichletDenominator,
					    DirichletDenominatorLinOpDagF,DirichletDenominatorLinOpDagD);

  PeriodicMxPCG DenominatorDSolver(ActionStoppingCondition,
				   MX_inner,
				   MaxCGIterations,
				   GridPtrF,
				   FrbGridF,
				   PeriDenominatorF,PeriDenominator,
				   PeriodicDenominatorLinOpF, PeriodicDenominatorLinOpD);

  PeriodicMxDagPCG DenominatorDdagSolver(ActionStoppingCondition,
				      MX_inner,
				      MaxCGIterations,
				      GridPtrF,
				      FrbGridF,
				      PeriDenominatorF,PeriDenominator,
				      PeriodicDenominatorLinOpDagF, PeriodicDenominatorLinOpDagD);

  SchurFactoredFermionOperator<DomainWallFermionR::Impl_t> BoundaryDenominator(PeriDenominator,
									       DirichletDenominator,
									       //ActionCG,ActionCG,
									       //ActionCG,ActionCG,
									       DenominatorOmegaSolver,DenominatorOmegaDagSolver,
									       DenominatorDSolver,DenominatorDdagSolver,
									       Block);
  Level1.push_back(new
		   DomainDecomposedBoundaryTwoFlavourRatioPseudoFermion<DomainWallFermionR::Impl_t>
		   (BoundaryNumerator,
		    BoundaryDenominator));

  /////////////////////////////////////////////////////////////
  // Gauge action
  /////////////////////////////////////////////////////////////
  Level2.push_back(&GaugeAction);
  TheHMC.TheAction.push_back(Level1);
  TheHMC.TheAction.push_back(Level2);
  std::cout << GridLogMessage << " Action complete "<< std::endl;

  /////////////////////////////////////////////////////////////
  // HMC parameters are serialisable
  std::cout << GridLogMessage << " Running the HMC "<< std::endl;
  TheHMC.Run();  // no smearing

  Grid_finalize();
} // main



