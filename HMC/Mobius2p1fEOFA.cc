/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: 

Copyright (C) 2015-2016

Author: Peter Boyle <pabobyle@ph.ed.ac.uk>
Author: Guido Cossu
Author: David Murphy

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

#define MIXED_PRECISION

namespace Grid{ 
  namespace QCD{

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
      GridBase * GridPtrD = FermOpD.Umu._grid;
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
      // Could test to make sure that LinOpF and LinOpD agree to single prec?
      ////////////////////////////////////////////////////////////////////////////////////
      /*
      GridBase *Fgrid = psi._grid;
      FieldD tmp2(Fgrid);
      FieldD tmp1(Fgrid);
      LinOpU.Op(src,tmp1);
      LinOpD.Op(src,tmp2);
      std::cout << " Double gauge field "<< norm2(FermOpD.Umu)<<std::endl;
      std::cout << " Single gauge field "<< norm2(FermOpF.Umu)<<std::endl;
      std::cout << " Test of operators "<<norm2(tmp1)<<std::endl;
      std::cout << " Test of operators "<<norm2(tmp2)<<std::endl;
      tmp1=tmp1-tmp2;
      std::cout << " Test of operators diff "<<norm2(tmp1)<<std::endl;
      */

      ////////////////////////////////////////////////////////////////////////////////////
      // Make a mixed precision conjugate gradient
      ////////////////////////////////////////////////////////////////////////////////////
      MixedPrecisionConjugateGradient<FieldD,FieldF> MPCG(Tolerance,MaxInnerIterations,MaxOuterIterations,SinglePrecGrid5,LinOpF,LinOpD);
      std::cout << GridLogMessage << "Calling mixed precision Conjugate Gradient" <<std::endl;
      MPCG(src,psi);
    }
  };
}};

int main(int argc, char **argv) {
  using namespace Grid;
  using namespace Grid::QCD;

  Grid_init(&argc, &argv);
  int threads = GridThread::GetThreads();
  // here make a routine to print all the relevant information on the run
  std::cout << GridLogMessage << "Grid is setup to use " << threads << " threads" << std::endl;

   // Typedefs to simplify notation
  typedef WilsonImplR FermionImplPolicy;
  typedef MobiusFermionR FermionAction;
  typedef MobiusFermionF FermionActionF;
  typedef MobiusEOFAFermionR FermionEOFAAction;
  typedef MobiusEOFAFermionF FermionEOFAActionF;
  typedef typename FermionAction::FermionField FermionField;
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
  MD.MDsteps = 6;
  MD.trajL   = 1.0;
  
  HMCparameters HMCparams;
  HMCparams.StartTrajectory  = 590;
  HMCparams.Trajectories     = 1000;
  HMCparams.NoMetropolisUntil=  0;
  //  "[HotStart, ColdStart, TepidStart, CheckpointStart]\n";
  //  HMCparams.StartingType     =std::string("ColdStart");
  HMCparams.StartingType     =std::string("CheckpointStart");
  HMCparams.MD = MD;
  HMCWrapper TheHMC(HMCparams);

  // Grid from the command line arguments --grid and --mpi
  TheHMC.Resources.AddFourDimGrid("gauge"); // use default simd lanes decomposition
  
  CheckpointerParameters CPparams;
  CPparams.config_prefix = "ckpoint_EODWF_lat";
  CPparams.rng_prefix    = "ckpoint_EODWF_rng";
  CPparams.saveInterval  = 10;
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

  const int Ls      = 16;
  Real beta         = 2.13;
  Real light_mass   = 0.01;
  Real strange_mass = 0.04;
  Real pv_mass      = 1.0;
  RealD M5  = 1.8;
  RealD b   = 1.0; 
  RealD c   = 0.0;

  std::vector<Real> hasenbusch({ 0.1, 0.3, 0.6 });

  auto GridPtr   = TheHMC.Resources.GetCartesian();
  auto GridRBPtr = TheHMC.Resources.GetRBCartesian();
  auto FGrid     = SpaceTimeGrid::makeFiveDimGrid(Ls,GridPtr);
  auto FrbGrid   = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,GridPtr);

  std::vector<int> latt  = GridDefaultLatt();
  std::vector<int> mpi   = GridDefaultMpi();
  std::vector<int> simdF = GridDefaultSimd(Nd,vComplexF::Nsimd());
  std::vector<int> simdD = GridDefaultSimd(Nd,vComplexD::Nsimd());
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
  FermionActionF::ImplParams ParamsF(boundary);
  
  double ActionStoppingCondition     = 1e-10;
  double DerivativeStoppingCondition = 1e-6;
  double MaxCGIterations = 30000;

  ////////////////////////////////////
  // Collect actions
  ////////////////////////////////////
  ActionLevel<HMCWrapper::Field> Level1(1);
  ActionLevel<HMCWrapper::Field> Level2(8);

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
  OFRp.lo       = 0.1;
  OFRp.hi       = 25.0;
  OFRp.MaxIter  = 10000;
  OFRp.tolerance= 1.0e-9;
  OFRp.degree   = 14;
  OFRp.precision= 50;

  
  MobiusEOFAFermionR Strange_Op_L (U , *FGrid , *FrbGrid , *GridPtr , *GridRBPtr , strange_mass, strange_mass, pv_mass, 0.0, -1, M5, b, c);
  MobiusEOFAFermionF Strange_Op_LF(UF, *FGridF, *FrbGridF, *GridPtrF, *GridRBPtrF, strange_mass, strange_mass, pv_mass, 0.0, -1, M5, b, c);
  MobiusEOFAFermionR Strange_Op_R (U , *FGrid , *FrbGrid , *GridPtr , *GridRBPtr , pv_mass, strange_mass,      pv_mass, -1.0, 1, M5, b, c);
  MobiusEOFAFermionF Strange_Op_RF(UF, *FGridF, *FrbGridF, *GridPtrF, *GridRBPtrF, pv_mass, strange_mass,      pv_mass, -1.0, 1, M5, b, c);

  ConjugateGradient<FermionField>      ActionCG(ActionStoppingCondition,MaxCGIterations);
  ConjugateGradient<FermionField>  DerivativeCG(DerivativeStoppingCondition,MaxCGIterations);
#ifdef MIXED_PRECISION
  const int MX_inner = 1000;
  // Mixed precision EOFA
  LinearOperatorEOFAD Strange_LinOp_L (Strange_Op_L);
  LinearOperatorEOFAD Strange_LinOp_R (Strange_Op_R);
  LinearOperatorEOFAF Strange_LinOp_LF(Strange_Op_LF);
  LinearOperatorEOFAF Strange_LinOp_RF(Strange_Op_RF);

  MxPCG_EOFA ActionCGL(ActionStoppingCondition,
		       MX_inner,
		       MaxCGIterations,
		       GridPtrF,
		       FrbGridF,
		       Strange_Op_LF,Strange_Op_L,
		       Strange_LinOp_LF,Strange_LinOp_L);

  MxPCG_EOFA DerivativeCGL(DerivativeStoppingCondition,
			   MX_inner,
			   MaxCGIterations,
			   GridPtrF,
			   FrbGridF,
			   Strange_Op_LF,Strange_Op_L,
			   Strange_LinOp_LF,Strange_LinOp_L);
  
  MxPCG_EOFA ActionCGR(ActionStoppingCondition,
		       MX_inner,
		       MaxCGIterations,
		       GridPtrF,
		       FrbGridF,
		       Strange_Op_RF,Strange_Op_R,
		       Strange_LinOp_RF,Strange_LinOp_R);
  
  MxPCG_EOFA DerivativeCGR(DerivativeStoppingCondition,
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
	 OFRp, true);
#else
  ExactOneFlavourRatioPseudoFermionAction<FermionImplPolicy> 
    EOFA(Strange_Op_L, Strange_Op_R, 
	 ActionCG, ActionCG,
	 DerivativeCG, DerivativeCG,
	 OFRp, true);
#endif
  Level1.push_back(&EOFA);

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
  std::vector<FermionAction *> Numerators;
  std::vector<FermionAction *> Denominators;
  std::vector<TwoFlavourEvenOddRatioPseudoFermionAction<FermionImplPolicy> *> Quotients;
  std::vector<MxPCG *> ActionMPCG;
  std::vector<MxPCG *> MPCG;
  std::vector<FermionActionF *> DenominatorsF;
  std::vector<LinearOperatorD *> LinOpD;
  std::vector<LinearOperatorF *> LinOpF; 

  for(int h=0;h<n_hasenbusch+1;h++){

    std::cout << GridLogMessage << " 2f quotient Action  "<< light_num[h] << " / " << light_den[h]<< std::endl;

    Numerators.push_back  (new FermionAction(U,*FGrid,*FrbGrid,*GridPtr,*GridRBPtr,light_num[h],M5,b,c, Params));
    Denominators.push_back(new FermionAction(U,*FGrid,*FrbGrid,*GridPtr,*GridRBPtr,light_den[h],M5,b,c, Params));

#ifdef MIXED_PRECISION
    ////////////////////////////////////////////////////////////////////////////
    // Mixed precision CG for 2f force
    ////////////////////////////////////////////////////////////////////////////

    DenominatorsF.push_back(new FermionActionF(UF,*FGridF,*FrbGridF,*GridPtrF,*GridRBPtrF,light_den[h],M5,b,c, ParamsF));
    LinOpD.push_back(new LinearOperatorD(*Denominators[h]));
    LinOpF.push_back(new LinearOperatorF(*DenominatorsF[h]));

    MPCG.push_back(new MxPCG(DerivativeStoppingCondition,
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

  std::cout << GridLogMessage << " Running the HMC "<< std::endl;
  TheHMC.Run();  // no smearing

  Grid_finalize();
} // main



