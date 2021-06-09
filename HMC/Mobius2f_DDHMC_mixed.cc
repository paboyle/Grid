/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

nnSource file: 

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

template<class Impl>
class DomainLocalTwoFlavourEvenOddRatioPseudoFermionAction
  : public TwoFlavourEvenOddRatioPseudoFermionAction<Impl>
{
public:
  INHERIT_IMPL_TYPES(Impl);
  Coordinate Block;
  DomainDecomposition Domains;
  DomainLocalTwoFlavourEvenOddRatioPseudoFermionAction(FermionOperator<Impl>  &_NumOp, 
						       FermionOperator<Impl>  &_DenOp, 
						       OperatorFunction<FermionField> & DS,
						       OperatorFunction<FermionField> & AS,
						       OperatorFunction<FermionField> & HS,
						       Coordinate &_Block ) :
    Block(_Block),
    Domains(_Block),
    TwoFlavourEvenOddRatioPseudoFermionAction<Impl>(_NumOp,_DenOp,DS,AS,HS)
    {};
  virtual void refreshRestrict(FermionField &eta)
  {
    Domains.ProjectDomain(eta,1);
    DumpSliceNorm("refresh Restrict eta",eta);
  };
};

#define MIXED_PRECISION

NAMESPACE_END(Grid);

int main(int argc, char **argv)
{
  using namespace Grid;

  Grid_init(&argc, &argv);
  int threads = GridThread::GetThreads();
  // here make a routine to print all the relevant information on the run
  std::cout << GridLogMessage << "Grid is setup to use " << threads << " threads" << std::endl;

   // Typedefs to simplify notation
  typedef WilsonImplR FimplD;
  typedef WilsonImplF FimplF;
  typedef FermionOperator<FimplF> FermionOperatorF;
  typedef FermionOperator<FimplD> FermionOperatorD;
  typedef MobiusFermionR FermionActionD;
  typedef MobiusFermionF FermionActionF;
  typedef DirichletFermionOperator<WilsonImplR> DirichletFermionD;
  typedef DirichletFermionOperator<WilsonImplF> DirichletFermionF;

  typedef MobiusEOFAFermionR FermionEOFAAction;
  typedef typename FermionActionD::FermionField FermionFieldD;
  typedef typename FermionActionF::FermionField FermionFieldF;
  
  typedef SchurDiagMooeeOperator<FermionOperator<FimplF>,FermionFieldF> LinearOperatorF;
  typedef SchurDiagMooeeOperator<FermionOperator<FimplD>,FermionFieldD> LinearOperatorD;
  typedef SchurDiagMooeeDagOperator<FermionOperator<FimplF>,FermionFieldF> LinearOperatorDagF;
  typedef SchurDiagMooeeDagOperator<FermionOperator<FimplD>,FermionFieldD> LinearOperatorDagD;
  
  typedef Grid::XmlReader       Serialiser;
  
  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  IntegratorParameters MD;
  //  typedef GenericHMCRunner<LeapFrog> HMCWrapper; 
  //  MD.name    = std::string("Leap Frog");
  //  typedef GenericHMCRunner<ForceGradient> HMCWrapper; 
  //  MD.name    = std::string("Force Gradient");
  typedef GenericHMCRunner<MinimumNorm2> HMCWrapper; 
  MD.name    = std::string("MinimumNorm2");
  MD.MDsteps = 4; // dH = 0.08
  //  MD.MDsteps = 3; // dH = 0.8
  MD.trajL   = 1.0;
  
  HMCparameters HMCparams;
  HMCparams.StartTrajectory  = 48;
  HMCparams.Trajectories     = 20;
  HMCparams.NoMetropolisUntil=  0;
  //  "[HotStart, ColdStart, TepidStart, CheckpointStart]\n";
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

  RNGModuleParameters RNGpar;
  RNGpar.serial_seeds = "1 2 3 4 5";
  RNGpar.parallel_seeds = "6 7 8 9 10";
  TheHMC.Resources.SetRNGSeeds(RNGpar);

  // Momentum Dirichlet
  Coordinate Block({0,0,0,24});
  
  TheHMC.Resources.SetMomentumFilter(new DDHMCFilter<WilsonImplR::Field>(Block));
  // Construct observables
  // here there is too much indirection 
  typedef PlaquetteMod<HMCWrapper::ImplPolicy> PlaqObs;
  TheHMC.Resources.AddObservable<PlaqObs>();
  //////////////////////////////////////////////

  const int Ls      = 16;
  Real beta         = 2.13;
  //  Real light_mass   = 0.04;
  Real light_mass   = 0.01;
  Real pv_mass      = 1.0;
  RealD M5  = 1.8;
  RealD b   = 1.0; 
  RealD c   = 0.0;

  std::vector<Real> hasenbusch({ 0.1, 0.4, 0.7 });
  
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
  FermionActionD::ImplParams Params(boundary);
  FermionActionD::ImplParams DirichletParams(boundary);
  DirichletParams.locally_periodic=true;
  
  double ActionStoppingCondition     = 1e-10;
  double DerivativeStoppingCondition = 1e-10;
  //  double BoundaryDerivativeStoppingCondition = 1e-6;
  double BoundaryDerivativeStoppingCondition = 1e-10;
  double MaxCGIterations = 30000;

  ////////////////////////////////////
  // Collect actions
  ////////////////////////////////////
  ActionLevel<HMCWrapper::Field> Level1(1);
  ActionLevel<HMCWrapper::Field> Level2(3);
  ActionLevel<HMCWrapper::Field> Level3(8);

  ConjugateGradient<FermionFieldD>      ActionCG(ActionStoppingCondition,MaxCGIterations);
  ConjugateGradient<FermionFieldD>  DerivativeCG(DerivativeStoppingCondition,MaxCGIterations);

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

  /////////////////////////////////////////////////
  // These are consumed/owned by the Dirichlet wrappers
  /////////////////////////////////////////////////
  std::vector<FermionActionD *> DNumeratorsD;
  std::vector<FermionActionF *> DNumeratorsF;
  std::vector<FermionActionD *> DDenominatorsD;
  std::vector<FermionActionF *> DDenominatorsF;

  /////////////////////////////////////////////////
  // Dirichlet wrappers
  /////////////////////////////////////////////////
  std::vector<DirichletFermionD *> DirichletNumeratorsD;
  std::vector<DirichletFermionF *> DirichletNumeratorsF;
  std::vector<DirichletFermionD *> DirichletDenominatorsD;
  std::vector<DirichletFermionF *> DirichletDenominatorsF;
  
  std::vector<DomainLocalTwoFlavourEvenOddRatioPseudoFermionAction<FimplD> *> Quotients;

  typedef MixedPrecisionConjugateGradientOperatorFunction<FermionOperatorD,
							  FermionOperatorF,
							  LinearOperatorD,
							  LinearOperatorF> MxPCG;
  std::vector<MxPCG *> ActionMPCG;
  std::vector<MxPCG *> MPCG;
  std::vector<LinearOperatorD *> LinOpD;
  std::vector<LinearOperatorF *> LinOpF; 

  const int MX_inner = 1000;
  const RealD MX_tol = 1.0e-8;

  for(int h=0;h<n_hasenbusch+1;h++){

    std::cout << GridLogMessage << " 2f quotient Action  "<< light_num[h] << " / " << light_den[h]<< std::endl;

    DNumeratorsD.push_back (new FermionActionD(U,*FGrid,*FrbGrid,*GridPtr,*GridRBPtr,light_num[h],M5,b,c, DirichletParams));
    DNumeratorsF.push_back (new FermionActionF(UF,*FGridF,*FrbGridF,*GridPtrF,*GridRBPtrF,light_num[h],M5,b,c, DirichletParams));
    DDenominatorsD.push_back(new FermionActionD(U,*FGrid,*FrbGrid,*GridPtr,*GridRBPtr,light_den[h],M5,b,c, DirichletParams));
    DDenominatorsF.push_back(new FermionActionF(UF,*FGridF,*FrbGridF,*GridPtrF,*GridRBPtrF,light_den[h],M5,b,c, DirichletParams));

    DirichletNumeratorsD.push_back  (new  DirichletFermionD(*DNumeratorsD[h],Block));
    DirichletNumeratorsF.push_back  (new  DirichletFermionF(*DNumeratorsF[h],Block));
    DirichletDenominatorsD.push_back(new  DirichletFermionD(*DDenominatorsD[h],Block));
    DirichletDenominatorsF.push_back(new  DirichletFermionF(*DDenominatorsF[h],Block));

    // Dirichlet Schur even odd MpsDagMpc operators on local domains
    LinOpD.push_back(new LinearOperatorD(*DirichletDenominatorsD[h]));
    LinOpF.push_back(new LinearOperatorF(*DirichletDenominatorsF[h]));

    // Derivative
    MPCG.push_back(new MxPCG(DerivativeStoppingCondition,MX_tol,
			     MX_inner,
			     MaxCGIterations,
			     FrbGridF,
			     *DirichletDenominatorsF[h],*DirichletDenominatorsD[h],
			     *LinOpF[h], *LinOpD[h]) );

    // Action
    ActionMPCG.push_back(new MxPCG(ActionStoppingCondition,MX_tol,
				   MX_inner,
				   MaxCGIterations,
				   FrbGridF,
				   *DirichletDenominatorsF[h],*DirichletDenominatorsD[h],
				   *LinOpF[h], *LinOpD[h]) );
    
    ////////////////////////////////////////////////////////////////////////////
    // Standard CG for 2f force
    ////////////////////////////////////////////////////////////////////////////
    Quotients.push_back   (new
			   DomainLocalTwoFlavourEvenOddRatioPseudoFermionAction<FimplD>
			   (*DirichletNumeratorsD[h],
			    *DirichletDenominatorsD[h],
			    *MPCG[h],
			    *ActionMPCG[h],
			    ActionCG,Block));

    Level2.push_back(Quotients[h]);
  }

  /////////////////////////////////////////////////////////////
  // Boundary action
  /////////////////////////////////////////////////////////////

  int l_idx = 0;
  int pv_idx = n_hasenbusch;
  RealD h_mass = 0.012;
  std::cout << GridLogMessage<<" Boundary action masses " <<light_num[l_idx]<<" / "<<light_den[pv_idx]<<std::endl;


  // OmegaBar cross domain boundary and is used in Boundary operator, so no locally_periodic hack in the boundary det
  // Dirichlet is applied in gauge link only. OmegaBar solve is too expensive. Monitor cost.
  FermionActionD    PeriNumeratorD  (U,*FGrid,*FrbGrid,*GridPtr,*GridRBPtr,light_num[pv_idx],M5,b,c, Params);
  FermionActionF    PeriNumeratorF  (UF,*FGridF,*FrbGridF,*GridPtrF,*GridRBPtrF,light_num[pv_idx],M5,b,c, Params);
  FermionActionD    DirichletNumeratorDD(U,*FGrid,*FrbGrid,*GridPtr,*GridRBPtr,light_num[pv_idx],M5,b,c, Params);
  FermionActionF    DirichletNumeratorFF(UF,*FGridF,*FrbGridF,*GridPtrF,*GridRBPtrF,light_num[pv_idx],M5,b,c, Params);
  DirichletFermionD DirichletNumeratorD  (DirichletNumeratorDD,Block);
  DirichletFermionF DirichletNumeratorF  (DirichletNumeratorFF,Block);

  FermionActionD    PeriDenominatorD(U,*FGrid,*FrbGrid,*GridPtr,*GridRBPtr,light_den[l_idx] ,M5,b,c, Params);
  FermionActionF    PeriDenominatorF(UF,*FGridF,*FrbGridF,*GridPtrF,*GridRBPtrF,light_den[l_idx] ,M5,b,c, Params);
  FermionActionD    DirichletDenominatorDD(U,*FGrid,*FrbGrid,*GridPtr,*GridRBPtr,light_den[l_idx] ,M5,b,c, Params);
  FermionActionF    DirichletDenominatorFF(UF,*FGridF,*FrbGridF,*GridPtrF,*GridRBPtrF,light_den[l_idx] ,M5,b,c, Params);
  DirichletFermionD DirichletDenominatorD(DirichletDenominatorDD,Block);
  DirichletFermionF DirichletDenominatorF(DirichletDenominatorFF,Block);

  FermionActionD    PeriHasenD  (U,*FGrid,*FrbGrid,*GridPtr,*GridRBPtr,h_mass ,M5,b,c, Params);
  FermionActionF    PeriHasenF  (UF,*FGridF,*FrbGridF,*GridPtrF,*GridRBPtrF,h_mass,M5,b,c, Params);
  FermionActionD    DHasenD(U,*FGrid,*FrbGrid,*GridPtr,*GridRBPtr,h_mass,M5,b,c, Params);
  FermionActionF    DHasenF(UF,*FGridF,*FrbGridF,*GridPtrF,*GridRBPtrF,h_mass,M5,b,c, Params);
  DirichletFermionD DirichletHasenD(DHasenD,Block);
  DirichletFermionF DirichletHasenF(DHasenF,Block);
  
  SchurFactoredFermionOperator<FimplD,FimplF> BoundaryNumerator(PeriNumeratorD,PeriNumeratorF,
								DirichletNumeratorD,DirichletNumeratorF,
								Block);

  SchurFactoredFermionOperator<FimplD,FimplF> BoundaryDenominator(PeriDenominatorD,PeriDenominatorF,
								  DirichletDenominatorD,DirichletDenominatorF,
								  Block);

  SchurFactoredFermionOperator<FimplD,FimplF> BoundaryHasen(PeriHasenD,PeriHasenF,
							    DirichletHasenD,DirichletHasenF,
							    Block);

  std::cout << GridLogMessage << " Boundary NO ratio "<< std::endl;
  Level1.push_back(new
		   DomainDecomposedBoundaryTwoFlavourPseudoFermion<FimplD,FimplF>
		   (BoundaryDenominator,
		    BoundaryDerivativeStoppingCondition,ActionStoppingCondition,MX_tol));
  Level1.push_back(new
		   DomainDecomposedBoundaryTwoFlavourBosonPseudoFermion<FimplD,FimplF>
		   (BoundaryNumerator,
		    BoundaryDerivativeStoppingCondition,ActionStoppingCondition,MX_tol));
  /*
  Level1.push_back(new
		   DomainDecomposedBoundaryTwoFlavourRatioPseudoFermion<FimplD,FimplF>
		   (BoundaryNumerator,
		    BoundaryDenominator,
		    BoundaryDerivativeStoppingCondition,ActionStoppingCondition));
  */

  /////////////////////////////////////////////////////////////
  // Gauge action
  /////////////////////////////////////////////////////////////
  Level3.push_back(&GaugeAction);
  TheHMC.TheAction.push_back(Level1);
  TheHMC.TheAction.push_back(Level2);
  TheHMC.TheAction.push_back(Level3);
  std::cout << GridLogMessage << " Action complete "<< std::endl;

  /////////////////////////////////////////////////////////////
  // HMC parameters are serialisable
  std::cout << GridLogMessage << " Running the HMC "<< std::endl;
  TheHMC.Run();  // no smearing

  Grid_finalize();
} // main



