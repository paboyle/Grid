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

//#include <Grid/qcd/action/pseudofermion/TwoFlavourRatioEO4DPseudoFermion.h>
//#include <Grid/qcd/action/pseudofermion/TwoFlavourRatio4DPseudoFermion.h>


NAMESPACE_BEGIN(Grid);

  /*
   * Need a plan for gauge field update for mixed precision in HMC                      (2x speed up)
   *    -- Store the single prec action operator.
   *    -- Clone the gauge field from the operator function argument.
   *    -- Build the mixed precision operator dynamically from the passed operator and single prec clone.
   */

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
  typedef typename FermionAction::Impl_t Fimpl;
  typedef DirichletFermionOperator<WilsonImplR> DirichletFermion;

  typedef MobiusEOFAFermionR FermionEOFAAction;
  typedef typename FermionAction::FermionField FermionField;

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
  HMCparams.Trajectories     = 1;
  HMCparams.NoMetropolisUntil=  0;
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

  // TheHMC.Resources.SetMomentumFilter(new DirichletFilter<WilsonImplR::Field>(Block));
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

  std::vector<Real> hasenbusch({ 0.04, 0.2 });

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

  // These lines are unecessary if BC are all periodic
  std::vector<Complex> boundary = {1,1,1,-1};
  FermionAction::ImplParams Params(boundary);
  FermionAction::ImplParams DirichletParams(boundary);
  DirichletParams.locally_periodic=true;

  double ActionStoppingCondition     = 1e-10;
  double DerivativeStoppingCondition = 1e-10;
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
  std::vector<FermionAction *> PeriNumerators;
  std::vector<FermionAction *> PeriDenominators;

  std::vector<FermionAction *> DNumerators;
  std::vector<FermionAction *> DDenominators;
  std::vector<DirichletFermion *> DirichletNumerators;
  std::vector<DirichletFermion *> DirichletDenominators;
  std::vector<SchurFactoredFermionOperator<Fimpl> *> BoundaryNumerators;
  std::vector<SchurFactoredFermionOperator<Fimpl> *> BoundaryDenominators;

  std::vector<TwoFlavourEvenOddRatioPseudoFermionAction<Fimpl> *> Quotients;
  std::vector<DomainDecomposedBoundaryTwoFlavourRatioPseudoFermion<Fimpl> *> BoundaryQuotients;
  std::vector<DomainDecomposedBoundaryTwoFlavourBosonPseudoFermion<Fimpl> *> BoundaryNums;
  std::vector<DomainDecomposedBoundaryTwoFlavourPseudoFermion<Fimpl> *> BoundaryDens;
  for(int h=0;h<n_hasenbusch+1;h++){

    PeriNumerators.push_back  (new FermionAction(U,*FGrid,*FrbGrid,*GridPtr,*GridRBPtr,light_num[h],M5,b,c, Params));
    PeriDenominators.push_back(new FermionAction(U,*FGrid,*FrbGrid,*GridPtr,*GridRBPtr,light_den[h],M5,b,c, Params));
       DNumerators.push_back  (new FermionAction(U,*FGrid,*FrbGrid,*GridPtr,*GridRBPtr,light_num[h],M5,b,c, DirichletParams));
       DDenominators.push_back(new FermionAction(U,*FGrid,*FrbGrid,*GridPtr,*GridRBPtr,light_den[h],M5,b,c, DirichletParams));

    DirichletNumerators.push_back  (new DirichletFermion(*DNumerators[h],Block));
    DirichletDenominators.push_back(new DirichletFermion(*DDenominators[h],Block));

    BoundaryNumerators.push_back (new SchurFactoredFermionOperator<Fimpl>
				  (*PeriNumerators[h],
				   *DirichletNumerators[h],
				   ActionCG,ActionCG,
				   ActionCG,ActionCG,
				   Block));

    BoundaryDenominators.push_back (new SchurFactoredFermionOperator<Fimpl>
				  (*PeriDenominators[h],
				   *DirichletDenominators[h],
				   ActionCG,ActionCG,
				   ActionCG,ActionCG,
				   Block));
    
    ////////////////////////////////////////////////////////////////////////////
    // Standard CG for 2f force
    ////////////////////////////////////////////////////////////////////////////
    std::cout << GridLogMessage << " 2f quotient Action  "<< light_num[h] << " / " << light_den[h]<< std::endl;

    Quotients.push_back   (new
			   TwoFlavourEvenOddRatioPseudoFermionAction<Fimpl>
			   (*DirichletNumerators[h],
			    *DirichletDenominators[h],
			    ActionCG,
			    ActionCG));

    std::cout << GridLogMessage << " 2f Boundary quotient Action  "<< light_num[h] << " / " << light_den[h]<< std::endl;
    
    BoundaryQuotients.push_back(new
				DomainDecomposedBoundaryTwoFlavourRatioPseudoFermion<Fimpl>
				(*BoundaryNumerators[h],
				 *BoundaryDenominators[h]));

  }

  for(int h=0;h<n_hasenbusch+1;h++){
    Level1.push_back(Quotients[h]);
    Level1.push_back(BoundaryQuotients[h]);
    //    Level1.push_back(BoundaryNums[h]);
    //    Level1.push_back(BoundaryDens[h]);
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



