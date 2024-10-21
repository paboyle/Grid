/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Copyright (C) 2023

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

#if Nc == 3
#include <Grid/qcd/smearing/GaugeConfigurationMasked.h>
#include <Grid/qcd/smearing/JacobianAction.h>
#endif

using namespace Grid;

int main(int argc, char **argv)
{
#if Nc != 3
#warning FTHMC2p1f_3GeV will not work for Nc != 3
  std::cout << "This program will currently only work for Nc == 3." << std::endl;
#else
  std::cout << std::setprecision(12);
  
  Grid_init(&argc, &argv);
  int threads = GridThread::GetThreads();
  // here make a routine to print all the relevant information on the run
  std::cout << GridLogMessage << "Grid is setup to use " << threads << " threads" << std::endl;

   // Typedefs to simplify notation
  typedef WilsonImplR FermionImplPolicy;
  typedef MobiusFermionD FermionAction;
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
  MD.MDsteps = 24;
  MD.trajL   = 1.0;

  HMCparameters HMCparams;
  HMCparams.StartTrajectory  = 0;
  HMCparams.Trajectories     = 200;
  HMCparams.NoMetropolisUntil=  20;
  // "[HotStart, ColdStart, TepidStart, CheckpointStart]\n";
  //  HMCparams.StartingType     =std::string("HotStart");
  HMCparams.StartingType     =std::string("ColdStart");
  //  HMCparams.StartingType     =std::string("CheckpointStart");
  HMCparams.MD = MD;
  HMCWrapper TheHMC(HMCparams);

  // Grid from the command line arguments --grid and --mpi
  TheHMC.Resources.AddFourDimGrid("gauge"); // use default simd lanes decomposition

  CheckpointerParameters CPparams;
  CPparams.config_prefix = "ckpoint_EODWF_lat";
  CPparams.smeared_prefix = "ckpoint_EODWF_lat_smr";
  CPparams.rng_prefix    = "ckpoint_EODWF_rng";
  CPparams.saveInterval  = 1;
  CPparams.saveSmeared   = true;
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
  Real beta         = 2.37;
  Real light_mass   = 0.0047;
  Real strange_mass = 0.0186;
  Real pv_mass      = 1.0;
  RealD M5  = 1.8;
  RealD b   = 1.0; // Scale factor one, Shamir
  RealD c   = 0.0;

  OneFlavourRationalParams OFRp;
  OFRp.lo       = 1.0e-2;
  OFRp.hi       = 64;
  OFRp.MaxIter  = 10000;
  OFRp.tolerance= 1.0e-10;
  OFRp.degree   = 14;
  OFRp.precision= 40;

  std::vector<Real> hasenbusch({ 0.05, 0.1, 0.25, 0.5 });

  auto GridPtr   = TheHMC.Resources.GetCartesian();
  auto GridRBPtr = TheHMC.Resources.GetRBCartesian();
  auto FGrid     = SpaceTimeGrid::makeFiveDimGrid(Ls,GridPtr);
  auto FrbGrid   = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,GridPtr);

  IwasakiGaugeActionR GaugeAction(beta);

  // temporarily need a gauge field
  LatticeGaugeField U(GridPtr);
  LatticeGaugeField Uhot(GridPtr);

  // These lines are unecessary if BC are all periodic
  std::vector<Complex> boundary = {1,1,1,-1};
  FermionAction::ImplParams Params(boundary);

  double StoppingCondition = 1e-10;
  double MaxCGIterations = 30000;
  ConjugateGradient<FermionField>  CG(StoppingCondition,MaxCGIterations);

  bool ApplySmearing = true;
  
  ////////////////////////////////////
  // Collect actions
  ////////////////////////////////////
  ActionLevel<HMCWrapper::Field> Level1(1);
  ActionLevel<HMCWrapper::Field> Level2(2);

  ////////////////////////////////////
  // Strange action
  ////////////////////////////////////

  MobiusEOFAFermionD Strange_Op_L (U , *FGrid , *FrbGrid , *GridPtr , *GridRBPtr , strange_mass, strange_mass, pv_mass, 0.0, -1, M5, b, c);
  MobiusEOFAFermionD Strange_Op_R (U , *FGrid , *FrbGrid , *GridPtr , *GridRBPtr , pv_mass, strange_mass,      pv_mass, -1.0, 1, M5, b, c);
  ExactOneFlavourRatioPseudoFermionAction<FermionImplPolicy> 
    EOFA(Strange_Op_L, Strange_Op_R, 
	 CG,
	 CG, CG,
	 CG, CG, 
	 OFRp, false);

  EOFA.is_smeared = ApplySmearing;
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

  std::vector<FermionAction *> Numerators;
  std::vector<FermionAction *> Denominators;
  std::vector<TwoFlavourEvenOddRatioPseudoFermionAction<FermionImplPolicy> *> Quotients;

  for(int h=0;h<n_hasenbusch+1;h++){
    std::cout << GridLogMessage << " 2f quotient Action  "<< light_num[h] << " / " << light_den[h]<< std::endl;
    Numerators.push_back  (new FermionAction(U,*FGrid,*FrbGrid,*GridPtr,*GridRBPtr,light_num[h],M5,b,c, Params));
    Denominators.push_back(new FermionAction(U,*FGrid,*FrbGrid,*GridPtr,*GridRBPtr,light_den[h],M5,b,c, Params));
    Quotients.push_back   (new TwoFlavourEvenOddRatioPseudoFermionAction<FermionImplPolicy>(*Numerators[h],*Denominators[h],CG,CG));
  }

  for(int h=0;h<n_hasenbusch+1;h++){
    Quotients[h]->is_smeared = ApplySmearing;
    Level1.push_back(Quotients[h]);
  }

  /////////////////////////////////////////////////////////////
  // lnDetJacobianAction
  /////////////////////////////////////////////////////////////
  double rho = 0.1;  // smearing parameter
  int Nsmear = 1;    // number of smearing levels - must be multiple of 2Nd
  int Nstep  = 8*Nsmear;    // number of smearing levels - must be multiple of 2Nd
  Smear_Stout<HMCWrapper::ImplPolicy> Stout(rho);
  SmearedConfigurationMasked<HMCWrapper::ImplPolicy> SmearingPolicy(GridPtr, Nstep, Stout);
  JacobianAction<HMCWrapper::ImplPolicy> Jacobian(&SmearingPolicy);
  if( ApplySmearing ) Level1.push_back(&Jacobian);
  std::cout << GridLogMessage << " Built the Jacobian "<< std::endl;


  /////////////////////////////////////////////////////////////
  // Gauge action
  /////////////////////////////////////////////////////////////
  GaugeAction.is_smeared = ApplySmearing;
  Level2.push_back(&GaugeAction);

  std::cout << GridLogMessage << " ************************************************"<< std::endl;
  std::cout << GridLogMessage << " Action complete -- NO FERMIONS FOR NOW -- FIXME"<< std::endl;
  std::cout << GridLogMessage << " ************************************************"<< std::endl;
  std::cout << GridLogMessage <<  std::endl;
  std::cout << GridLogMessage <<  std::endl;


  std::cout << GridLogMessage << " Running the FT HMC "<< std::endl;
  TheHMC.TheAction.push_back(Level1);
  TheHMC.TheAction.push_back(Level2);

  TheHMC.ReadCommandLine(argc,argv);  // params on CML or from param file
  TheHMC.initializeGaugeFieldAndRNGs(U);

  TheHMC.Run(SmearingPolicy); // for smearing

  Grid_finalize();
#endif
} // main



