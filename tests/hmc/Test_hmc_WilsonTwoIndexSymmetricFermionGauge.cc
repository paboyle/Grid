/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./tests/Test_hmc_WilsonAdjointFermionGauge.cc

Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: neo <cossu@post.kek.jp>

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
#include "Grid/Grid.h"

int main(int argc, char **argv) {
  using namespace Grid;
   ;

  // Here change the allowed (higher) representations
  typedef Representations< FundamentalRepresentation, TwoIndexSymmetricRepresentation > TheRepresentations;

  Grid_init(&argc, &argv);
  int threads = GridThread::GetThreads();
  // here make a routine to print all the relevant information on the run
  std::cout << GridLogMessage << "Grid is setup to use " << threads << " threads" << std::endl;

   // Typedefs to simplify notation
  typedef GenericHMCRunnerHirep<TheRepresentations, MinimumNorm2> HMCWrapper;

  typedef WilsonTwoIndexSymmetricImplR FermionImplPolicy; // gauge field implemetation for the pseudofermions
  typedef WilsonTwoIndexSymmetricFermionR FermionAction; // type of lattice fermions (Wilson, DW, ...)
  typedef typename FermionAction::FermionField FermionField;

  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  HMCWrapper TheHMC;

  // Grid from the command line
  TheHMC.Resources.AddFourDimGrid("gauge");
  // Possibile to create the module by hand 
  // hardcoding parameters or using a Reader


  // Checkpointer definition
  CheckpointerParameters CPparams;  
  CPparams.config_prefix = "ckpoint_lat";
  CPparams.rng_prefix = "ckpoint_rng";
  CPparams.saveInterval = 5;
  CPparams.format = "IEEE64BIG";
  
  TheHMC.Resources.LoadNerscCheckpointer(CPparams);

  RNGModuleParameters RNGpar;
  RNGpar.serial_seeds = "1 2 3 4 5";
  RNGpar.parallel_seeds = "6 7 8 9 10";
  TheHMC.Resources.SetRNGSeeds(RNGpar);

  // Construct observables
  typedef PlaquetteMod<HMCWrapper::ImplPolicy> PlaqObs;
  TheHMC.Resources.AddObservable<PlaqObs>();
  //////////////////////////////////////////////

  /////////////////////////////////////////////////////////////
  // Collect actions, here use more encapsulation
  // need wrappers of the fermionic classes 
  // that have a complex construction
  // standard
  RealD beta = 2.25 ;
  WilsonGaugeActionR Waction(beta);
    
  auto GridPtr = TheHMC.Resources.GetCartesian();
  auto GridRBPtr = TheHMC.Resources.GetRBCartesian();

  // temporarily need a gauge field
  TwoIndexSymmetricRepresentation::LatticeField U(GridPtr);

  Real mass = -0.95;

  // Can we define an overloaded operator that does not need U and initialises
  // it with zeroes?
  FermionAction FermOp(U, *GridPtr, *GridRBPtr, mass);

  ConjugateGradient<FermionField> CG(1.0e-8, 2000, false);

  TwoFlavourPseudoFermionAction<FermionImplPolicy> Nf2(FermOp, CG, CG);

  // Set smearing (true/false), default: false
  Nf2.is_smeared = false;


    // Collect actions
  ActionLevel<LatticeGaugeField, TheRepresentations > Level1(1);
  Level1.push_back(&Nf2);

  ActionLevel<LatticeGaugeField, TheRepresentations > Level2(4);
  Level2.push_back(&Waction);

  TheHMC.TheAction.push_back(Level1);
  TheHMC.TheAction.push_back(Level2);

  // HMC parameters are serialisable 
  TheHMC.Parameters.MD.MDsteps = 20;
  TheHMC.Parameters.MD.trajL   = 1.0;

  TheHMC.ReadCommandLine(argc, argv); // these can be parameters from file
  TheHMC.Run();  // no smearing
  // TheHMC.Run(SmearingPolicy); // for smearing

  Grid_finalize();

} // main


