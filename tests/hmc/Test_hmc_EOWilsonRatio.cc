    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_hmc_EOWilsonRatio.cc

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: paboyle <paboyle@ph.ed.ac.uk>

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

    See the full license in the file "LICENSE" in the top level distribution directory
    *************************************************************************************/
    /*  END LEGAL */
#include <Grid/Grid.h>



int main(int argc, char **argv) {
  using namespace Grid;
  using namespace Grid::QCD;

  Grid_init(&argc, &argv);
  int threads = GridThread::GetThreads();
  // here make a routine to print all the relevant information on the run
  std::cout << GridLogMessage << "Grid is setup to use " << threads << " threads" << std::endl;

   // Typedefs to simplify notation
  typedef GenericHMCRunner<MinimumNorm2> HMCWrapper;  // Uses the default minimum norm
  typedef WilsonImplR FermionImplPolicy;
  typedef WilsonFermionR FermionAction;
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
  
  TheHMC.Resources.LoadBinaryCheckpointer(CPparams);

  RNGModuleParameters RNGpar;
  RNGpar.SerialSeed = {1,2,3,4,5};
  RNGpar.ParallelSeed = {6,7,8,9,10};
  TheHMC.Resources.SetRNGSeeds(RNGpar);

  // Construct observables
  // here there is too much indirection 
  PlaquetteObsParameters PlPar;
  PlPar.output_prefix = "Plaquette";
  PlaquetteMod<HMCWrapper::ImplPolicy> PlaqModule(PlPar);
  TheHMC.Resources.AddObservable(&PlaqModule);
  //////////////////////////////////////////////

  /////////////////////////////////////////////////////////////
  // Collect actions, here use more encapsulation
  // need wrappers of the fermionic classes 
  // that have a complex construction
  // standard
  RealD beta = 5.6 ;
  WilsonGaugeActionR Waction(beta);
    
  auto GridPtr = TheHMC.Resources.GetCartesian();
  auto GridRBPtr = TheHMC.Resources.GetRBCartesian();

  // temporarily need a gauge field
  LatticeGaugeField U(GridPtr);

  Real mass = -0.77;
  RealD pv  =  0.0;

  // Can we define an overloaded operator that does not need U and initialises
  // it with zeroes?
  FermionAction DenOp(U, *GridPtr, *GridRBPtr, mass);
  FermionAction NumOp(U, *GridPtr, *GridRBPtr, pv);

  ConjugateGradient<FermionField> CG(1.0e-8, 2000);

  TwoFlavourEvenOddRatioPseudoFermionAction<FermionImplPolicy> Nf2(NumOp, DenOp,CG,CG);

  // Set smearing (true/false), default: false
  Nf2.is_smeared = false;


    // Collect actions
  ActionLevel<HMCWrapper::Field> Level1(1);
  Level1.push_back(&Nf2);

  ActionLevel<HMCWrapper::Field> Level2(4);
  Level2.push_back(&Waction);

  TheHMC.TheAction.push_back(Level1);
  TheHMC.TheAction.push_back(Level2);
  /////////////////////////////////////////////////////////////
  // HMC parameters are serialisable 
  TheHMC.Parameters.MD.MDsteps = 20;
  TheHMC.Parameters.MD.trajL   = 1.0;

  TheHMC.ReadCommandLine(argc, argv); // these can be parameters from file
  TheHMC.Run();  // no smearing

  Grid_finalize();

} // main


