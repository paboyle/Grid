/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./tests/Test_hmc_WilsonFermionGauge.cc

Copyright (C) 2016

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

See the full license in the file "LICENSE" in the top level distribution directory
*************************************************************************************/
/*  END LEGAL */
#include <Grid/Grid.h>

int main(int argc, char **argv) {
  using namespace Grid;
   ;

  Grid_init(&argc, &argv);
  int threads = GridThread::GetThreads();
  // here make a routine to print all the relevant information on the run
  std::cout << GridLogMessage << "Grid is setup to use " << threads << " threads" << std::endl;

   // Typedefs to simplify notation
  typedef ScalarGenericHMCRunner HMCWrapper;  // Uses the default minimum norm, real scalar fields
  //typedef Representations<EmptyRep<typename ScalarMatrixImplTypes<vComplex, 3>::Field> > ScalarMatrixFields;
  //typedef HMCWrapperTemplate<ScalarMatrixImplTypes<vComplex, 3>, MinimumNorm2, ScalarMatrixFields> HMCWrapper;
  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  HMCWrapper TheHMC;

  // Grid from the command line
  GridModule ScalarGrid;
  ScalarGrid.set_full( SpaceTimeGrid::makeFourDimGrid(
        GridDefaultLatt(), GridDefaultSimd(Nd, vReal::Nsimd()),
        GridDefaultMpi()));
  ScalarGrid.set_rb(SpaceTimeGrid::makeFourDimRedBlackGrid(ScalarGrid.get_full()));
  TheHMC.Resources.AddGrid("scalar", ScalarGrid);
  // Possibile to create the module by hand 
  // hardcoding parameters or using a Reader

  // Checkpointer definition
  CheckpointerParameters CPparams;  
  CPparams.config_prefix = "ckpoint_scalar_lat";
  CPparams.rng_prefix = "ckpoint_scalar_rng";
  CPparams.saveInterval = 5;
  CPparams.format = "IEEE64BIG";
  
  TheHMC.Resources.LoadBinaryCheckpointer(CPparams);

  RNGModuleParameters RNGpar;
  RNGpar.serial_seeds = "1 2 3 4 5";
  RNGpar.parallel_seeds = "6 7 8 9 10";
  TheHMC.Resources.SetRNGSeeds(RNGpar);

  //////////////////////////////////////////////

  /////////////////////////////////////////////////////////////
  // Collect actions, here use more encapsulation
  // need wrappers of the fermionic classes 
  // that have a complex construction
  // standard

  // Real Scalar action
  ScalarActionR Saction(0.11,0.);
  
  //typedef ScalarAction<ScalarMatrixImplTypes<vComplex, 3>> ScalarMatrixActionR;
  //ScalarMatrixActionR Saction(0.11,0.);

  // Collect actions
  ActionLevel<ScalarActionR::Field, ScalarFields> Level1(1);
  Level1.push_back(&Saction);

  TheHMC.TheAction.push_back(Level1);
  /////////////////////////////////////////////////////////////


  // HMC parameters are serialisable 
  TheHMC.Parameters.MD.MDsteps = 20;
  TheHMC.Parameters.MD.trajL   = 1.0;

  TheHMC.ReadCommandLine(argc, argv); // these can be parameters from file
  TheHMC.Run(); 

  Grid_finalize();

} // main


