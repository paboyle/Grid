/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./tests/Test_hmc_shGordonAction.cc

Copyright (C) 2018

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

namespace Grid {
class ScalarActionParameters : Serializable {
 public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(ScalarActionParameters,
    double, mass_squared,
    double, g);

    template <class ReaderClass >
  ScalarActionParameters(Reader<ReaderClass>& Reader){
    read(Reader, "ScalarAction", *this);
  }

};
}




int main(int argc, char **argv) {
  using namespace Grid;
  using namespace Grid::QCD;

  typedef Grid::JSONReader       Serialiser;

  Grid_init(&argc, &argv);
  int threads = GridThread::GetThreads();
  // here make a routine to print all the relevant information on the run
  std::cout << GridLogMessage << "Grid is setup to use " << threads << " threads" << std::endl;

   // Typedefs to simplify notation
  typedef ScalarGenericHMCRunner HMCWrapper;  // Uses the default minimum norm, real scalar fields
  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  HMCWrapper TheHMC;
  TheHMC.ReadCommandLine(argc, argv);

  if (TheHMC.ParameterFile.empty()){
    std::cout << "Input file not specified."
              << "Use --ParameterFile option in the command line.\nAborting" 
              << std::endl;
    exit(1);
  }
  Serialiser Reader(TheHMC.ParameterFile);
  // Grid from the command line
  constexpr int Ndimensions = 2;
  GridModule ScalarGrid;
  if (GridDefaultLatt().size() != Ndimensions){
    std::cout << "Incorrect dimension of the grid\n. Expected dim="<< Ndimensions << std::endl;
    exit(1);
  }
  if (GridDefaultMpi().size() != Ndimensions){
    std::cout << "Incorrect dimension of the mpi grid\n. Expected dim="<< Ndimensions << std::endl;
    exit(1);
  }
  ScalarGrid.set_full(new GridCartesian(GridDefaultLatt(),GridDefaultSimd(Ndimensions, vComplex::Nsimd()),GridDefaultMpi()));
  ScalarGrid.set_rb(new GridRedBlackCartesian(ScalarGrid.get_full()));
  TheHMC.Resources.AddGrid("scalar", ScalarGrid);
  std::cout << "Lattice size : " << GridDefaultLatt() << std::endl;

  CheckpointerParameters CPparams(Reader);
  TheHMC.Resources.LoadBinaryCheckpointer(CPparams);

  RNGModuleParameters RNGpar(Reader);
  TheHMC.Resources.SetRNGSeeds(RNGpar);

  // Some online observable measurements
  typedef ExpScalarMod<HMCWrapper::ImplPolicy> ExpObs;
  ExpScalarParameters ExpParams(Reader);
  TheHMC.Resources.AddObservable<ExpObs>(ExpParams);
  ///////////////////////////////////////////

  // Real Scalar sh-Gordon action
  ScalarActionParameters SPar(Reader);
  shGordonActionR Saction(SPar.mass_squared, SPar.g);

  // Collect actions
  ActionLevel<shGordonActionR::Field, ScalarFields> Level1(1);
  Level1.push_back(&Saction);

  TheHMC.TheAction.push_back(Level1);
  /////////////////////////////////////////////////////////////

  TheHMC.Parameters.initialize(Reader);
  TheHMC.Run(); 

  Grid_finalize();

} // main


