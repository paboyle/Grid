/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./tests/Test_hmc_WilsonFermionGauge.cc

Copyright (C) 2015

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

See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
/*  END LEGAL */
#include <Grid/Grid.h>

namespace Grid{
  struct ActionParameters: Serializable {
    GRID_SERIALIZABLE_CLASS_MEMBERS(ActionParameters,
				    double, beta)

    ActionParameters() = default;

    template <class ReaderClass >
    ActionParameters(Reader<ReaderClass>& Reader){
      read(Reader, "Action", *this);
    }

  };

}


int main(int argc, char **argv) {
  using namespace Grid;
  using namespace Grid::QCD;

  Grid_init(&argc, &argv);
  GridLogLayout();

  // Typedefs to simplify notation
  typedef GenericHMCRunner<MinimumNorm2> HMCWrapper;  // Uses the default minimum norm
  HMCWrapper TheHMC;
  typedef Grid::JSONReader       Serialiser;

  // Grid from the command line
  TheHMC.Resources.AddFourDimGrid("gauge");
  TheHMC.ReadCommandLine(argc, argv); // these can be parameters from file
  // Reader, file should come from command line
  if (TheHMC.ParameterFile.empty()){
    std::cout << "Input file not specified."
	      << "Use --ParameterFile option in the command line.\nAborting"
	      << std::endl;
    exit(1);
  }
  Serialiser Reader(TheHMC.ParameterFile);

  // Read parameters from input file
  ActionParameters WilsonPar(Reader);

  // Checkpointer definition
  CheckpointerParameters CPparams(Reader);
  //TheHMC.Resources.LoadNerscCheckpointer(CPparams);

  // Store metadata in the Scidac checkpointer
  TheHMC.Resources.LoadScidacCheckpointer(CPparams, WilsonPar);

  RNGModuleParameters RNGpar(Reader);
  TheHMC.Resources.SetRNGSeeds(RNGpar);

  // Construct observables
  // here there is too much indirection
  typedef PlaquetteMod<HMCWrapper::ImplPolicy> PlaqObs;
  typedef TopologicalChargeMod<HMCWrapper::ImplPolicy> QObs;
  TheHMC.Resources.AddObservable<PlaqObs>();
  TopologyObsParameters TopParams(Reader);
  TheHMC.Resources.AddObservable<QObs>(TopParams);
  //////////////////////////////////////////////

  /////////////////////////////////////////////////////////////
  // Collect actions, here use more encapsulation
  // need wrappers of the fermionic classes
  // that have a complex construction
  // standard
  WilsonGaugeActionR Waction(WilsonPar.beta);

  ActionLevel<HMCWrapper::Field> Level1(1);
  Level1.push_back(&Waction);
  //Level1.push_back(WGMod.getPtr());
  TheHMC.TheAction.push_back(Level1);
  /////////////////////////////////////////////////////////////

  // HMC parameters are serialisable
  TheHMC.Parameters.initialize(Reader);

  //TheHMC.Parameters.MD.MDsteps = 17;
  //TheHMC.Parameters.MD.trajL   = 1.0;

  TheHMC.Run();  // no smearing

  Grid_finalize();

} // main
