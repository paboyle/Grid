/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./tests/Test_hmc_WilsonFermionGauge.cc

Copyright (C) 2015

Author: Peter Boyle <pabobyle@ph.ed.ac.uk>
Author: neo <cossu@post.kek.jp>
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

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

namespace Grid {
  namespace QCD {

    //Change here the type of reader
    typedef Grid::TextReader InputFileReader; 


    class HMCRunnerParameters : Serializable {
    public:
      GRID_SERIALIZABLE_CLASS_MEMBERS(HMCRunnerParameters,
        double, beta,
        int, MDsteps,
        double, TrajectorLength,
        int, SaveInterval,
        std::string, format,
        std::string, conf_prefix,
        std::string, rng_prefix,
        std::string, serial_seeds,
        std::string, parallel_seeds,
        );

      HMCRunnerParameters() {}
    };
  }
}

int main(int argc, char **argv) {
  Grid_init(&argc, &argv);
   // Typedefs to simplify notation
  typedef BinaryHmcRunner<MinimumNorm2> HMCWrapper;// Uses the default minimum norm
  typedef WilsonGaugeActionR GaugeAction;

  int threads = GridThread::GetThreads();
  std::cout << GridLogMessage << "Grid is setup to use " << threads << " threads" << std::endl;

  //////////////////////////////////////////////////////////////
  // Input file section 
  // make input file name general
  HMCRunnerParameters HMCPar;
  InputFileReader Reader("input.wilson_gauge.params");
  read(Reader, "HMC", HMCPar);
  std::cout << GridLogMessage << HMCPar << std::endl;

  // Seeds for the random number generators
  // generalise
  std::vector<int> SerSeed = strToVec<int>(HMCPar.serial_seeds);
  std::vector<int> ParSeed = strToVec<int>(HMCPar.parallel_seeds);

    // Add observables
    // options for checkpointers
    // this can be moved outside the BuildTheAction
    //BinaryHmcCheckpointer
    //ILDGHmcCheckpointer
    //NerscHmcCheckpointer
  BinaryHmcCheckpointer<HMCWrapper::ImplPolicy> Checkpoint(HMCPar.conf_prefix, HMCPar.rng_prefix, 
                                                           HMCPar.SaveInterval, HMCPar.format);
    // Can implement also a specific function in the hmcrunner
    // AddCheckpoint (...) that takes the same parameters + a string/tag
    // defining the type of the checkpointer
    // with tags can be implemented by overloading and no ifs
    // Then force all checkpoint to have few common functions
    // return an object that is then passed to the Run function

  /////////////////////////////////////////////////////////////
  HMCWrapper TheHMC;
  TheHMC.Resources.AddFourDimGrid("gauge");

  /////////////////////////////////////////////////////////////
  // Collect actions, here use more encapsulation

  // Gauge action
  std::cout << GridLogMessage << "Beta: " << HMCPar.beta << std::endl;
  GaugeAction Waction(HMCPar.beta);

  ActionLevel<HMCWrapper::Field> Level1(1);
  Level1.push_back(&Waction);
  TheHMC.TheAction.push_back(Level1);
  /////////////////////////////////////////////////////////////

  // Construct observables 
  PlaquetteLogger<HMCWrapper::ImplPolicy> PlaqLog("Plaquette");
  TheHMC.ObservablesList.push_back(&PlaqLog);
  TheHMC.ObservablesList.push_back(&Checkpoint);
  //////////////////////////////////////////////

  // Fill resources
  TheHMC.Resources.AddRNGSeeds(SerSeed, ParSeed);
  TheHMC.MDparameters.set(HMCPar.MDsteps, HMCPar.TrajectorLength);

  // eventually smearing here
  ////////////////////////////////////////////////////////////////


  TheHMC.ReadCommandLine(argc, argv);
  TheHMC.Run(Checkpoint);  // no smearing

  Grid_finalize();
}
