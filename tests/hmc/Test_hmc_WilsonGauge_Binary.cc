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

namespace Grid {
  namespace QCD {

    //Change here the type of reader
    typedef Grid::XmlReader InputFileReader; 


    class HMCRunnerParameters : Serializable {
    public:
      GRID_SERIALIZABLE_CLASS_MEMBERS(HMCRunnerParameters,
        double, beta,
        int, MDsteps,
        double, TrajectoryLength,
        //int, SaveInterval,
        //std::string, format,
        //std::string, conf_prefix,
        //std::string, rng_prefix,
        std::string, serial_seeds,
        std::string, parallel_seeds,
        );

      HMCRunnerParameters() {}
    };
  }
}

int main(int argc, char **argv) {
  using namespace Grid;
  using namespace Grid::QCD;

  Grid_init(&argc, &argv);
  int threads = GridThread::GetThreads();

   // Typedefs to simplify notation
  typedef GenericHMCRunner<MinimumNorm2> HMCWrapper;  // Uses the default minimum norm

  // here make a routine to print all the relevant information on the run
  std::cout << GridLogMessage << "Grid is setup to use " << threads << " threads" << std::endl;

  //////////////////////////////////////////////////////////////
  // Input file section 
  // make input file name general
  // now working with the text reader but I should drop this support
  // i need a structured format where every object is able
  // to locate the required data: XML, JSON, YAML.
  InputFileReader Reader("input.wilson_gauge.params.xml");
  HMCRunnerParameters HMCPar;
  read(Reader, "HMC", HMCPar);

  // Seeds for the random number generators
  // generalise, ugly now
  std::vector<int> SerSeed = strToVec<int>(HMCPar.serial_seeds);
  std::vector<int> ParSeed = strToVec<int>(HMCPar.parallel_seeds);

  HMCWrapper TheHMC;
  TheHMC.Resources.AddFourDimGrid("gauge");
  TheHMC.Resources.LoadBinaryCheckpointer(Reader);

  /////////////////////////////////////////////////////////////
  // Collect actions, here use more encapsulation
  // need wrappers of the fermionic classes 
  // that have a complex construction

  // Gauge action
  WilsonGaugeActionR Waction(HMCPar.beta);

  ActionLevel<HMCWrapper::Field> Level1(1);
  Level1.push_back(&Waction);
  TheHMC.TheAction.push_back(Level1);
  /////////////////////////////////////////////////////////////

  // Construct observables 
  PlaquetteLogger<HMCWrapper::ImplPolicy> PlaqLog("Plaquette");
  TheHMC.ObservablesList.push_back(&PlaqLog);
  //////////////////////////////////////////////

  // Fill resources
  // here we can simplify a lot if the input file is structured
  // just pass the input file reader 
  TheHMC.Resources.AddRNGSeeds(SerSeed, ParSeed);
  TheHMC.MDparameters.set(HMCPar.MDsteps, HMCPar.TrajectoryLength);

  // eventually smearing here
  // ...
  ////////////////////////////////////////////////////////////////

  TheHMC.ReadCommandLine(argc, argv); // these must be parameters from file
  TheHMC.Run();  // no smearing

  Grid_finalize();

} // main
