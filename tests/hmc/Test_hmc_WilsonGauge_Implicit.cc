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

namespace Grid{
struct RMHMCActionParameters: Serializable {
  GRID_SERIALIZABLE_CLASS_MEMBERS(RMHMCActionParameters,
				  double, gauge_beta)

  template <class ReaderClass >
  RMHMCActionParameters(Reader<ReaderClass>& Reader){
    read(Reader, "Action", *this);
  }
};
}

int main(int argc, char **argv) {
  using namespace Grid;
//  using namespace Grid::QCD;

  Grid_init(&argc, &argv);
  GridLogIntegrator.Active(1);
  int threads = GridThread::GetThreads();
  // here make a routine to print all the relevant information on the run
  std::cout << GridLogMessage << "Grid is setup to use " << threads << " threads" << std::endl;

   // Typedefs to simplify notation
  typedef GenericHMCRunner<ImplicitMinimumNorm2> HMCWrapper;  // Uses the default minimum norm
   // Serialiser
//  typedef Grid::JSONReader       Serialiser;
  typedef Grid::XmlReader       Serialiser;


  HMCWrapper TheHMC;
  TheHMC.ReadCommandLine(argc, argv); // these can be parameters from file
  
  // Reader, file should come from command line
  if (TheHMC.ParameterFile.empty()){
    std::cout << "Input file not specified."
              << "Use --ParameterFile option in the command line.\nAborting" 
              << std::endl;
    exit(1);
  }
  Serialiser Reader(TheHMC.ParameterFile);

  RMHMCActionParameters ActionParams(Reader);  

  // Grid from the command line
  TheHMC.Resources.AddFourDimGrid("gauge");
 

  // Checkpointer definition
  CheckpointerParameters CPparams(Reader);
//  TheHMC.Resources.LoadBinaryCheckpointer(CPparams);
  TheHMC.Resources.LoadNerscCheckpointer(CPparams);

  RNGModuleParameters RNGpar(Reader);
  TheHMC.Resources.SetRNGSeeds(RNGpar);

  // Construct observables
  typedef PlaquetteMod<HMCWrapper::ImplPolicy> PlaqObs;
  typedef TopologicalChargeMod<HMCWrapper::ImplPolicy> QObs;
  TheHMC.Resources.AddObservable<PlaqObs>();
  TopologyObsParameters TopParams(Reader);
  TheHMC.Resources.AddObservable<QObs>(TopParams);
  /////////////////////////////////////////////////////////////
  // Collect actions
  WilsonGaugeActionR Waction(ActionParams.gauge_beta);

  
  ActionLevel<HMCWrapper::Field> Level1(1);
  Level1.push_back(&Waction);
  TheHMC.TheAction.push_back(Level1);
  /////////////////////////////////////////////////////////////
  TheHMC.Parameters.initialize(Reader);
  TheHMC.Run(); 

  Grid_finalize();

} // main

/* Examples for input files

JSON

{
    "Checkpointer": {
	"config_prefix": "ckpoint_json_lat",
	"rng_prefix": "ckpoint_json_rng",
	"saveInterval": 10,
	"format": "IEEE64BIG"
    },
    "RandomNumberGenerator": {
	"serial_seeds": "1 2 3 4 6",
	"parallel_seeds": "55 7 8 9 11"
    },
    "Action":{
	"gauge_beta": 5.8
    },
    "TopologyMeasurement":{
	"interval": 1,
	"do_smearing": true,
	"Smearing":{
	    "steps": 200,
	    "step_size": 0.01,
	    "meas_interval": 50,
	    "maxTau": 2.0
	}
    },
    "HMC":{
	"StartTrajectory": 0,
	"Trajectories": 10,
	"MetropolisTest": true,
	"NoMetropolisUntil": 10,
	"StartingType": "HotStart",
	"MD":{
	    "name": "MinimumNorm2",
	    "MDsteps": 40,
	    "trajL": 1.0
	}
    }
}


XML example not provided yet

*/
