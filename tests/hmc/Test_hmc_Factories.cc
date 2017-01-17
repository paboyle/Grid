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

int main(int argc, char **argv) {
  using namespace Grid;
  using namespace Grid::QCD;

  Grid_init(&argc, &argv);
  int threads = GridThread::GetThreads();
  // here make a routine to print all the relevant information on the run
  std::cout << GridLogMessage << "Grid is setup to use " << threads << " threads" << std::endl;

   // Typedefs to simplify notation 
  typedef GenericHMCRunner<MinimumNorm2> HMCWrapper;  // Uses the default minimum norm
  typedef Grid::XmlReader InputFileReader; 

  // Reader, file should come from command line
  InputFileReader Reader("input.wilson_gauge.params.xml");

  HMCWrapper TheHMC;

  TheHMC.Parameters.initialize(Reader);
  TheHMC.Resources.initialize(Reader);

  // Construct observables
  // here there is too much indirection 
  PlaquetteLogger<HMCWrapper::ImplPolicy> PlaqLog("Plaquette");
  TheHMC.ObservablesList.push_back(&PlaqLog);
  //////////////////////////////////////////////

  /////////////////////////////////////////////////////////////
  // Collect actions, here use more encapsulation
  // need wrappers of the fermionic classes 
  // that have a complex construction

  // standard
  RealD beta = 5.6 ;
  WilsonGaugeActionR Waction(beta);

  ActionLevel<HMCWrapper::Field> Level1(1);
  Level1.push_back(&Waction);
  TheHMC.TheAction.push_back(Level1);
  /////////////////////////////////////////////////////////////

  // eventually smearing here
  // ...
  ////////////////////////////////////////////////////////////////

  TheHMC.Run();  // no smearing

  Grid_finalize();

} // main
