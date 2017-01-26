/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./tests/Test_hmc_Factories.cc

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

See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
/*  END LEGAL */
#include <Grid/Grid.h>

namespace Grid{

// Put this section in a separate header
// ifdefs ?? Local makefile suggestion , policy as make parameter
typedef QCD::PeriodicGimplR   ImplementationPolicy;
typedef QCD::WilsonImplR      FermionImplementationPolicy;
typedef QCD::NoHirep          RepresentationPolicy;
typedef Grid::XmlReader       Serialiser;

// Register all object names
#include "Grid/qcd/modules/Registration.h"

} // Grid

int main(int argc, char **argv) {
  using namespace Grid;
  using namespace Grid::QCD;

  Grid_init(&argc, &argv);

  int threads = GridThread::GetThreads();
  // here make a routine to print all the relevant information on the run
  std::cout << GridLogMessage << "Grid is setup to use " << threads << " threads" << std::endl;

   // Typedefs to simplify notation 
  //typedef XmlReader InputFileReader; 

  // Reader, file should come from command line
  Serialiser Reader("input.wilson_gauge.params.xml");

  // Test HMC factory (put in an external file)
  auto &HMCfactory = HMCModuleFactory::getInstance();
  // Simplify this step (IntergratorName field?)
  HMCparameters HMCpar(Reader);
  
  // Construct the module
  auto HMCmodule = HMCfactory.create(HMCpar.MD.name, Reader);

  HMCmodule->getPtr()->initialize(Reader);
  HMCmodule->getPtr()->Run();

  Grid_finalize();
  return 0;

} // main

