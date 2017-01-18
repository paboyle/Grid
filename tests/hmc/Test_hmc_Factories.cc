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
  // ifdefs ?? Local makefile suggestion , policy as make parameter
typedef QCD::PeriodicGimplR ImplementationPolicy;
typedef QCD::NoHirep RepresentationPolicy;

static Registrar< HMCLeapFrog<ImplementationPolicy, RepresentationPolicy, XmlReader>      , HMCRunnerModuleFactory<hmc_string, XmlReader> > __HMCLFmodXMLInit("LeapFrog");
static Registrar< HMCMinimumNorm2<ImplementationPolicy, RepresentationPolicy, XmlReader>  , HMCRunnerModuleFactory<hmc_string, XmlReader> > __HMCMN2modXMLInit("MinimumNorm2");
static Registrar< HMCForceGradient<ImplementationPolicy, RepresentationPolicy, XmlReader> , HMCRunnerModuleFactory<hmc_string, XmlReader> > __HMCFGmodXMLInit("ForceGradient");
}

int main(int argc, char **argv) {
  using namespace Grid;
  using namespace Grid::QCD;

  Grid_init(&argc, &argv);
  int threads = GridThread::GetThreads();
  // here make a routine to print all the relevant information on the run
  std::cout << GridLogMessage << "Grid is setup to use " << threads << " threads" << std::endl;

   // Typedefs to simplify notation 
  typedef Grid::XmlReader InputFileReader; 

  // Reader, file should come from command line
  InputFileReader Reader("input.wilson_gauge.params.xml");

  // Test HMC factory (put in an external file)
  auto &HMCfactory = HMCRunnerModuleFactory<hmc_string, InputFileReader >::getInstance();
  // Simplify this step (IntergratorName field?)
  HMCparameters HMCpar(Reader);
  
  // Construct the module
  auto myHMCmodule = HMCfactory.create(HMCpar.MD.name, Reader);

  myHMCmodule->getPtr()->initialize(Reader);
  myHMCmodule->getPtr()->Run();


  Grid_finalize();

} // main
