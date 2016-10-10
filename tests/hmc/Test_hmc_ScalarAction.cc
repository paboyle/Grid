/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./tests/Test_hmc_WilsonFermionGauge.cc

Copyright (C) 2015

Author: Peter Boyle <pabobyle@ph.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: neo <cossu@post.kek.jp>
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

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

namespace Grid {
namespace QCD {

// Derive from the BinaryHmcRunner (templated for gauge fields)
class HmcRunner : public ScalarBinaryHmcRunner {
 public:
  void BuildTheAction(int argc, char **argv)

  {
    // Notice that the Grid is for reals now
    UGrid = SpaceTimeGrid::makeFourDimGrid(
        GridDefaultLatt(), GridDefaultSimd(Nd, vReal::Nsimd()),
        GridDefaultMpi());
    UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);

    // Real Scalar action
    ScalarActionR Saction(0.11,0.);

    // Collect actions
    ActionLevel<Field, ScalarFields> Level1(1);
    Level1.push_back(&Saction);

    TheAction.push_back(Level1);


    // Add observables and checkpointers
    int SaveInterval = 1;
    std::string format = std::string("IEEE64BIG");
    std::string conf_prefix = std::string("ckpoint_scalar_lat");
    std::string rng_prefix = std::string("ckpoint_scalar_rng");
    BinaryHmcCheckpointer<ScalarBinaryHmcRunner::ImplPolicy> Checkpoint(conf_prefix, rng_prefix, SaveInterval, format);
    ObservablesList.push_back(&Checkpoint);

    Run(argc, argv, Checkpoint);
  };  
};
}
}

int main(int argc, char **argv) {
  Grid_init(&argc, &argv);

  int threads = GridThread::GetThreads();
  std::cout << GridLogMessage << "Grid is setup to use " << threads
            << " threads" << std::endl;

  HmcRunner TheHMC;

  TheHMC.BuildTheAction(argc, argv);
}
