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

class HMCRunnerParameters : Serializable {
 public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(HMCRunnerParameters,
                                  double, beta,
                                  double, mass,
                                  int, MaxCGIterations,
                                  double, StoppingCondition,
                                  bool, smearedAction,
                                  int, SaveInterval,
                                  std::string, format,
                                  std::string, conf_prefix,
                                  std::string, rng_prefix,
                                  double, rho,
                                  int, SmearingLevels,
                                  );

  HMCRunnerParameters() {}
};

// Derive from the BinaryHmcRunner (templated for gauge fields)
class HmcRunner : public BinaryHmcRunner {
 public:
  void BuildTheAction(int argc, char **argv)

  {
    int Ndim=5;
    typedef WilsonImplR ImplPolicy;
    typedef WilsonFermionR FermionAction;
    typedef typename FermionAction::FermionField FermionField;

    std::vector<int> simd = GridDefaultSimd(Ndim-1,vComplex::Nsimd());
    simd.push_back(1);


    //UGrid = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), simd, GridDefaultMpi());
    UGrid   =  new GridCartesian(GridDefaultLatt(),simd,GridDefaultMpi());

    // Gauge action
    int Ls = UGrid->_fdimensions[Nd - 1];
    std::vector<RealD> betat(Ls,6.0);
    std::vector<RealD> betas(Ls,5.6);
    betat[Ls-1]= 0.0;
    betas={5.2,5.5,5.8,6,6,5.8,5.5,5.2};
    std:cout << GridLogMessage << "Betas: " << betas << std::endl;
    VariableWilsonGaugeActionR Waction(betas, betat, UGrid);
    //WilsonGaugeActionR Waction(5.6);

    // Collect actions
    ActionLevel<Field> Level1(1);
    Level1.push_back(&Waction);
    TheAction.push_back(Level1);

    // Add observables
    int SaveInterval = 1;
    std::string format = std::string("IEEE64BIG");
    std::string conf_prefix = std::string("ckpoint_lat");
    std::string rng_prefix = std::string("ckpoint_rng");
    BinaryHmcCheckpointer<BinaryHmcRunner::ImplPolicy> Checkpoint(
        conf_prefix, rng_prefix, SaveInterval, format);
    // Can implement also a specific function in the hmcrunner
    // AddCheckpoint (...) that takes the same parameters + a string/tag
    // defining the type of the checkpointer
    // with tags can be implemented by overloading and no ifs
    // Then force all checkpoint to have few common functions
    // return an object that is then passed to the Run function

    PlaquetteLogger<BinaryHmcRunner::ImplPolicy> PlaqLog(
        std::string("Plaquette"));
    ObservablesList.push_back(&PlaqLog);
    ObservablesList.push_back(&Checkpoint);

    Run(argc, argv, Checkpoint);  // no smearing
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

  // Seeds for the random number generators
  std::vector<int> SerSeed({1, 2, 3, 4, 5});
  std::vector<int> ParSeed({6, 7, 8, 9, 10});
  TheHMC.RNGSeeds(SerSeed, ParSeed);

  TheHMC.MDparameters.set(20, 1.0);// MDsteps, traj length

  TheHMC.BuildTheAction(argc, argv);

  Grid_finalize();
}
