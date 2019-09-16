/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/solver/Test_wilson_mg.cc

    Copyright (C) 2015-2018

    Author: Daniel Richtmann <daniel.richtmann@ur.de>

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
#include <Test_multigrid_common.h>

using namespace std;
using namespace Grid;

int main(int argc, char **argv) {

  Grid_init(&argc, &argv);

  GridCartesian *        FGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd, vComplex::Nsimd()), GridDefaultMpi());
  GridRedBlackCartesian *FrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(FGrid);

  std::vector<int> fSeeds({1, 2, 3, 4});
  GridParallelRNG  fPRNG(FGrid);
  fPRNG.SeedFixedIntegers(fSeeds);

  // clang-format off
  LatticeFermion    src(FGrid); gaussian(fPRNG, src);
  LatticeFermion result(FGrid); result = Zero();
  LatticeGaugeField Umu(FGrid); SU3::HotConfiguration(fPRNG, Umu);
  // clang-format on

  RealD mass = -0.25;

  MultiGridParams mgParams;
  std::string     inputXml{"./mg_params.xml"};

  if(GridCmdOptionExists(argv, argv + argc, "--inputxml")) {
    inputXml = GridCmdOptionPayload(argv, argv + argc, "--inputxml");
    assert(inputXml.length() != 0);
  }

  {
    XmlWriter writer("mg_params_template.xml");
    write(writer, "Params", mgParams);
    std::cout << GridLogMessage << "Written mg_params_template.xml" << std::endl;

    XmlReader reader(inputXml);
    read(reader, "Params", mgParams);
    std::cout << GridLogMessage << "Read in " << inputXml << std::endl;
  }

  checkParameterValidity(mgParams);
  std::cout << mgParams << std::endl;

  LevelInfo levelInfo(FGrid, mgParams);

  // Note: We do chiral doubling, so actually only nbasis/2 full basis vectors are used
  const int nbasis = 40;

  WilsonFermionR Dw(Umu, *FGrid, *FrbGrid, mass);

  MdagMLinearOperator<WilsonFermionR, LatticeFermion> MdagMOpDw(Dw);

  std::cout << GridLogMessage << "**************************************************" << std::endl;
  std::cout << GridLogMessage << "Testing Multigrid for Wilson" << std::endl;
  std::cout << GridLogMessage << "**************************************************" << std::endl;

  TrivialPrecon<LatticeFermion> TrivialPrecon;
  auto MGPreconDw = createMGInstance<vSpinColourVector, vTComplex, nbasis, WilsonFermionR>(mgParams, levelInfo, Dw, Dw);

  MGPreconDw->setup();

  if(GridCmdOptionExists(argv, argv + argc, "--runchecks")) {
    RealD toleranceForMGChecks = (getPrecision<LatticeFermion>::value == 1) ? 1e-6 : 1e-13;
    MGPreconDw->runChecks(toleranceForMGChecks);
  }

  std::vector<std::unique_ptr<OperatorFunction<LatticeFermion>>> solversDw;

  solversDw.emplace_back(new ConjugateGradient<LatticeFermion>(1.0e-12, 50000, false));
  solversDw.emplace_back(new FlexibleGeneralisedMinimalResidual<LatticeFermion>(1.0e-12, 50000, TrivialPrecon, 100, false));
  solversDw.emplace_back(new FlexibleGeneralisedMinimalResidual<LatticeFermion>(1.0e-12, 50000, *MGPreconDw, 100, false));

  for(auto const &solver : solversDw) {
    std::cout << std::endl << "Starting with a new solver" << std::endl;
    result = zero;
    (*solver)(MdagMOpDw, src, result);
  }

  MGPreconDw->reportTimings();

  Grid_finalize();
}
