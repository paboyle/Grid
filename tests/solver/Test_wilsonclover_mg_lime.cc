/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/solver/Test_wilsonclover_mg_mp.cc

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
/* */


#include <Grid/Grid.h>
#include <Test_multigrid_common.h>

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

int main(int argc, char **argv) {


  Grid_init(&argc, &argv);

  // clang-format off
  GridCartesian         *FGrid_d   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd, vComplexD::Nsimd()), GridDefaultMpi());
  GridCartesian         *FGrid_f   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd, vComplexF::Nsimd()), GridDefaultMpi());
  GridRedBlackCartesian *FrbGrid_d = SpaceTimeGrid::makeFourDimRedBlackGrid(FGrid_d);
  GridRedBlackCartesian *FrbGrid_f = SpaceTimeGrid::makeFourDimRedBlackGrid(FGrid_f);
  // clang-format on

  std::vector<int> fSeeds({1, 2, 3, 4});
  GridParallelRNG  fPRNG(FGrid_d);
  fPRNG.SeedFixedIntegers(fSeeds);

  // clang-format off
  LatticeFermionD       src_d(FGrid_d); gaussian(fPRNG, src_d);
  LatticeFermionD resultMGD_d(FGrid_d); resultMGD_d = zero;
  LatticeFermionD resultMGF_d(FGrid_d); resultMGF_d = zero;
  LatticeGaugeFieldD    Umu_d(FGrid_d);

#if 0
  {
    FieldMetaData header;
    std::string file("./qcdsf.769.00399.lime");
    std::cout <<GridLogMessage<<"**************************************"<<std::endl;
    std::cout <<GridLogMessage<<"** Reading back ILDG conf    *********"<<std::endl;
    std::cout <<GridLogMessage<<"**************************************"<<std::endl;
    IldgReader _IldgReader;
    _IldgReader.open(file);
    _IldgReader.readConfiguration(Umu_d,header);
    _IldgReader.close();

  }
#else
{
  FieldMetaData header;
  std::string file("./ckpoint_lat.IEEE64BIG.1100");
  NerscIO::readConfiguration(Umu_d,header,file);
}
#endif
  // SU3::HotConfiguration(fPRNG, Umu_d);

  LatticeGaugeFieldF    Umu_f(FGrid_f); precisionChange(Umu_f, Umu_d);
  // clang-format on

  RealD mass  = -0.25;
  RealD csw_r = 1.0;
  RealD csw_t = 1.0;

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

  LevelInfo levelInfo_d(FGrid_d, mgParams);
  LevelInfo levelInfo_f(FGrid_f, mgParams);

  // Note: We do chiral doubling, so actually only nbasis/2 full basis vectors are used
  const int nbasis = 40;

  WilsonCloverFermionD Dwc_d(Umu_d, *FGrid_d, *FrbGrid_d, mass, csw_r, csw_t);
  WilsonCloverFermionF Dwc_f(Umu_f, *FGrid_f, *FrbGrid_f, mass, csw_r, csw_t);

  MdagMLinearOperator<WilsonCloverFermionD, LatticeFermionD> MdagMOpDwc_d(Dwc_d);
  MdagMLinearOperator<WilsonCloverFermionF, LatticeFermionF> MdagMOpDwc_f(Dwc_f);

  std::cout << GridLogMessage << "**************************************************" << std::endl;
  std::cout << GridLogMessage << "Testing single-precision Multigrid for Wilson Clover" << std::endl;
  std::cout << GridLogMessage << "**************************************************" << std::endl;

  auto MGPreconDwc_f = createMGInstance<vSpinColourVectorF, vTComplexF, nbasis, WilsonCloverFermionF>(mgParams, levelInfo_f, Dwc_f, Dwc_f);

  MGPreconDwc_f->setup();

  if(GridCmdOptionExists(argv, argv + argc, "--runchecks")) {
    MGPreconDwc_f->runChecks(1e-6);
  }

  MixedPrecisionFlexibleGeneralisedMinimalResidual<LatticeFermionD, LatticeFermionF> MPFGMRESPREC(
    1.0e-12, 50000, FGrid_f, *MGPreconDwc_f, 100, false);

  std::cout << std::endl << "Starting with a new solver" << std::endl;
  MPFGMRESPREC(MdagMOpDwc_d, src_d, resultMGF_d);

  MGPreconDwc_f->reportTimings();

  if(GridCmdOptionExists(argv, argv + argc, "--docomparison")) {

    std::cout << GridLogMessage << "**************************************************" << std::endl;
    std::cout << GridLogMessage << "Testing double-precision Multigrid for Wilson Clover" << std::endl;
    std::cout << GridLogMessage << "**************************************************" << std::endl;

    auto MGPreconDwc_d = createMGInstance<vSpinColourVectorD, vTComplexD, nbasis, WilsonCloverFermionD>(mgParams, levelInfo_d, Dwc_d, Dwc_d);

    MGPreconDwc_d->setup();

    if(GridCmdOptionExists(argv, argv + argc, "--runchecks")) {
      MGPreconDwc_d->runChecks(1e-13);
    }

    FlexibleGeneralisedMinimalResidual<LatticeFermionD> FGMRESPREC(1.0e-12, 50000, *MGPreconDwc_d, 100, false);

    std::cout << std::endl << "Starting with a new solver" << std::endl;
    FGMRESPREC(MdagMOpDwc_d, src_d, resultMGD_d);

    MGPreconDwc_d->reportTimings();

    std::cout << GridLogMessage << "**************************************************" << std::endl;
    std::cout << GridLogMessage << "Comparing single-precision Multigrid with double-precision one for Wilson Clover" << std::endl;
    std::cout << GridLogMessage << "**************************************************" << std::endl;

    LatticeFermionD diffFullSolver(FGrid_d);

    RealD deviationFullSolver = axpy_norm(diffFullSolver, -1.0, resultMGF_d, resultMGD_d);

    // clang-format off
    LatticeFermionF src_f(FGrid_f);    precisionChange(src_f, src_d);
    LatticeFermionF resMGF_f(FGrid_f); resMGF_f = zero;
    LatticeFermionD resMGD_d(FGrid_d); resMGD_d = zero;
    // clang-format on

    (*MGPreconDwc_f)(src_f, resMGF_f);
    (*MGPreconDwc_d)(src_d, resMGD_d);

    LatticeFermionD diffOnlyMG(FGrid_d);
    LatticeFermionD resMGF_d(FGrid_d);
    precisionChange(resMGF_d, resMGF_f);

    RealD deviationOnlyPrec = axpy_norm(diffOnlyMG, -1.0, resMGF_d, resMGD_d);

    // clang-format off
    std::cout << GridLogMessage << "Absolute difference between FGMRES preconditioned by double and single precicision MG: " << deviationFullSolver                      << std::endl;
    std::cout << GridLogMessage << "Relative deviation  between FGMRES preconditioned by double and single precicision MG: " << deviationFullSolver / norm2(resultMGD_d) << std::endl;
    std::cout << GridLogMessage << "Absolute difference between one iteration of MG Prec in double and single precision:   " << deviationOnlyPrec                        << std::endl;
    std::cout << GridLogMessage << "Relative deviation  between one iteration of MG Prec in double and single precision:   " << deviationOnlyPrec / norm2(resMGD_d)      << std::endl;
    // clang-format on
  }

  Grid_finalize();
}
