/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./examples/Example_gamma5_block_lanczos.cc

    Copyright (C) 2026

    Author: Chulwoo Jung <chulwoo@bnl.gov>

    γ5-Block Lanczos example for the Wilson Dirac operator.

    Reads a gauge configuration from "config" (NERSC format) and Lanczos
    parameters from "LanParams.xml".  Runs Gamma5BlockLanczos to
    compute eigenvalues of D_W directly (not H_W = γ5 D_W).

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

#include <cstdlib>

#include <Grid/Grid.h>
#include <Grid/lattice/PaddedCell.h>
#include <Grid/stencil/GeneralLocalStencil.h>

using namespace std;
using namespace Grid;

namespace Grid {

struct LanczosParameters : Serializable {
  GRID_SERIALIZABLE_CLASS_MEMBERS(LanczosParameters,
    RealD,   mass,
    Integer, Nstop,
    Integer, Nk,
    Integer, Np,
    Integer, maxIter,
    Integer, reorthog,
    Integer, verify,
    Integer, ReadEvec,
    RealD,   resid)

  LanczosParameters()
    : mass(-0.5), Nstop(10), Nk(20), maxIter(100),
      reorthog(1), verify(0), ReadEvec(0), resid(1e-8)
  {}

  template<class ReaderClass>
  LanczosParameters(Reader<ReaderClass>& r) { initialize(r); }

  template<class ReaderClass>
  void initialize(Reader<ReaderClass>& r) {
    read(r, "LanczosParameters", *this);
  }
};

} // namespace Grid

template<class T>
void writeField(T& in, std::string const& fname) {
#ifdef HAVE_LIME
  std::cout << GridLogMessage << "Writing to: " << fname << std::endl;
  Grid::emptyUserRecord record;
  Grid::ScidacWriter WR(in.Grid()->IsBoss());
  WR.open(fname);
  WR.writeScidacFieldRecord(in, record, 0);
  WR.close();
#endif
}

typedef WilsonFermionD WilsonOp;
typedef typename WilsonFermionD::FermionField FermionField;

int main(int argc, char** argv)
{
  Grid_init(&argc, &argv);

  GridCartesian* UGrid = SpaceTimeGrid::makeFourDimGrid(
      GridDefaultLatt(),
      GridDefaultSimd(Nd, vComplex::Nsimd()),
      GridDefaultMpi());
  GridRedBlackCartesian* UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);

  std::vector<int> seeds4({1, 2, 3, 4});
  GridParallelRNG RNG4(UGrid);
  RNG4.SeedFixedIntegers(seeds4);

  // Read gauge configuration
  LatticeGaugeField Umu(UGrid);
  FieldMetaData header;
  std::string configFile("config");
  NerscIO::readConfiguration(Umu, header, configFile);
  std::cout << GridLogMessage << "Loaded gauge configuration: " << configFile << std::endl;

  // Read Lanczos parameters
  LanczosParameters Params;
  {
    XmlReader rd("LanParams.xml");
    read(rd, "LanczosParameters", Params);
  }
  std::cout << GridLogMessage << Params << std::endl;

  // Build Wilson Dirac operator and wrap in a non-Hermitian linear operator
  std::vector<Complex> boundary = {1, 1, 1, -1};
  WilsonOp::ImplParams WilsonParams(boundary);
  WilsonOp Dwilson(Umu, *UGrid, *UrbGrid, Params.mass, WilsonParams);
  NonHermitianLinearOperator<WilsonOp, FermionField> DLinOp(Dwilson);

  // γ5 functor: for 4D Wilson fermions γ5 is Gamma(Gamma5)
  Gamma G5(Gamma::Algebra::Gamma5);
  auto gamma5 = [&G5](const FermionField& in, FermionField& out) {
    out = G5 * in;
  };

  // Starting vectors: two independent random vectors
  FermionField src(UGrid), src2(UGrid);
  random(RNG4, src);
  random(RNG4, src2);
  std::cout << GridLogMessage << "Using two random starting vectors" << std::endl;

  std::cout << GridLogMessage << std::endl;
  std::cout << GridLogMessage << "*******************************************" << std::endl;
  std::cout << GridLogMessage << " Running γ5-Block Lanczos" << std::endl;
  std::cout << GridLogMessage << "   mass     = " << Params.mass     << std::endl;
  std::cout << GridLogMessage << "   Nk       = " << Params.Nk       << std::endl;
  std::cout << GridLogMessage << "   maxIter  = " << Params.maxIter  << std::endl;
  std::cout << GridLogMessage << "   Nstop    = " << Params.Nstop    << std::endl;
  std::cout << GridLogMessage << "   resid    = " << Params.resid    << std::endl;
  std::cout << GridLogMessage << "   reorthog = " << Params.reorthog << std::endl;
  std::cout << GridLogMessage << "   verify   = " << Params.verify   << std::endl;
  std::cout << GridLogMessage << "*******************************************" << std::endl;
  std::cout << GridLogMessage << std::endl;

  Gamma5BlockLanczos<FermionField> G5BL(DLinOp, UGrid, gamma5, Params.resid);
  G5BL.doEvalCheck = (Params.verify != 0);
  G5BL.doVerify = (Params.verify != 0);
//  G5BL(src, src2, Params.maxIter, Params.Nstop, Params.reorthog != 0);
  G5BL.restart(src, src2, Params.maxIter, Params.Nk+Params.Np, Params.Nk, Params.Nstop, Params.reorthog != 0, EvalNormSmall);
//  G5BL.implicitRestart(src, src2, Params.maxIter, Params.Nk+Params.Np, Params.Nk, Params.Nstop, Params.reorthog != 0, EvalNormSmall);
  if (Params.verify ) G5BL.verify("after restart");

  // Summary of results
  Eigen::VectorXcd        evals     = G5BL.getEvals();
  std::vector<RealD>      residuals = G5BL.getResiduals();
  std::vector<FermionField> evecs   = G5BL.getEvecs();

  int Nout = (int)evals.size();
  std::cout << GridLogMessage << std::endl;
  std::cout << GridLogMessage << "*******************************************" << std::endl;
  std::cout << GridLogMessage << " γ5-Block Lanczos: " << Nout << " Ritz pairs" << std::endl;
  std::cout << GridLogMessage << "*******************************************" << std::endl;
  for (int i = 0; i < Nout; i++) {
    std::cout << GridLogMessage
              << "  [" << std::setw(3) << i << "]"
              << "  lambda = " << evals(i)
              << "  |res| = " << residuals[i] << std::endl;
  }

  // Write the first Nstop eigenvectors (in SCIDAC format when LIME is available)
  int Nwrite = std::min((int)Params.Nstop, Nout);
  for (int i = 0; i < Nwrite; i++) {
    std::string fname = "./g5bl_evec_m" + std::to_string(Params.mass)
                        + "_" + std::to_string(i);
    writeField(evecs[i], fname);
  }

  std::cout << GridLogMessage << std::endl;
  std::cout << GridLogMessage << "Done" << std::endl;

  Grid_finalize();
  return 0;
}
