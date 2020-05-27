/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./benchmarks/Benchmark_wilson_clover_mrhs.cc

    Copyright (C) 2015 - 2020

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
#include <Benchmark_helpers.h>

using namespace Grid;

struct WorkPerSite {
  double flop;
  double elem;
  double word = sizeof(Complex); // good default with possibility to override

  double byte() const { return elem * word; }
  double intensity() const { return elem == 0. ? 0. : flop/byte(); }
};

inline double cMatMulFlop(int n) { return ((8 - 2 / (double)n) * n * n); }

#define doComparison(name4d, name5d, vecs4d, vecs5d, tolerance) \
  do { \
    GridBase* grid = vecs4d[0].Grid(); \
    LatticeFermion vec_tmp(grid); vec_tmp.Checkerboard() = vecs4d[0].Checkerboard(); \
    LatticeFermion diff(grid); \
    \
    for(int i = 0; i < nrhs; i++) { \
      ExtractSlice(vec_tmp, vecs5d, i, 0); \
      diff = vecs4d[i] - vec_tmp; \
      auto diff2 = norm2(diff); \
      std::cout << GridLogMessage << "vector " << i << ": norm2(" << #name4d << "* v) = " << norm2(vecs4d[i]) \
                << " norm2(" << #name5d << " * v) = " << norm2(vec_tmp) << " diff = " << diff2 << std::endl; \
      assert(diff2 == 0.); \
    } \
  } while(0)

int main(int argc, char** argv) {
  Grid_init(&argc, &argv);

  /////////////////////////////////////////////////////////////////////////////
  //                          Read from command line                         //
  /////////////////////////////////////////////////////////////////////////////

  // clang-format off
  int      nrhs           = Commandline::readInt(&argc, &argv, "--nrhs", 20);
  uint64_t nIterMin       = Commandline::readInt(&argc, &argv, "--miniter", 1000);
  uint64_t nSecMin        = Commandline::readInt(&argc, &argv, "--minsec", 5);
  bool     countPerformed = Commandline::readToggle(&argc, &argv, "--performed"); // calculate performance with actually performed traffic rather than minimum required
  // clang-format on

  /////////////////////////////////////////////////////////////////////////////
  //                              General setup                              //
  /////////////////////////////////////////////////////////////////////////////

  // clang-format off
  GridCartesian*         UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd, vComplex::Nsimd()), GridDefaultMpi());
  GridRedBlackCartesian* UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian*         FGrid   = SpaceTimeGrid::makeFiveDimGrid(nrhs, UGrid);
  GridRedBlackCartesian* FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(nrhs, UGrid);
  // clang-format on

  UGrid->show_decomposition();
  UrbGrid->show_decomposition();
  FGrid->show_decomposition();
  FrbGrid->show_decomposition();

  GridParallelRNG UPRNG(UGrid);
  GridParallelRNG FPRNG(FGrid);

  std::vector<int> seeds({1, 2, 3, 4});

  UPRNG.SeedFixedIntegers(seeds);
  FPRNG.SeedFixedIntegers(seeds);

  RealD tol = getPrecision<vComplex>::value == 2 ? 1e-15 : 1e-7;

  /////////////////////////////////////////////////////////////////////////////
  //                    Setup of Dirac Matrix and Operator                   //
  /////////////////////////////////////////////////////////////////////////////

  LatticeGaugeField Umu(UGrid); SU3::HotConfiguration(UPRNG, Umu);

  RealD mass = 0.5;
  RealD csw  = 1.0;

  typename WilsonCloverFermionR::ImplParams implParams;
  WilsonAnisotropyCoefficients              anisParams;

  std::vector<Complex> boundary_phases(Nd, 1.);
  if(Commandline::readToggle(&argc, &argv, "--antiperiodic")) boundary_phases[Nd - 1] = -1.;
  implParams.boundary_phases = boundary_phases;

  WilsonFermionR       Dw(Umu, *UGrid, *UrbGrid, mass, implParams, anisParams);
  WilsonCloverFermionR Dwc(Umu, *UGrid, *UrbGrid, mass, csw, csw, anisParams, implParams);

  WilsonMRHSFermionR       Dw5(Umu, *FGrid, *FrbGrid, *UGrid, *UrbGrid, mass, implParams);
  WilsonCloverMRHSFermionR Dwc5(Umu, *FGrid, *FrbGrid, *UGrid, *UrbGrid, mass, csw, csw, anisParams, implParams);

  /////////////////////////////////////////////////////////////////////////////
  //             Calculate numbers needed for performance figures            //
  /////////////////////////////////////////////////////////////////////////////

  double volumeUGrid   = std::accumulate(UGrid->_fdimensions.begin(),   UGrid->_fdimensions.end(),   1, std::multiplies<double>());
  double volumeUrbGrid = std::accumulate(UrbGrid->_fdimensions.begin(), UrbGrid->_fdimensions.end(), 1, std::multiplies<double>());
  double volumeFGrid   = std::accumulate(FGrid->_fdimensions.begin(),   FGrid->_fdimensions.end(),   1, std::multiplies<double>());
  double volumeFrbGrid = std::accumulate(FrbGrid->_fdimensions.begin(), FrbGrid->_fdimensions.end(), 1, std::multiplies<double>());

  WorkPerSite wilson_dhop, clover_dhop;
  WorkPerSite wilson_diag, clover_diag;
  WorkPerSite wilson_diag_necessary, clover_diag_necessary;
  WorkPerSite wilson_diag_performed, clover_diag_performed;
  WorkPerSite wilson_dir, clover_dir;
  WorkPerSite wilson_full, clover_full;
  WorkPerSite spinor, gaugemat, clovmat, clovmat_packed;

  spinor.elem         = Ns * Nc;
  gaugemat.elem       = Nc * Nc;
  clovmat.elem        = spinor.elem * spinor.elem; // full spincolor matrix
  clovmat_packed.elem = Nhs * (Nc + (Nhs*Nc) * (Nhs*Nc - 1) / 2); // 2 * (6 diag reals (= 3 complex) + 15 complex upper triang)

  wilson_dhop.flop = 8 * spinor.elem + 8*Nhs*cMatMulFlop(Nc) + 7 * (2*spinor.elem); // spin-proj + su3 x 2 half spinors + accum output
  wilson_dhop.elem = 8 * (gaugemat.elem + spinor.elem) + spinor.elem;

  clover_dhop.flop = wilson_dhop.flop;
  clover_dhop.elem = wilson_dhop.elem;

  wilson_diag_necessary.flop = 6 * spinor.elem;
  wilson_diag_performed.flop = 6 * spinor.elem;
  wilson_diag_necessary.elem = 2 * spinor.elem;
  wilson_diag_performed.elem = 2 * spinor.elem;

  clover_diag_necessary.flop = 2 * (cMatMulFlop(Nhs*Nc) - 4*Nhs*Nc); // 2 6x6 matmuls - diff between compl mul and real mul on diagonal
  clover_diag_performed.flop = cMatMulFlop(spinor.elem);             // clovmat x spinor

  clover_diag_necessary.elem = clovmat_packed.elem + 2 * spinor.elem;
  clover_diag_performed.elem = clovmat.elem + 2 * spinor.elem;

  wilson_dir.flop = spinor.elem + Nhs*cMatMulFlop(Nc); // spin project + su3 x half spinor + TODO reconstruct?
  wilson_dir.elem = gaugemat.elem + 2 * spinor.elem; // gauge mat + neigh spinor + output spinor

  clover_dir.flop = wilson_dir.flop;
  clover_dir.elem = wilson_dir.elem;

  if(countPerformed) { // count traffic that is actually performed by the hardware
    std::cout << GridLogMessage << "Using traffic value that is actually transfered for the full matrices" << std::endl;
    wilson_diag = wilson_diag_performed;
    clover_diag = clover_diag_performed;

    wilson_full.flop = wilson_dhop.flop + 8 * spinor.elem; // axpy = 8 flop
    wilson_full.elem = wilson_dhop.elem + 3 * spinor.elem; // axpy = 3 spinors

    clover_full.flop = clover_dhop.flop + clover_diag.flop + 2 * spinor.elem; // += = 2 flop
    clover_full.elem = clover_dhop.elem + clover_diag.elem + 3 * spinor.elem; // += = 3 spinors
  } else { // count minimal traffic that is actually necessary (e.g., don't count the extra traffic for the spinors since it wouldn't be needed in a fused impl for M)
    std::cout << GridLogMessage << "Using minimal traffic value that is actually necessary for the full matrices" << std::endl;
    wilson_diag = wilson_diag_necessary;
    clover_diag = clover_diag_necessary;

    wilson_full.flop = wilson_dhop.flop + wilson_diag.flop;
    wilson_full.elem = wilson_dhop.elem + wilson_diag.elem - 2*spinor.elem; // only one spinor load + store required, not 2

    clover_full.flop = clover_dhop.flop + clover_diag.flop;
    clover_full.elem = clover_dhop.elem + clover_diag.elem - 2*spinor.elem; // only one spinor load + store required, not 2
  }

  // NOTE flop values per site from output of quda's dslash_test binary:
  // wilson = 1320
  // mobius = 1320
  // clover = 1824
  // twisted-mass = 1392
  // twisted-clover = 1872
  // domain-wall-4d = 1320
  // domain-wall = 1419
  // clover-hasenbusch-twist = 1824

  // clang-format off
  std::cout << GridLogPerformance << "4d volume: " << volumeUGrid << std::endl;
  std::cout << GridLogPerformance << "5d volume: " << volumeFGrid << std::endl;
  std::cout << GridLogPerformance << "Dw.Dhop flop/site, byte/site, flop/byte: "   << wilson_dhop.flop << " " << wilson_dhop.byte() << " " << wilson_dhop.intensity() << std::endl;
  std::cout << GridLogPerformance << "Dwc.Dhop flop/site, byte/site, flop/byte: "  << clover_dhop.flop << " " << clover_dhop.byte() << " " << clover_dhop.intensity() << std::endl;
  std::cout << GridLogPerformance << "Dw.Mdiag flop/site, byte/site, flop/byte: "  << wilson_diag.flop << " " << wilson_diag.byte() << " " << wilson_diag.intensity() << std::endl;
  std::cout << GridLogPerformance << "Dwc.Mdiag flop/site, byte/site, flop/byte: " << clover_diag.flop << " " << clover_diag.byte() << " " << clover_diag.intensity() << std::endl;
  std::cout << GridLogPerformance << "Dw.Mdir flop/site, byte/site, flop/byte: "   << wilson_dir.flop  << " " << wilson_dir.byte()  << " " << wilson_dir.intensity()  << std::endl;
  std::cout << GridLogPerformance << "Dwc.Mdir flop/site, byte/site, flop/byte: "  << clover_dir.flop  << " " << clover_dir.byte()  << " " << clover_dir.intensity()  << std::endl;
  std::cout << GridLogPerformance << "Dw.M flop/site, byte/site, flop/byte: "      << wilson_full.flop << " " << wilson_full.byte() << " " << wilson_full.intensity() << std::endl;
  std::cout << GridLogPerformance << "Dwc.M flop/site, byte/site, flop/byte: "     << clover_full.flop << " " << clover_full.byte() << " " << clover_full.intensity() << std::endl;
  // clang-format on

  /////////////////////////////////////////////////////////////////////////////
  //                             Setup of vectors                            //
  /////////////////////////////////////////////////////////////////////////////

  std::vector<LatticeFermion> vecs_src_4d(nrhs, UGrid);
  std::vector<LatticeFermion> vecs_res_4d(nrhs, UGrid);

  LatticeFermion vecs_src_5d(FGrid);
  LatticeFermion vecs_res_5d(FGrid);

  for(int rhs=0; rhs<nrhs; rhs++) {
    random(UPRNG, vecs_src_4d[rhs]);
    InsertSlice(vecs_src_4d[rhs], vecs_src_5d, rhs, 0);
  }

  /////////////////////////////////////////////////////////////////////////////
  //                      Benchmarks for Wilson fermions                     //
  /////////////////////////////////////////////////////////////////////////////

  // Wilson hopping term = "dslash" = "Dhop" //////////////////////////////////

  {
    for(int rhs=0; rhs<nrhs; rhs++) vecs_res_4d[rhs] = Zero();
    vecs_res_5d = Zero();

    double flop = volumeUGrid * nrhs * wilson_dhop.flop;
    double byte = volumeUGrid * nrhs * wilson_dhop.byte();

    BenchmarkFunctionMRHS(Dw.Dhop,
                          flop, byte, nIterMin, nSecMin, nrhs,
                          vecs_src_4d[rhs], vecs_res_4d[rhs], 0);

    BenchmarkFunction(Dw5.Dhop,
                      flop, byte, nIterMin, nSecMin,
                      vecs_src_5d, vecs_res_5d, 0);

    doComparison(Dw.Dhop, Dw5.Dhop, vecs_res_4d, vecs_res_5d, tol);
  }

  // Wilson half cb hopping term = "Meooe" ////////////////////////////////////

  {
    for(int rhs=0; rhs<nrhs; rhs++) vecs_res_4d[rhs] = Zero();
    vecs_res_5d = Zero();

    std::vector<LatticeFermion> vecs_src_4d_e(nrhs, UrbGrid);
    std::vector<LatticeFermion> vecs_src_4d_o(nrhs, UrbGrid);
    std::vector<LatticeFermion> vecs_res_4d_e(nrhs, UrbGrid);
    std::vector<LatticeFermion> vecs_res_4d_o(nrhs, UrbGrid);

    for(int rhs=0; rhs<nrhs; rhs++) {
      pickCheckerboard(Even, vecs_src_4d_e[rhs], vecs_src_4d[rhs]);
      pickCheckerboard(Odd,  vecs_src_4d_o[rhs], vecs_src_4d[rhs]);
      pickCheckerboard(Even, vecs_res_4d_e[rhs], vecs_res_4d[rhs]);
      pickCheckerboard(Odd,  vecs_res_4d_o[rhs], vecs_res_4d[rhs]);
    }

    // clang-format off
    LatticeFermion vecs_src_5d_e(FrbGrid); pickCheckerboard(Even, vecs_src_5d_e, vecs_src_5d);
    LatticeFermion vecs_src_5d_o(FrbGrid); pickCheckerboard(Odd,  vecs_src_5d_o, vecs_src_5d);
    LatticeFermion vecs_res_5d_e(FrbGrid); pickCheckerboard(Even, vecs_res_5d_e, vecs_res_5d);
    LatticeFermion vecs_res_5d_o(FrbGrid); pickCheckerboard(Odd,  vecs_res_5d_o, vecs_res_5d);
    // clang-format on

    double flop = volumeUrbGrid * nrhs * wilson_dhop.flop;
    double byte = volumeUrbGrid * nrhs * wilson_dhop.byte();

    BenchmarkFunctionMRHS(Dw.Meooe,
                          flop, byte, nIterMin, nSecMin, nrhs,
                          vecs_src_4d_e[rhs], vecs_res_4d_o[rhs]);
    BenchmarkFunctionMRHS(Dw.Meooe,
                          flop, byte, nIterMin, nSecMin, nrhs,
                          vecs_src_4d_o[rhs], vecs_res_4d_e[rhs]);

    BenchmarkFunction(Dw5.Meooe,
                      flop, byte, nIterMin, nSecMin,
                      vecs_src_5d_e, vecs_res_5d_o);
    BenchmarkFunction(Dw5.Meooe,
                      flop, byte, nIterMin, nSecMin,
                      vecs_src_5d_o, vecs_res_5d_e);

    for(int rhs=0; rhs<nrhs; rhs++) {
      setCheckerboard(vecs_res_4d[rhs], vecs_res_4d_e[rhs]);
      setCheckerboard(vecs_res_4d[rhs], vecs_res_4d_o[rhs]);
    }
    setCheckerboard(vecs_res_5d, vecs_res_5d_e);
    setCheckerboard(vecs_res_5d, vecs_res_5d_o);

    doComparison(Dw.Meooe, Dw5.Meooe, vecs_res_4d, vecs_res_5d, tol);
  }

  // Wilson half cb hopping term daggered = "MeooeDag" ////////////////////////

  {
    for(int rhs=0; rhs<nrhs; rhs++) vecs_res_4d[rhs] = Zero();
    vecs_res_5d = Zero();

    std::vector<LatticeFermion> vecs_src_4d_e(nrhs, UrbGrid);
    std::vector<LatticeFermion> vecs_src_4d_o(nrhs, UrbGrid);
    std::vector<LatticeFermion> vecs_res_4d_e(nrhs, UrbGrid);
    std::vector<LatticeFermion> vecs_res_4d_o(nrhs, UrbGrid);

    for(int rhs=0; rhs<nrhs; rhs++) {
      pickCheckerboard(Even, vecs_src_4d_e[rhs], vecs_src_4d[rhs]);
      pickCheckerboard(Odd,  vecs_src_4d_o[rhs], vecs_src_4d[rhs]);
      pickCheckerboard(Even, vecs_res_4d_e[rhs], vecs_res_4d[rhs]);
      pickCheckerboard(Odd,  vecs_res_4d_o[rhs], vecs_res_4d[rhs]);
    }

    // clang-format off
    LatticeFermion vecs_src_5d_e(FrbGrid); pickCheckerboard(Even, vecs_src_5d_e, vecs_src_5d);
    LatticeFermion vecs_src_5d_o(FrbGrid); pickCheckerboard(Odd,  vecs_src_5d_o, vecs_src_5d);
    LatticeFermion vecs_res_5d_e(FrbGrid); pickCheckerboard(Even, vecs_res_5d_e, vecs_res_5d);
    LatticeFermion vecs_res_5d_o(FrbGrid); pickCheckerboard(Odd,  vecs_res_5d_o, vecs_res_5d);
    // clang-format on

    double flop = volumeUrbGrid * nrhs * wilson_dhop.flop;
    double byte = volumeUrbGrid * nrhs * wilson_dhop.byte();

    BenchmarkFunctionMRHS(Dw.MeooeDag,
                          flop, byte, nIterMin, nSecMin, nrhs,
                          vecs_src_4d_e[rhs], vecs_res_4d_o[rhs]);
    BenchmarkFunctionMRHS(Dw.MeooeDag,
                          flop, byte, nIterMin, nSecMin, nrhs,
                          vecs_src_4d_o[rhs], vecs_res_4d_e[rhs]);

    BenchmarkFunction(Dw5.MeooeDag,
                      flop, byte, nIterMin, nSecMin,
                      vecs_src_5d_e, vecs_res_5d_o);
    BenchmarkFunction(Dw5.MeooeDag,
                      flop, byte, nIterMin, nSecMin,
                      vecs_src_5d_o, vecs_res_5d_e);

    for(int rhs=0; rhs<nrhs; rhs++) {
      setCheckerboard(vecs_res_4d[rhs], vecs_res_4d_e[rhs]);
      setCheckerboard(vecs_res_4d[rhs], vecs_res_4d_o[rhs]);
    }
    setCheckerboard(vecs_res_5d, vecs_res_5d_e);
    setCheckerboard(vecs_res_5d, vecs_res_5d_o);

    doComparison(Dw.MeooeDag, Dw5.MeooeDag, vecs_res_4d, vecs_res_5d, tol);
  }

  // Wilson diagonal term = "Mooee" ///////////////////////////////////////////

  {
    for(int rhs=0; rhs<nrhs; rhs++) vecs_res_4d[rhs] = Zero();
    vecs_res_5d = Zero();

    double flop = volumeUGrid * nrhs * wilson_diag.flop;
    double byte = volumeUGrid * nrhs * wilson_diag.byte();

    BenchmarkFunctionMRHS(Dw.Mooee,
                          flop, byte, nIterMin, nSecMin, nrhs,
                          vecs_src_4d[rhs], vecs_res_4d[rhs]);

    BenchmarkFunction(Dw5.Mooee,
                      flop, byte, nIterMin, nSecMin,
                      vecs_src_5d, vecs_res_5d);

    doComparison(Dw.Mooee, Dw5.Mooee, vecs_res_4d, vecs_res_5d, tol);
  }

  // Wilson diagonal term inverse = "MooeeInv" ////////////////////////////////

  {
    for(int rhs=0; rhs<nrhs; rhs++) vecs_res_4d[rhs] = Zero();
    vecs_res_5d = Zero();

    double flop = volumeUGrid * nrhs * wilson_diag.flop;
    double byte = volumeUGrid * nrhs * wilson_diag.byte();

    BenchmarkFunctionMRHS(Dw.MooeeInv,
                          flop, byte, nIterMin, nSecMin, nrhs,
                          vecs_src_4d[rhs], vecs_res_4d[rhs]);

    BenchmarkFunction(Dw5.MooeeInv,
                      flop, byte, nIterMin, nSecMin,
                      vecs_src_5d, vecs_res_5d);

    doComparison(Dw.MooeeInv, Dw5.MooeeInv, vecs_res_4d, vecs_res_5d, tol);
  }

  // Wilson diagonal term daggered = "MooeeDag" ///////////////////////////////

  {
    for(int rhs=0; rhs<nrhs; rhs++) vecs_res_4d[rhs] = Zero();
    vecs_res_5d = Zero();

    double flop = volumeUGrid * nrhs * wilson_diag.flop;
    double byte = volumeUGrid * nrhs * wilson_diag.byte();

    BenchmarkFunctionMRHS(Dw.MooeeDag,
                          flop, byte, nIterMin, nSecMin, nrhs,
                          vecs_src_4d[rhs], vecs_res_4d[rhs]);

    BenchmarkFunction(Dw5.MooeeDag,
                      flop, byte, nIterMin, nSecMin,
                      vecs_src_5d, vecs_res_5d);

    doComparison(Dw.MooeeDag, Dw5.MooeeDag, vecs_res_4d, vecs_res_5d, tol);
  }

  // Wilson diagonal term inverse daggered = "MooeeInvDag" ////////////////////

  {
    for(int rhs=0; rhs<nrhs; rhs++) vecs_res_4d[rhs] = Zero();
    vecs_res_5d = Zero();

    double flop = volumeUGrid * nrhs * wilson_diag.flop;
    double byte = volumeUGrid * nrhs * wilson_diag.byte();

    BenchmarkFunctionMRHS(Dw.MooeeInvDag,
                          flop, byte, nIterMin, nSecMin, nrhs,
                          vecs_src_4d[rhs], vecs_res_4d[rhs]);

    BenchmarkFunction(Dw5.MooeeInvDag,
                      flop, byte, nIterMin, nSecMin,
                      vecs_src_5d, vecs_res_5d);

    doComparison(Dw.MooeeInvDag, Dw5.MooeeInvDag, vecs_res_4d, vecs_res_5d, tol);
  }

  // Wilson directional term = "Mdir" /////////////////////////////////////////

  {
    for(int rhs=0; rhs<nrhs; rhs++) vecs_res_4d[rhs] = Zero();
    vecs_res_5d = Zero();

    double flop = volumeUGrid * nrhs * wilson_dir.flop;
    double byte = volumeUGrid * nrhs * wilson_dir.byte();

    int dir = 1;
    int disp = +1;

    BenchmarkFunctionMRHS(Dw.Mdir,
                          flop, byte, nIterMin, nSecMin, nrhs,
                          vecs_src_4d[rhs], vecs_res_4d[rhs], dir, disp);

    BenchmarkFunction(Dw5.Mdir,
                      flop, byte, nIterMin, nSecMin,
                      vecs_src_5d, vecs_res_5d, dir, disp);

    doComparison(Dw.Mdir, Dw5.Mdir, vecs_res_4d, vecs_res_5d, tol);
  }

  // Wilson full matrix = "M" /////////////////////////////////////////////////

  {
    for(int rhs=0; rhs<nrhs; rhs++) vecs_res_4d[rhs] = Zero();
    vecs_res_5d = Zero();

    double flop = volumeUGrid * nrhs * wilson_full.flop;
    double byte = volumeUGrid * nrhs * wilson_full.byte();

    BenchmarkFunctionMRHS(Dw.M,
                          flop, byte, nIterMin, nSecMin, nrhs,
                          vecs_src_4d[rhs], vecs_res_4d[rhs]);

    BenchmarkFunction(Dw5.M,
                      flop, byte, nIterMin, nSecMin,
                      vecs_src_5d, vecs_res_5d);

    doComparison(Dw.M, Dw5.M, vecs_res_4d, vecs_res_5d, tol);
  }

  // Wilson full matrix daggered = "Mdag" /////////////////////////////////////

  {
    for(int rhs=0; rhs<nrhs; rhs++) vecs_res_4d[rhs] = Zero();
    vecs_res_5d = Zero();

    double flop = volumeUGrid * nrhs * wilson_full.flop;
    double byte = volumeUGrid * nrhs * wilson_full.byte();

    BenchmarkFunctionMRHS(Dw.Mdag,
                          flop, byte, nIterMin, nSecMin, nrhs,
                          vecs_src_4d[rhs], vecs_res_4d[rhs]);

    BenchmarkFunction(Dw5.Mdag,
                      flop, byte, nIterMin, nSecMin,
                      vecs_src_5d, vecs_res_5d);

    doComparison(Dw.Mdag, Dw5.Mdag, vecs_res_4d, vecs_res_5d, tol);
  }

  /////////////////////////////////////////////////////////////////////////////
  //                      Benchmarks for Clover fermions                     //
  /////////////////////////////////////////////////////////////////////////////

  // Clover hopping term = "dslash" = "Dhop" //////////////////////////////////

  {
    for(int rhs=0; rhs<nrhs; rhs++) vecs_res_4d[rhs] = Zero();
    vecs_res_5d = Zero();

    double flop = volumeUGrid * nrhs * clover_dhop.flop;
    double byte = volumeUGrid * nrhs * clover_dhop.byte();

    BenchmarkFunctionMRHS(Dwc.Dhop,
                          flop, byte, nIterMin, nSecMin, nrhs,
                          vecs_src_4d[rhs], vecs_res_4d[rhs], 0);

    BenchmarkFunction(Dwc5.Dhop,
                      flop, byte, nIterMin, nSecMin,
                      vecs_src_5d, vecs_res_5d, 0);

    doComparison(Dwc.Dhop, Dwc5.Dhop, vecs_res_4d, vecs_res_5d, tol);
  }

  // Clover half cb hopping term = "Meooe" ////////////////////////////////////

  {
    for(int rhs=0; rhs<nrhs; rhs++) vecs_res_4d[rhs] = Zero();
    vecs_res_5d = Zero();

    std::vector<LatticeFermion> vecs_src_4d_e(nrhs, UrbGrid);
    std::vector<LatticeFermion> vecs_src_4d_o(nrhs, UrbGrid);
    std::vector<LatticeFermion> vecs_res_4d_e(nrhs, UrbGrid);
    std::vector<LatticeFermion> vecs_res_4d_o(nrhs, UrbGrid);

    for(int rhs=0; rhs<nrhs; rhs++) {
      pickCheckerboard(Even, vecs_src_4d_e[rhs], vecs_src_4d[rhs]);
      pickCheckerboard(Odd,  vecs_src_4d_o[rhs], vecs_src_4d[rhs]);
      pickCheckerboard(Even, vecs_res_4d_e[rhs], vecs_res_4d[rhs]);
      pickCheckerboard(Odd,  vecs_res_4d_o[rhs], vecs_res_4d[rhs]);
    }

    // clang-format off
    LatticeFermion vecs_src_5d_e(FrbGrid); pickCheckerboard(Even, vecs_src_5d_e, vecs_src_5d);
    LatticeFermion vecs_src_5d_o(FrbGrid); pickCheckerboard(Odd,  vecs_src_5d_o, vecs_src_5d);
    LatticeFermion vecs_res_5d_e(FrbGrid); pickCheckerboard(Even, vecs_res_5d_e, vecs_res_5d);
    LatticeFermion vecs_res_5d_o(FrbGrid); pickCheckerboard(Odd,  vecs_res_5d_o, vecs_res_5d);
    // clang-format on

    double flop = volumeUrbGrid * nrhs * clover_dhop.flop;
    double byte = volumeUrbGrid * nrhs * clover_dhop.byte();

    BenchmarkFunctionMRHS(Dwc.Meooe,
                          flop, byte, nIterMin, nSecMin, nrhs,
                          vecs_src_4d_e[rhs], vecs_res_4d_o[rhs]);
    BenchmarkFunctionMRHS(Dwc.Meooe,
                          flop, byte, nIterMin, nSecMin, nrhs,
                          vecs_src_4d_o[rhs], vecs_res_4d_e[rhs]);

    BenchmarkFunction(Dwc5.Meooe,
                      flop, byte, nIterMin, nSecMin,
                      vecs_src_5d_e, vecs_res_5d_o);
    BenchmarkFunction(Dwc5.Meooe,
                      flop, byte, nIterMin, nSecMin,
                      vecs_src_5d_o, vecs_res_5d_e);

    for(int rhs=0; rhs<nrhs; rhs++) {
      setCheckerboard(vecs_res_4d[rhs], vecs_res_4d_e[rhs]);
      setCheckerboard(vecs_res_4d[rhs], vecs_res_4d_o[rhs]);
    }
    setCheckerboard(vecs_res_5d, vecs_res_5d_e);
    setCheckerboard(vecs_res_5d, vecs_res_5d_o);

    doComparison(Dwc.Meooe, Dwc5.Meooe, vecs_res_4d, vecs_res_5d, tol);
  }

  // Clover half cb hopping term daggered = "MeooeDag" ////////////////////////

  {
    for(int rhs=0; rhs<nrhs; rhs++) vecs_res_4d[rhs] = Zero();
    vecs_res_5d = Zero();

    std::vector<LatticeFermion> vecs_src_4d_e(nrhs, UrbGrid);
    std::vector<LatticeFermion> vecs_src_4d_o(nrhs, UrbGrid);
    std::vector<LatticeFermion> vecs_res_4d_e(nrhs, UrbGrid);
    std::vector<LatticeFermion> vecs_res_4d_o(nrhs, UrbGrid);

    for(int rhs=0; rhs<nrhs; rhs++) {
      pickCheckerboard(Even, vecs_src_4d_e[rhs], vecs_src_4d[rhs]);
      pickCheckerboard(Odd,  vecs_src_4d_o[rhs], vecs_src_4d[rhs]);
      pickCheckerboard(Even, vecs_res_4d_e[rhs], vecs_res_4d[rhs]);
      pickCheckerboard(Odd,  vecs_res_4d_o[rhs], vecs_res_4d[rhs]);
    }

    // clang-format off
    LatticeFermion vecs_src_5d_e(FrbGrid); pickCheckerboard(Even, vecs_src_5d_e, vecs_src_5d);
    LatticeFermion vecs_src_5d_o(FrbGrid); pickCheckerboard(Odd,  vecs_src_5d_o, vecs_src_5d);
    LatticeFermion vecs_res_5d_e(FrbGrid); pickCheckerboard(Even, vecs_res_5d_e, vecs_res_5d);
    LatticeFermion vecs_res_5d_o(FrbGrid); pickCheckerboard(Odd,  vecs_res_5d_o, vecs_res_5d);
    // clang-format on

    double flop = volumeUrbGrid * nrhs * clover_dhop.flop;
    double byte = volumeUrbGrid * nrhs * clover_dhop.byte();

    BenchmarkFunctionMRHS(Dwc.MeooeDag,
                          flop, byte, nIterMin, nSecMin, nrhs,
                          vecs_src_4d_e[rhs], vecs_res_4d_o[rhs]);
    BenchmarkFunctionMRHS(Dwc.MeooeDag,
                          flop, byte, nIterMin, nSecMin, nrhs,
                          vecs_src_4d_o[rhs], vecs_res_4d_e[rhs]);

    BenchmarkFunction(Dwc5.MeooeDag,
                      flop, byte, nIterMin, nSecMin,
                      vecs_src_5d_e, vecs_res_5d_o);
    BenchmarkFunction(Dwc5.MeooeDag,
                      flop, byte, nIterMin, nSecMin,
                      vecs_src_5d_o, vecs_res_5d_e);

    for(int rhs=0; rhs<nrhs; rhs++) {
      setCheckerboard(vecs_res_4d[rhs], vecs_res_4d_e[rhs]);
      setCheckerboard(vecs_res_4d[rhs], vecs_res_4d_o[rhs]);
    }
    setCheckerboard(vecs_res_5d, vecs_res_5d_e);
    setCheckerboard(vecs_res_5d, vecs_res_5d_o);

    doComparison(Dwc.MeooeDag, Dcw5.MeooeDag, vecs_res_4d, vecs_res_5d, tol);
  }

  // Clover diagonal term = "Mooee" ///////////////////////////////////////////

  {
    for(int rhs=0; rhs<nrhs; rhs++) vecs_res_4d[rhs] = Zero();
    vecs_res_5d = Zero();

    double flop = volumeUGrid * nrhs * clover_diag.flop;
    double byte = volumeUGrid * nrhs * clover_diag.byte();

    BenchmarkFunctionMRHS(Dwc.Mooee,
                          flop, byte, nIterMin, nSecMin, nrhs,
                          vecs_src_4d[rhs], vecs_res_4d[rhs]);

    BenchmarkFunction(Dwc5.Mooee,
                      flop, byte, nIterMin, nSecMin,
                      vecs_src_5d, vecs_res_5d);

    doComparison(Dwc.Mooee, Dwc5.Mooee, vecs_res_4d, vecs_res_5d, tol);
  }

  // Clover diagonal term inverse = "MooeeInv" ////////////////////////////////

  {
    for(int rhs=0; rhs<nrhs; rhs++) vecs_res_4d[rhs] = Zero();
    vecs_res_5d = Zero();

    double flop = volumeUGrid * nrhs * clover_diag.flop;
    double byte = volumeUGrid * nrhs * clover_diag.byte();

    BenchmarkFunctionMRHS(Dwc.MooeeInv,
                          flop, byte, nIterMin, nSecMin, nrhs,
                          vecs_src_4d[rhs], vecs_res_4d[rhs]);

    BenchmarkFunction(Dwc5.MooeeInv,
                      flop, byte, nIterMin, nSecMin,
                      vecs_src_5d, vecs_res_5d);

    doComparison(Dwc.MooeeInv, Dwc5.MooeeInv, vecs_res_4d, vecs_res_5d, tol);
  }

  // Clover diagonal term daggered = "MooeeDag" ///////////////////////////////

  {
    for(int rhs=0; rhs<nrhs; rhs++) vecs_res_4d[rhs] = Zero();
    vecs_res_5d = Zero();

    double flop = volumeUGrid * nrhs * clover_diag.flop;
    double byte = volumeUGrid * nrhs * clover_diag.byte();

    BenchmarkFunctionMRHS(Dwc.MooeeDag,
                          flop, byte, nIterMin, nSecMin, nrhs,
                          vecs_src_4d[rhs], vecs_res_4d[rhs]);

    BenchmarkFunction(Dwc5.MooeeDag,
                      flop, byte, nIterMin, nSecMin,
                      vecs_src_5d, vecs_res_5d);

    doComparison(Dwc.MooeeDag, Dwc5.MooeeDag, vecs_res_4d, vecs_res_5d, tol);
  }

  // Clover diagonal term inverse daggered = "MooeeInvDag" ////////////////////

  {
    for(int rhs=0; rhs<nrhs; rhs++) vecs_res_4d[rhs] = Zero();
    vecs_res_5d = Zero();

    double flop = volumeUGrid * nrhs * clover_diag.flop;
    double byte = volumeUGrid * nrhs * clover_diag.byte();

    BenchmarkFunctionMRHS(Dwc.MooeeInvDag,
                          flop, byte, nIterMin, nSecMin, nrhs,
                          vecs_src_4d[rhs], vecs_res_4d[rhs]);

    BenchmarkFunction(Dwc5.MooeeInvDag,
                      flop, byte, nIterMin, nSecMin,
                      vecs_src_5d, vecs_res_5d);

    doComparison(Dwc.MooeeInvDag, Dwc5.MooeeInvDag, vecs_res_4d, vecs_res_5d, tol);
  }

  // Clover directional term = "Mdir" /////////////////////////////////////////

  {
    for(int rhs=0; rhs<nrhs; rhs++) vecs_res_4d[rhs] = Zero();
    vecs_res_5d = Zero();

    double flop = volumeUGrid * nrhs * clover_dir.flop;
    double byte = volumeUGrid * nrhs * clover_dir.byte();

    int dir = 1;
    int disp = +1;

    BenchmarkFunctionMRHS(Dwc.Mdir,
                          flop, byte, nIterMin, nSecMin, nrhs,
                          vecs_src_4d[rhs], vecs_res_4d[rhs], dir, disp);

    BenchmarkFunction(Dwc5.Mdir,
                      flop, byte, nIterMin, nSecMin,
                      vecs_src_5d, vecs_res_5d, dir, disp);

    doComparison(Dwc.Mdir, Dwc5.Mdir, vecs_res_4d, vecs_res_5d, tol);
  }

  // Clover full matrix = "M" /////////////////////////////////////////////////

  {
    for(int rhs=0; rhs<nrhs; rhs++) vecs_res_4d[rhs] = Zero();
    vecs_res_5d = Zero();

    double flop = volumeUGrid * nrhs * clover_full.flop;
    double byte = volumeUGrid * nrhs * clover_full.byte();

    BenchmarkFunctionMRHS(Dwc.M,
                          flop, byte, nIterMin, nSecMin, nrhs,
                          vecs_src_4d[rhs], vecs_res_4d[rhs]);

    BenchmarkFunction(Dwc5.M,
                      flop, byte, nIterMin, nSecMin,
                      vecs_src_5d, vecs_res_5d);

    doComparison(Dwc.M, Dwc5.M, vecs_res_4d, vecs_res_5d, tol);
  }

  // Clover full matrix daggered = "Mdag" /////////////////////////////////////

  {
    for(int rhs = 0; rhs < nrhs; rhs++) vecs_res_4d[rhs] = Zero();
    vecs_res_5d = Zero();

    double flop = volumeUGrid * nrhs * clover_full.flop;
    double byte = volumeUGrid * nrhs * clover_full.byte();

    BenchmarkFunctionMRHS(Dwc.Mdag,
                          flop, byte, nIterMin, nSecMin, nrhs,
                          vecs_src_4d[rhs], vecs_res_4d[rhs]);

    BenchmarkFunction(Dwc5.Mdag,
                      flop, byte, nIterMin, nSecMin,
                      vecs_src_5d, vecs_res_5d);

    doComparison(Dwc.Mdag, Dwc5.Mdag, vecs_res_4d, vecs_res_5d, tol);
  }

  Grid_finalize();
}
