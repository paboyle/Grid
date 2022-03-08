/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./tests/core/Test_compact_wilson_clover_speedup.cc

    Copyright (C) 2020 - 2022

    Author: Daniel Richtmann <daniel.richtmann@gmail.com>
    Author: Nils Meyer       <nils.meyer@ur.de>

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

using namespace Grid;

NAMESPACE_BEGIN(CommandlineHelpers);

static bool checkPresent(int* argc, char*** argv, const std::string& option) {
  return GridCmdOptionExists(*argv, *argv + *argc, option);
}

static std::string getContent(int* argc, char*** argv, const std::string& option) {
  return GridCmdOptionPayload(*argv, *argv + *argc, option);
}

static int readInt(int* argc, char*** argv, std::string&& option, int defaultValue) {
  std::string arg;
  int         ret = defaultValue;
  if(checkPresent(argc, argv, option)) {
    arg = getContent(argc, argv, option);
    GridCmdOptionInt(arg, ret);
  }
  return ret;
}

static float readFloat(int* argc, char*** argv, std::string&& option, float defaultValue) {
  std::string arg;
  float       ret = defaultValue;
  if(checkPresent(argc, argv, option)) {
    arg = getContent(argc, argv, option);
    GridCmdOptionFloat(arg, ret);
  }
  return ret;
}

NAMESPACE_END(CommandlineHelpers);


#define _grid_printf(LOGGER, ...) \
  { \
    if((LOGGER).isActive()) { /* this makes it safe to put, e.g., norm2 in the calling code w.r.t. performance */ \
      char _printf_buf[1024]; \
      std::sprintf(_printf_buf, __VA_ARGS__); \
      std::cout << (LOGGER) << _printf_buf; \
      fflush(stdout); \
    } \
  }
#define grid_printf_msg(...) _grid_printf(GridLogMessage, __VA_ARGS__)


template<typename Field>
bool resultsAgree(const Field& ref, const Field& res, const std::string& name) {
  RealD checkTolerance = (getPrecision<Field>::value == 2) ? 1e-15 : 1e-7;
  Field diff(ref.Grid());
  diff = ref - res;
  auto absDev = norm2(diff);
  auto relDev = absDev / norm2(ref);
  std::cout << GridLogMessage
            << "norm2(reference), norm2(" << name << "), abs. deviation, rel. deviation: " << norm2(ref) << " "
            << norm2(res) << " " << absDev << " " << relDev << " -> check "
            << ((relDev < checkTolerance) ? "passed" : "failed") << std::endl;

  return relDev <= checkTolerance;
}


template<typename vCoeff_t>
void runBenchmark(int* argc, char*** argv) {
  // read from command line
  const int   nIter        = CommandlineHelpers::readInt(     argc, argv, "--niter", 1000);
  const RealD mass         = CommandlineHelpers::readFloat(   argc, argv, "--mass",  0.5);
  const RealD csw          = CommandlineHelpers::readFloat(   argc, argv, "--csw",   1.0);
  const RealD cF           = CommandlineHelpers::readFloat(   argc, argv, "--cF",    1.0);
  const bool  antiPeriodic = CommandlineHelpers::checkPresent(argc, argv, "--antiperiodic");

  // precision
  static_assert(getPrecision<vCoeff_t>::value == 2 || getPrecision<vCoeff_t>::value == 1, "Incorrect precision"); // double or single
  std::string precision = (getPrecision<vCoeff_t>::value == 2 ? "double" : "single");

  // setup grids
  GridCartesian*         UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd, vCoeff_t::Nsimd()), GridDefaultMpi());
  GridRedBlackCartesian* UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  // clang-format on

  // setup rng
  std::vector<int> seeds({1, 2, 3, 4});
  GridParallelRNG  pRNG(UGrid);
  pRNG.SeedFixedIntegers(seeds);

  // type definitions
  typedef WilsonImpl<vCoeff_t, FundamentalRepresentation, CoeffReal> WImpl;
  typedef WilsonCloverFermion<WImpl> WilsonCloverOperator;
  typedef CompactWilsonCloverFermion<WImpl> CompactWilsonCloverOperator;
  typedef typename WilsonCloverOperator::FermionField Fermion;
  typedef typename WilsonCloverOperator::GaugeField Gauge;

  // setup fields
  Fermion src(UGrid); random(pRNG, src);
  Fermion ref(UGrid); ref = Zero();
  Fermion res(UGrid); res = Zero();
  Fermion hop(UGrid); hop = Zero();
  Fermion diff(UGrid); diff = Zero();
  Gauge   Umu(UGrid); SU3::HotConfiguration(pRNG, Umu);

  // setup boundary phases
  typename WilsonCloverOperator::ImplParams implParams;
  std::vector<Complex> boundary_phases(Nd, 1.);
  if(antiPeriodic) boundary_phases[Nd-1] = -1.;
  implParams.boundary_phases = boundary_phases;
  WilsonAnisotropyCoefficients anisParams;

  // misc stuff needed for benchmarks
  double volume=1.0; for(int mu=0; mu<Nd; mu++) volume*=UGrid->_fdimensions[mu];

  // setup fermion operators
  WilsonCloverOperator        Dwc(        Umu, *UGrid, *UrbGrid, mass, csw, csw,     anisParams, implParams);
  CompactWilsonCloverOperator Dwc_compact(Umu, *UGrid, *UrbGrid, mass, csw, csw, cF, anisParams, implParams);

  // now test the conversions
  typename CompactWilsonCloverOperator::CloverField         tmp_ref(UGrid);  tmp_ref  = Dwc.CloverTerm;
  typename CompactWilsonCloverOperator::CloverField         tmp_res(UGrid);  tmp_res  = Zero();
  typename CompactWilsonCloverOperator::CloverField         tmp_diff(UGrid); tmp_diff = Zero();
  typename CompactWilsonCloverOperator::CloverDiagonalField diagonal(UGrid); diagonal = Zero();
  typename CompactWilsonCloverOperator::CloverTriangleField triangle(UGrid); diagonal = Zero();
  CompactWilsonCloverOperator::CompactHelpers::ConvertLayout(tmp_ref, diagonal, triangle);
  CompactWilsonCloverOperator::CompactHelpers::ConvertLayout(diagonal, triangle, tmp_res);
  tmp_diff = tmp_ref - tmp_res;
  std::cout << GridLogMessage << "conversion: ref, res, diff, eps"
            << " " << norm2(tmp_ref)
            << " " << norm2(tmp_res)
            << " " << norm2(tmp_diff)
            << " " << norm2(tmp_diff) / norm2(tmp_ref)
            << std::endl;

  // performance per site (use minimal values necessary)
  double hop_flop_per_site            = 1320; // Rich's Talk + what Peter uses
  double hop_byte_per_site            = (8 * 9 + 9 * 12) * 2 * getPrecision<vCoeff_t>::value * 4;
  double clov_flop_per_site           = 504; // Rich's Talk and 1412.2629
  double clov_byte_per_site           = (2 * 18 + 12 + 12) * 2 * getPrecision<vCoeff_t>::value * 4;
  double clov_flop_per_site_performed = 1128;
  double clov_byte_per_site_performed = (12 * 12 + 12 + 12) * 2 * getPrecision<vCoeff_t>::value * 4;

  // total performance numbers
  double hop_gflop_total            = volume * nIter * hop_flop_per_site / 1e9;
  double hop_gbyte_total            = volume * nIter * hop_byte_per_site / 1e9;
  double clov_gflop_total           = volume * nIter * clov_flop_per_site / 1e9;
  double clov_gbyte_total           = volume * nIter * clov_byte_per_site / 1e9;
  double clov_gflop_performed_total = volume * nIter * clov_flop_per_site_performed / 1e9;
  double clov_gbyte_performed_total = volume * nIter * clov_byte_per_site_performed / 1e9;

  // warmup + measure dhop
  for(auto n : {1, 2, 3, 4, 5}) Dwc.Dhop(src, hop, 0);
  double t0 = usecond();
  for(int n = 0; n < nIter; n++) Dwc.Dhop(src, hop, 0);
  double t1 = usecond();
  double secs_hop = (t1-t0)/1e6;
  grid_printf_msg("Performance(%35s, %s): %2.4f s, %6.0f GFlop/s, %6.0f GByte/s, speedup vs ref = %.2f, fraction of hop = %.2f\n",
              "hop", precision.c_str(), secs_hop, hop_gflop_total/secs_hop, hop_gbyte_total/secs_hop, 0.0, secs_hop/secs_hop);

#define BENCH_CLOVER_KERNEL(KERNEL) { \
  /* warmup + measure reference clover */ \
  for(auto n : {1, 2, 3, 4, 5}) Dwc.KERNEL(src, ref); \
  double t2 = usecond(); \
  for(int n = 0; n < nIter; n++) Dwc.KERNEL(src, ref); \
  double t3 = usecond(); \
  double secs_ref = (t3-t2)/1e6; \
  grid_printf_msg("Performance(%35s, %s): %2.4f s, %6.0f GFlop/s, %6.0f GByte/s, speedup vs ref = %.2f, fraction of hop = %.2f\n", \
                  "reference_"#KERNEL, precision.c_str(), secs_ref, clov_gflop_total/secs_ref, clov_gbyte_total/secs_ref, secs_ref/secs_ref, secs_ref/secs_hop); \
  grid_printf_msg("Performance(%35s, %s): %2.4f s, %6.0f GFlop/s, %6.0f GByte/s, speedup vs ref = %.2f, fraction of hop = %.2f\n", /* to see how well the ET performs */  \
                  "reference_"#KERNEL"_performed", precision.c_str(), secs_ref, clov_gflop_performed_total/secs_ref, clov_gbyte_performed_total/secs_ref, secs_ref/secs_ref, secs_ref/secs_hop); \
\
  /* warmup + measure compact clover */ \
  for(auto n : {1, 2, 3, 4, 5}) Dwc_compact.KERNEL(src, res); \
  double t4 = usecond(); \
  for(int n = 0; n < nIter; n++) Dwc_compact.KERNEL(src, res); \
  double t5 = usecond(); \
  double secs_res = (t5-t4)/1e6; \
  grid_printf_msg("Performance(%35s, %s): %2.4f s, %6.0f GFlop/s, %6.0f GByte/s, speedup vs ref = %.2f, fraction of hop = %.2f\n", \
                  "compact_"#KERNEL, precision.c_str(), secs_res, clov_gflop_total/secs_res, clov_gbyte_total/secs_res, secs_ref/secs_res, secs_res/secs_hop); \
  assert(resultsAgree(ref, res, #KERNEL)); \
}

  BENCH_CLOVER_KERNEL(Mooee);
  BENCH_CLOVER_KERNEL(MooeeDag);
  BENCH_CLOVER_KERNEL(MooeeInv);
  BENCH_CLOVER_KERNEL(MooeeInvDag);

  grid_printf_msg("finalize %s\n", precision.c_str());
}

int main(int argc, char** argv) {
  Grid_init(&argc, &argv);

  runBenchmark<vComplexD>(&argc, &argv);
  runBenchmark<vComplexF>(&argc, &argv);

  Grid_finalize();
}
