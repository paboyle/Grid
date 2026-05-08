/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: tests/core/Test_plaquette_stats.cc

Copyright (C) 2015

Author: Chulwoo Jung <chulwoo@bnl.gov>

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

/**
 * Test_plaquette_stats
 *
 * Measure every plaquette Re Tr[U_mu(x) U_nu(x+mu) U_mu†(x+nu) U_nu†(x)] / Nc
 * across all sites and all (mu,nu) planes and report max, min, and average.
 *
 * Usage:
 *   ./Test_plaquette_stats [Grid options] [--file <nersc_config>] [--hot]
 *                          [--threshold <val>]
 *
 *   --file <path>       Read gauge field from NERSC-format file
 *   --hot               Use random (hot) SU(3) start (default: cold/unit start)
 *   --threshold <val>   Print coordinates of every plaquette with Re Tr/Nc < val
 *
 * Grid size defaults to 8^3 x 16; override with --grid (e.g. --grid 4.4.4.8).
 */

#include <Grid/Grid.h>
#include <limits>

using namespace std;
using namespace Grid;

// Return the plane label string for output
static std::string planeName(int mu, int nu)
{
  const char dirs[] = "xyzt";
  std::string s;
  s += dirs[mu];
  s += dirs[nu];
  return s;
}

int main(int argc, char** argv)
{
  Grid_init(&argc, &argv);

  // ---- Lattice geometry ----
  Coordinate latt_size   = GridDefaultLatt();
  if (latt_size.size() == 0) {
    latt_size = Coordinate(std::vector<int>{8, 8, 8, 16});
  }
  Coordinate simd_layout = GridDefaultSimd(Nd, vComplex::Nsimd());
  Coordinate mpi_layout  = GridDefaultMpi();

  GridCartesian grid(latt_size, simd_layout, mpi_layout);

  std::cout << GridLogMessage << "Lattice: ";
  for (int d = 0; d < Nd; d++)
    std::cout << latt_size[d] << (d < Nd-1 ? "x" : "\n");

  // ---- Gauge field ----
  LatticeGaugeField Umu(&grid);

  // Check for --file argument
  std::string config_file = "";
  for (int i = 1; i < argc - 1; i++) {
    if (std::string(argv[i]) == "--file") {
      config_file = argv[i+1];
      break;
    }
  }

  bool doHot = false;
  for (int i = 1; i < argc; i++) {
    if (std::string(argv[i]) == "--hot") { doHot = true; break; }
  }

  double threshold = std::numeric_limits<double>::quiet_NaN();
  bool   doThresh  = false;
  for (int i = 1; i < argc - 1; i++) {
    if (std::string(argv[i]) == "--threshold") {
      threshold = std::stod(argv[i+1]);
      doThresh  = true;
      break;
    }
  }

  if (!config_file.empty()) {
    std::cout << GridLogMessage << "Reading gauge field from " << config_file << std::endl;
    FieldMetaData header;
    NerscIO::readConfiguration(Umu, header, config_file);
  } else {
    std::vector<int> seeds({1, 2, 3, 4});
    GridParallelRNG pRNG(&grid);
    pRNG.SeedFixedIntegers(seeds);
    if (doHot) {
      std::cout << GridLogMessage << "Generating hot (random SU(3)) start" << std::endl;
      SU<Nc>::HotConfiguration(pRNG, Umu);
    } else {
      std::cout << GridLogMessage << "Using cold (unit) gauge start" << std::endl;
      SU<Nc>::ColdConfiguration(pRNG, Umu);
    }
  }

  // ---- Extract link matrices ----
  std::vector<LatticeColourMatrix> U(Nd, &grid);
  for (int mu = 0; mu < Nd; mu++)
    U[mu] = PeekIndex<LorentzIndex>(Umu, mu);

  // ---- Per-plane plaquette statistics ----
  //
  // For each (mu, nu) plane (mu > nu) compute
  //   P_munu(x) = Re Tr[plaquette] / Nc
  // then report max, min, mean over all sites.
  //
  // WilsonLoops::traceDirPlaquette gives Tr[U_mu U_nu U_mu† U_nu†] (complex).

  double vol = grid.gSites();

  // Accumulate site-average plaquette (sum over planes / Nplanes / Nc)
  LatticeComplex plaq_all(&grid);
  plaq_all = Zero();
  int Nplanes = 0;

  std::cout << GridLogMessage
            << "======== Per-plane plaquette statistics ========" << std::endl;
  std::cout << GridLogMessage
            << std::setw(6)  << "plane"
            << std::setw(20) << "max"
            << std::setw(20) << "min"
            << std::setw(20) << "average"
            << std::endl;

  for (int mu = 1; mu < Nd; mu++) {
    for (int nu = 0; nu < mu; nu++) {

      // Per-site trace of plaquette in (mu,nu) plane
      LatticeComplex sitePlaq(&grid);
      ColourWilsonLoops::traceDirPlaquette(sitePlaq, U, mu, nu);

      plaq_all = plaq_all + sitePlaq;
      Nplanes++;

      // --- global average via sum() ---
      TComplex Tsum  = sum(sitePlaq);
      ComplexD csum  = TensorRemove(Tsum);
      RealD    avg   = csum.real() / vol / Nc;

      // --- global max and min via unvectorize + GlobalMax/GlobalMin ---
      std::vector<TComplex> sv;
      unvectorizeToLexOrdArray(sv, sitePlaq);

      RealD local_max = -1e38, local_min = 1e38;
      for (auto& tc : sv) {
        RealD r = TensorRemove(tc).real() / Nc;
        if (r > local_max) local_max = r;
        if (r < local_min) local_min = r;
      }

      // Reduce across MPI ranks (no GlobalMin; negate to use GlobalMax)
      grid.GlobalMax(local_max);
      local_min = -local_min;
      grid.GlobalMax(local_min);
      local_min = -local_min;

      std::cout << GridLogMessage
                << std::setw(6)  << planeName(mu, nu)
                << std::setw(20) << std::setprecision(10) << local_max
                << std::setw(20) << std::setprecision(10) << local_min
                << std::setw(20) << std::setprecision(10) << avg
                << std::endl;

      if (doThresh) {
        for (int s = 0; s < (int)sv.size(); s++) {
          RealD val = TensorRemove(sv[s]).real() / Nc;
          if (val < threshold) {
            Coordinate lc(Nd), gc(Nd);
            Lexicographic::CoorFromIndex(lc, s, grid._ldimensions);
            for (int d = 0; d < Nd; d++)
              gc[d] = grid._processor_coor[d] * grid._ldimensions[d] + lc[d];
            std::cout << "BELOW_THRESHOLD  plane=" << planeName(mu, nu)
                      << "  site=(" << gc[0] << "," << gc[1] << "," << gc[2] << "," << gc[3] << ")"
                      << "  P=" << std::setprecision(10) << val << std::endl;
          }
        }
      }
    }
  }

  // ---- Overall (averaged over all planes) statistics ----
  // plaq_all = sum of Tr[...] over all 6 (mu,nu) planes
  // Normalise to Re Tr / Nc per plane
  plaq_all = plaq_all * (1.0 / Nc / Nplanes);

  TComplex Tsum_all = sum(plaq_all);
  RealD    avg_all  = TensorRemove(Tsum_all).real() / vol;

  std::vector<TComplex> sv_all;
  unvectorizeToLexOrdArray(sv_all, plaq_all);

  RealD max_all = -1e38, min_all = 1e38;
  for (auto& tc : sv_all) {
    RealD r = TensorRemove(tc).real();
    if (r > max_all) max_all = r;
    if (r < min_all) min_all = r;
  }
  grid.GlobalMax(max_all);
  min_all = -min_all;
  grid.GlobalMax(min_all);
  min_all = -min_all;

  // Cross-check with built-in avgPlaquette
  RealD avg_builtin = ColourWilsonLoops::avgPlaquette(Umu);

  std::cout << GridLogMessage
            << "======== Overall plaquette statistics (all planes) ========" << std::endl;
  std::cout << GridLogMessage << "  max     = " << max_all  << std::endl;
  std::cout << GridLogMessage << "  min     = " << min_all  << std::endl;
  std::cout << GridLogMessage << "  average = " << avg_all  << std::endl;
  std::cout << GridLogMessage << "  avgPlaquette (builtin check) = " << avg_builtin << std::endl;

  Grid_finalize();
  return 0;
}
