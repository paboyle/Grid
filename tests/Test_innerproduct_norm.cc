/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./tests/Test_innerproduct_norm.cc

Copyright (C) 2015

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

using namespace Grid;

int main(int argc, char** argv) {
  Grid_init(&argc, &argv);

  const int nIter = 100;

  // clang-format off
  GridCartesian *Grid_d = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd, vComplexD::Nsimd()), GridDefaultMpi());
  GridCartesian *Grid_f = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd, vComplexF::Nsimd()), GridDefaultMpi());
  // clang-format on

  GridParallelRNG pRNG_d(Grid_d);
  GridParallelRNG pRNG_f(Grid_f);

  std::vector<int> seeds_d({1, 2, 3, 4});
  std::vector<int> seeds_f({5, 6, 7, 8});

  pRNG_d.SeedFixedIntegers(seeds_d);
  pRNG_f.SeedFixedIntegers(seeds_f);

  // clang-format off
  LatticeFermionD x_d(Grid_d); random(pRNG_d, x_d);
  LatticeFermionD y_d(Grid_d); random(pRNG_d, y_d);
  LatticeFermionF x_f(Grid_f); random(pRNG_f, x_f);
  LatticeFermionF y_f(Grid_f); random(pRNG_f, y_f);
  // clang-format on

  GridStopWatch sw_ref;
  GridStopWatch sw_res;

  { // double precision
    ComplexD ip_d_ref, ip_d_res, diff_ip_d;
    RealD    norm2_d_ref, norm2_d_res, diff_norm2_d;

    sw_ref.Reset();
    sw_ref.Start();
    for(int i = 0; i < nIter; ++i) {
      ip_d_ref    = innerProduct(x_d, y_d);
      norm2_d_ref = norm2(x_d);
    }
    sw_ref.Stop();

    sw_res.Reset();
    sw_res.Start();
    for(int i = 0; i < nIter; ++i) { innerProduct_norm(ip_d_res, norm2_d_res, x_d, y_d); }
    sw_res.Stop();

    diff_ip_d    = ip_d_ref - ip_d_res;
    diff_norm2_d = norm2_d_ref - norm2_d_res;

    // clang-format off
    std::cout << GridLogMessage << "Double: ip_ref = " << ip_d_ref << " ip_res = " << ip_d_res << " diff = " << diff_ip_d << std::endl;
    std::cout << GridLogMessage << "Double: norm2_ref = " << norm2_d_ref << " norm2_res = " << norm2_d_res << " diff = " << diff_norm2_d << std::endl;
    std::cout << GridLogMessage << "Double: time_ref = " << sw_ref.Elapsed() << " time_res = " << sw_res.Elapsed() << std::endl;
    // clang-format on

    assert(diff_ip_d == 0.);
    assert(diff_norm2_d == 0.);

    std::cout << GridLogMessage << "Double: all checks passed" << std::endl;
  }

  { // single precision
    ComplexD ip_f_ref, ip_f_res, diff_ip_f;
    RealD    norm2_f_ref, norm2_f_res, diff_norm2_f;

    sw_ref.Reset();
    sw_ref.Start();
    for(int i = 0; i < nIter; ++i) {
      ip_f_ref    = innerProduct(x_f, y_f);
      norm2_f_ref = norm2(x_f);
    }
    sw_ref.Stop();

    sw_res.Reset();
    sw_res.Start();
    for(int i = 0; i < nIter; ++i) { innerProduct_norm(ip_f_res, norm2_f_res, x_f, y_f); }
    sw_res.Stop();

    diff_ip_f    = ip_f_ref - ip_f_res;
    diff_norm2_f = norm2_f_ref - norm2_f_res;

    // clang-format off
    std::cout << GridLogMessage << "Single: ip_ref = " << ip_f_ref << " ip_res = " << ip_f_res << " diff = " << diff_ip_f << std::endl;
    std::cout << GridLogMessage << "Single: norm2_ref = " << norm2_f_ref << " norm2_res = " << norm2_f_res << " diff = " << diff_norm2_f << std::endl;
    std::cout << GridLogMessage << "Single: time_ref = " << sw_ref.Elapsed() << " time_res = " << sw_res.Elapsed() << std::endl;
    // clang-format on

    assert(diff_ip_f == 0.);
    assert(diff_norm2_f == 0.);

    std::cout << GridLogMessage << "Single: all checks passed" << std::endl;
  }

  Grid_finalize();
}
