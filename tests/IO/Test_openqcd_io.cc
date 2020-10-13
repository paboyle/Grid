/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./tests/io/Test_openqcd_io.cc

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



using namespace Grid;

int main(int argc, char** argv) {
#if !defined(GRID_COMMS_NONE)
  Grid_init(&argc, &argv);

  auto simd_layout = GridDefaultSimd(Nd, vComplex::Nsimd());
  auto mpi_layout  = GridDefaultMpi();
  auto latt_size   = GridDefaultLatt();

  GridCartesian grid(latt_size, simd_layout, mpi_layout);

  GridParallelRNG pRNG(&grid);

  pRNG.SeedFixedIntegers(std::vector<int>({45, 12, 81, 9}));

  LatticeGaugeField Umu_ref(&grid);
  LatticeGaugeField Umu_me(&grid);
  LatticeGaugeField Umu_diff(&grid);

  FieldMetaData header_ref;
  FieldMetaData header_me;

  Umu_ref = Zero();
  Umu_me  = Zero();

  std::string file("/home/daniel/configs/openqcd/test_16x8_pbcn6");

  if(GridCmdOptionExists(argv, argv + argc, "--config")) {
    file = GridCmdOptionPayload(argv, argv + argc, "--config");
    std::cout << "file: " << file << std::endl;
    assert(!file.empty());
  }

  OpenQcdIOChromaReference::readConfiguration(Umu_ref, header_ref, file);
  OpenQcdIO::readConfiguration(Umu_me, header_me, file);

  std::cout << GridLogMessage << header_ref << std::endl;
  std::cout << GridLogMessage << header_me << std::endl;

  Umu_diff = Umu_ref - Umu_me;

  // clang-format off
  std::cout << GridLogMessage
            << "norm2(Umu_ref) = " << norm2(Umu_ref)
            << " norm2(Umu_me) = " << norm2(Umu_me)
            << " norm2(Umu_diff) = " << norm2(Umu_diff) << std::endl;
  // clang-format on

  Grid_finalize();
#endif
}
