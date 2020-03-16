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
  Grid_init(&argc, &argv);

  GridCartesian* grid = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(),
                                                       GridDefaultSimd(Nd, vComplexD::Nsimd()),
                                                       GridDefaultMpi());

  LatticeGaugeField Umu(grid);

  FieldMetaData header;

  if(!Grid::GridCmdOptionExists(argv, argv + argc, "--config")) {
    std::cout << GridLogError << "You need to use --config /path/to/openqcd_config" << std::endl;
    abort();
  }

  std::string file = Grid::GridCmdOptionPayload(argv, argv + argc, "--config");
  assert(!file.empty());

  OpenQcdIO::readConfiguration(Umu, header, file);

  Grid_finalize();
}
