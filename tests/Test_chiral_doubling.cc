/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./tests/Test_chiral_doubling.cc

    Copyright (C) 2015 - 2020

    Author: Daniel Richtmann <daniel.richtmann@gmail.com>

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

  /////////////////////////////////////////////////////////////////////////////
  //                              General setup                              //
  /////////////////////////////////////////////////////////////////////////////

  GridCartesian* Grid = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd, vComplex::Nsimd()), GridDefaultMpi());

  std::cout << GridLogMessage << "Grid:" << std::endl; Grid->show_decomposition();

  GridParallelRNG pRNG(Grid);
  std::vector<int> seeds({1, 2, 3, 4});
  pRNG.SeedFixedIntegers(seeds);

  /////////////////////////////////////////////////////////////////////////////
  //                             The actual tests                            //
  /////////////////////////////////////////////////////////////////////////////

  const int nBasis = 40;
  const int nB     = nBasis / 2;

  std::vector<LatticeFermion> basisVecsRef(nB, Grid);
  std::vector<LatticeFermion> basisVecsRes(nBasis, Grid);

  for(int n=0; n<nB; n++) {
    random(pRNG, basisVecsRef[n]);
    basisVecsRes[n] = basisVecsRef[n];
  }

  performChiralDoublingG5C(basisVecsRes);
  undoChiralDoublingG5C(basisVecsRes);

  LatticeFermion diff(Grid);
  for(int n=0; n<nB; n++) {
    diff       = basisVecsRef[n] - basisVecsRes[n];
    auto diff2 = norm2(diff);
    std::cout << Grid::GridLogMessage << "Vector " << n << ", diff = " << diff2 << std::endl;
    assert(diff2 == 0.);
  }

  std::cout << Grid::GridLogMessage << "Test passed" << std::endl;

  Grid_finalize();
}
