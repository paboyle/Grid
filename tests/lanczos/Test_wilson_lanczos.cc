/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./tests/Test_dwf_lanczos.cc

Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>

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
#include <Grid/Grid.h>

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

typedef WilsonFermionR FermionOp;
typedef typename WilsonFermionR::FermionField FermionField;


RealD AllZero(RealD x) { return 0.; }

int main(int argc, char** argv) {
  Grid_init(&argc, &argv);

  GridCartesian* UGrid = SpaceTimeGrid::makeFourDimGrid(
      GridDefaultLatt(), GridDefaultSimd(Nd, vComplex::Nsimd()),
      GridDefaultMpi());
  GridRedBlackCartesian* UrbGrid =
      SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian* FGrid = UGrid;
  GridRedBlackCartesian* FrbGrid = UrbGrid;
  printf("UGrid=%p UrbGrid=%p FGrid=%p FrbGrid=%p\n", UGrid, UrbGrid, FGrid,
         FrbGrid);

  std::vector<int> seeds4({1, 2, 3, 4});
  std::vector<int> seeds5({5, 6, 7, 8});
  GridParallelRNG RNG5(FGrid);
  RNG5.SeedFixedIntegers(seeds5);
  GridParallelRNG RNG4(UGrid);
  RNG4.SeedFixedIntegers(seeds4);
  GridParallelRNG RNG5rb(FrbGrid);
  RNG5.SeedFixedIntegers(seeds5);

  LatticeGaugeField Umu(UGrid);
  SU3::HotConfiguration(RNG4, Umu);

/*
  std::vector<LatticeColourMatrix> U(4, UGrid);
  for (int mu = 0; mu < Nd; mu++) {
    U[mu] = PeekIndex<LorentzIndex>(Umu, mu);
  }
*/

  RealD mass = -0.1;
  RealD M5 = 1.8;
  RealD mob_b = 1.5;
  FermionOp WilsonOperator(Umu,*FGrid,*FrbGrid,mass);
  MdagMLinearOperator<FermionOp,LatticeFermion> HermOp(WilsonOperator); /// <-----
  //SchurDiagTwoOperator<FermionOp,FermionField> HermOp(WilsonOperator);

  const int Nstop = 20;
  const int Nk = 60;
  const int Np = 60;
  const int Nm = Nk + Np;
  const int MaxIt = 10000;
  RealD resid = 1.0e-6;

  std::vector<double> Coeffs{0, 1.};
  Polynomial<FermionField> PolyX(Coeffs);
  Chebyshev<FermionField> Cheby(0.0, 10., 12);

  FunctionHermOp<FermionField> OpCheby(Cheby,HermOp);
     PlainHermOp<FermionField> Op     (HermOp);

  ImplicitlyRestartedLanczos<FermionField> IRL(OpCheby, Op, Nstop, Nk, Nm, resid, MaxIt);

  std::vector<RealD> eval(Nm);
  FermionField src(FGrid);
  gaussian(RNG5, src);
  std::vector<FermionField> evec(Nm, FGrid);
  for (int i = 0; i < 1; i++) {
    std::cout << i << " / " << Nm << " grid pointer " << evec[i]._grid
              << std::endl;
  };

  int Nconv;
  IRL.calc(eval, evec, src, Nconv);

  std::cout << eval << std::endl;

  Grid_finalize();
}
