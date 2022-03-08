    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_dwf_cr_unprec.cc

    Copyright (C) 2019

Author: Peter Boyle <pboyle@bnl.gov>

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
#include <Grid/algorithms/iterative/QuasiMinimalResidual.h>

using namespace std;
using namespace Grid;


int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  GridCartesian         * Grid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * rbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(Grid);


  std::vector<int> seeds4({1,2,3,4});
  GridParallelRNG          RNG4(Grid);  RNG4.SeedFixedIntegers(seeds4);

  LatticeFermion    src(Grid); random(RNG4,src);
  LatticeFermion result(Grid); result=Zero();
  LatticeGaugeField Umu(Grid); SU<Nc>::HotConfiguration(RNG4,Umu);

  std::vector<LatticeColourMatrix> U(4,Grid);

  for(int mu=0;mu<Nd;mu++){
    U[mu] = PeekIndex<LorentzIndex>(Umu,mu);
  }

  QuasiMinimalResidual<LatticeFermion> QMR(1.0e-8,10000);
  
  RealD mass=0.0;
  WilsonFermionR Dw(Umu,*Grid,*rbGrid,mass);

  NonHermitianLinearOperator<WilsonFermionR,LatticeFermion> NonHermOp(Dw);
  QMR(NonHermOp,src,result);

  Grid_finalize();
}
