    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_wilson_cg_schur.cc

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

    See the full license in the file "LICENSE" in the top level distribution directory
    *************************************************************************************/
    /*  END LEGAL */
#include <Grid/Grid.h>

using namespace std;
using namespace Grid;
 ;

template<class d>
struct scal {
  d internal;
};

  Gamma::Algebra Gmu [] = {
    Gamma::Algebra::GammaX,
    Gamma::Algebra::GammaY,
    Gamma::Algebra::GammaZ,
    Gamma::Algebra::GammaT
  };

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);


  Coordinate latt_size   = GridDefaultLatt();
  Coordinate simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
  Coordinate mpi_layout  = GridDefaultMpi();
  GridCartesian               Grid(latt_size,simd_layout,mpi_layout);
  GridRedBlackCartesian     RBGrid(&Grid);

  std::vector<int> seeds({1,2,3,4});
  GridParallelRNG          pRNG(&Grid);  pRNG.SeedFixedIntegers(seeds);

  LatticeGaugeField Umu(&Grid); SU<Nc>::HotConfiguration(pRNG,Umu);

  LatticeFermion    src(&Grid); random(pRNG,src);
  LatticeFermion result(&Grid);
  LatticeFermion resid(&Grid);

  ConjugateGradient<LatticeFermion> CG(1.0e-8,10000);
  SchurRedBlackDiagMooeeSolve<LatticeFermion> SchurSolver(CG);

  RealD mass = -0.1;
  RealD csw_r = 1.0;
  RealD csw_t = 1.0;
  RealD cF = 1.0;

  std::cout << GridLogMessage << "Testing Wilson Clover" << std::endl;
  WilsonCloverFermionR Dw(Umu, Grid, RBGrid, mass, csw_r, csw_t);
  result=Zero();
  SchurSolver(Dw,src,result);

  std::cout << GridLogMessage << "Testing Compact Wilson Clover" << std::endl;
  CompactWilsonCloverFermionR Dw_compact(Umu, Grid, RBGrid, mass, csw_r, csw_t, 0.0);
  result=Zero();
  SchurSolver(Dw_compact,src,result);

  std::cout << GridLogMessage << "Testing Wilson Exp Clover" << std::endl;
  WilsonExpCloverFermionR Dwe(Umu, Grid, RBGrid, mass, csw_r, csw_t);
  result=Zero();
  SchurSolver(Dwe,src,result);

  std::cout << GridLogMessage << "Testing Compact Wilson Exp Clover" << std::endl;
  CompactWilsonExpCloverFermionR Dwe_compact(Umu, Grid, RBGrid, mass, csw_r, csw_t, 0.0);
  result=Zero();
  SchurSolver(Dwe_compact,src,result);
  
  Grid_finalize();
}
