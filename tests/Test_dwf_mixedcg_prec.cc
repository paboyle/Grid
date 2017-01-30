    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_dwf_cg_prec.cc

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
using namespace Grid::QCD;

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

  const int Ls=8;

  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplexD::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  GridCartesian         * UGrid_f   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplexF::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid_f = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid_f);
  GridCartesian         * FGrid_f   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid_f);
  GridRedBlackCartesian * FrbGrid_f = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid_f);
  
  std::vector<int> seeds4({1,2,3,4});
  std::vector<int> seeds5({5,6,7,8});
  GridParallelRNG          RNG5(FGrid);  RNG5.SeedFixedIntegers(seeds5);
  GridParallelRNG          RNG4(UGrid);  RNG4.SeedFixedIntegers(seeds4);

  LatticeFermionD    src(FGrid); random(RNG5,src);
  LatticeFermionD result(FGrid); result=zero;
  LatticeGaugeFieldD Umu(UGrid);
  LatticeGaugeFieldF Umu_f(UGrid_f); 
  
  SU3::HotConfiguration(RNG4,Umu);

  precisionChange(Umu_f,Umu);
  
  RealD mass=0.1;
  RealD M5=1.8;
  DomainWallFermionD Ddwf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);
  DomainWallFermionF Ddwf_f(Umu_f,*FGrid_f,*FrbGrid_f,*UGrid_f,*UrbGrid_f,mass,M5);

  LatticeFermionD    src_o(FrbGrid);
  LatticeFermionD result_o(FrbGrid);
  LatticeFermionD result_o_2(FrbGrid);
  pickCheckerboard(Odd,src_o,src);
  result_o.checkerboard = Odd;
  result_o = zero;
  result_o_2.checkerboard = Odd;
  result_o_2 = zero;

  SchurDiagMooeeOperator<DomainWallFermionD,LatticeFermionD> HermOpEO(Ddwf);
  SchurDiagMooeeOperator<DomainWallFermionF,LatticeFermionF> HermOpEO_f(Ddwf_f);

  std::cout << "Starting mixed CG" << std::endl;
  MixedPrecisionConjugateGradient<LatticeFermionD,LatticeFermionF> mCG(1.0e-8, 10000, 50, FrbGrid_f, HermOpEO_f, HermOpEO);
  mCG(src_o,result_o);

  std::cout << "Starting regular CG" << std::endl;
  ConjugateGradient<LatticeFermionD> CG(1.0e-8,10000);
  CG(HermOpEO,src_o,result_o_2);

  LatticeFermionD diff_o(FrbGrid);
  RealD diff = axpy_norm(diff_o, -1.0, result_o, result_o_2);

  std::cout << "Diff between mixed and regular CG: " << diff << std::endl;

  
  Grid_finalize();
}
