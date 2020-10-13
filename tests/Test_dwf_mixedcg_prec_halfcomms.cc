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

  const int Ls=24;

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
  LatticeFermionD result(FGrid); result=Zero();
  LatticeGaugeFieldD Umu(UGrid);
  LatticeGaugeFieldF Umu_f(UGrid_f); 
  
  SU<Nc>::HotConfiguration(RNG4,Umu);

  precisionChange(Umu_f,Umu);
  
  RealD mass=0.1;
  RealD M5=1.8;
  DomainWallFermionD Ddwf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);
  DomainWallFermionFH Ddwf_f(Umu_f,*FGrid_f,*FrbGrid_f,*UGrid_f,*UrbGrid_f,mass,M5);


  LatticeFermionD    src_o(FrbGrid);
  LatticeFermionD result_cg(FrbGrid);
  pickCheckerboard(Odd,src_o,src);
  result_cg.Checkerboard() = Odd;
  result_cg = Zero();
  LatticeFermionD result_mcg(result_cg);
  LatticeFermionD result_rlcg(result_cg);

  SchurDiagMooeeOperator<DomainWallFermionD,LatticeFermionD> HermOpEO(Ddwf);
  SchurDiagMooeeOperator<DomainWallFermionFH,LatticeFermionF> HermOpEO_f(Ddwf_f);

  //#define DO_MIXED_CG
#define DO_RLUP_CG

#ifdef DO_MIXED_CG
  std::cout << "Starting mixed CG" << std::endl;
  MixedPrecisionConjugateGradient<LatticeFermionD,LatticeFermionF> mCG(1.0e-8, 10000, 50, FrbGrid_f, HermOpEO_f, HermOpEO);
  mCG.InnerTolerance = 3.0e-5;
  mCG(src_o,result_mcg);
#endif

#ifdef DO_RLUP_CG
  std::cout << "Starting reliable update CG" << std::endl;
  ConjugateGradientReliableUpdate<LatticeFermionD,LatticeFermionF> rlCG(1.e-8, 10000, 0.1, FrbGrid_f, HermOpEO_f, HermOpEO);
  rlCG(src_o,result_rlcg);
#endif
  
  std::cout << "Starting regular CG" << std::endl;
  ConjugateGradient<LatticeFermionD> CG(1.0e-8,10000);
  CG(HermOpEO,src_o,result_cg);

#ifdef DO_MIXED_CG
  LatticeFermionD diff_mcg(FrbGrid);
  RealD vdiff_mcg = axpy_norm(diff_mcg, -1.0, result_cg, result_mcg);
  std::cout << "Diff between mixed and regular CG: " << vdiff_mcg << std::endl;
#endif

#ifdef DO_RLUP_CG
  LatticeFermionD diff_rlcg(FrbGrid);
  RealD vdiff_rlcg = axpy_norm(diff_rlcg, -1.0, result_cg, result_rlcg);
  std::cout << "Diff between reliable update and regular CG: " << vdiff_rlcg << std::endl;
#endif
  
  Grid_finalize();
}
