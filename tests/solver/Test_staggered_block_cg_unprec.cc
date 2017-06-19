    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_wilson_cg_unprec.cc

    Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
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
#include <Grid/algorithms/iterative/BlockConjugateGradient.h>

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
  typedef typename ImprovedStaggeredFermion5DR::FermionField FermionField; 
  typedef typename ImprovedStaggeredFermion5DR::ComplexField ComplexField; 
  typename ImprovedStaggeredFermion5DR::ImplParams params; 

  const int Ls=4;

  Grid_init(&argc,&argv);

  std::vector<int> latt_size   = GridDefaultLatt();
  std::vector<int> simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
  std::vector<int> mpi_layout  = GridDefaultMpi();

  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  std::vector<int> seeds({1,2,3,4});
  GridParallelRNG pRNG(UGrid );  pRNG.SeedFixedIntegers(seeds);
  GridParallelRNG pRNG5(FGrid);  pRNG5.SeedFixedIntegers(seeds);

  FermionField src(FGrid); random(pRNG5,src);
  FermionField result(FGrid); result=zero;
  RealD nrm = norm2(src);

  LatticeGaugeField Umu(UGrid); SU3::HotConfiguration(pRNG,Umu);

  RealD mass=0.003;
  ImprovedStaggeredFermion5DR Ds(Umu,Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass);
  MdagMLinearOperator<ImprovedStaggeredFermion5DR,FermionField> HermOp(Ds);

  ConjugateGradient<FermionField> CG(1.0e-8,10000);
  int blockDim = 0;
  BlockConjugateGradient<FermionField>    BCG(blockDim,1.0e-8,10000);
  MultiRHSConjugateGradient<FermionField> mCG(blockDim,1.0e-8,10000);

  std::cout << GridLogMessage << "************************************************************************ "<<std::endl;
  std::cout << GridLogMessage << " Calling 4d CG "<<std::endl;
  std::cout << GridLogMessage << "************************************************************************ "<<std::endl;
  ImprovedStaggeredFermionR Ds4d(Umu,Umu,*UGrid,*UrbGrid,mass);
  MdagMLinearOperator<ImprovedStaggeredFermionR,FermionField> HermOp4d(Ds4d);
  FermionField src4d(UGrid); random(pRNG,src4d);
  FermionField result4d(UGrid); result4d=zero;
  CG(HermOp4d,src4d,result4d);
  std::cout << GridLogMessage << "************************************************************************ "<<std::endl;


  std::cout << GridLogMessage << "************************************************************************ "<<std::endl;
  std::cout << GridLogMessage << " Calling 5d CG for "<<Ls <<" right hand sides" <<std::endl;
  std::cout << GridLogMessage << "************************************************************************ "<<std::endl;
  result=zero;
  CG(HermOp,src,result);
  std::cout << GridLogMessage << "************************************************************************ "<<std::endl;

  std::cout << GridLogMessage << "************************************************************************ "<<std::endl;
  std::cout << GridLogMessage << " Calling multiRHS CG for "<<Ls <<" right hand sides" <<std::endl;
  std::cout << GridLogMessage << "************************************************************************ "<<std::endl;
  result=zero;
  mCG(HermOp,src,result);
  std::cout << GridLogMessage << "************************************************************************ "<<std::endl;

  std::cout << GridLogMessage << "************************************************************************ "<<std::endl;
  std::cout << GridLogMessage << " Calling Block CG for "<<Ls <<" right hand sides" <<std::endl;
  std::cout << GridLogMessage << "************************************************************************ "<<std::endl;
  result=zero;
  BCG(HermOp,src,result);
  std::cout << GridLogMessage << "************************************************************************ "<<std::endl;


  Grid_finalize();
}
