    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_dwf_fpgcr.cc

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
#include <Grid.h>
#include <algorithms/iterative/PrecGeneralisedConjugateResidual.h>
using namespace std;
using namespace Grid;
using namespace Grid::QCD;


template<class d>
struct scal {
  d internal;
};

  Gamma::GammaMatrix Gmu [] = {
    Gamma::GammaX,
    Gamma::GammaY,
    Gamma::GammaZ,
    Gamma::GammaT
  };

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);



  const int Ls=8;

  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);


  std::vector<int> seeds4({1,2,3,4});
  std::vector<int> seeds5({5,6,7,8});
  GridParallelRNG          RNG5(FGrid);  RNG5.SeedFixedIntegers(seeds5);
  GridParallelRNG          RNG4(UGrid);  RNG4.SeedFixedIntegers(seeds4);

  LatticeFermion    src(FGrid); random(RNG5,src);
  LatticeFermion result(FGrid); result=zero;
  LatticeGaugeField Umu(UGrid); 

  SU3::HotConfiguration(RNG4,Umu);

  TrivialPrecon<LatticeFermion> simple;

  PrecGeneralisedConjugateResidual<LatticeFermion> PGCR(1.0e-6,10000,simple,4,160);

  ConjugateResidual<LatticeFermion> CR(1.0e-6,10000);

  ConjugateGradient<LatticeFermion> CG(1.0e-6,10000);
  
  RealD mass=0.5;
  RealD M5=1.8;
  DomainWallFermionR Ddwf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);

  std::cout<<GridLogMessage<<"*********************************************************"<<std::endl;
  std::cout<<GridLogMessage<<"* Solving with MdagM VPGCR "<<std::endl;
  std::cout<<GridLogMessage<<"*********************************************************"<<std::endl;
  MdagMLinearOperator<DomainWallFermionR,LatticeFermion> HermOp(Ddwf);
  result=zero;
  PGCR(HermOp,src,result);

  std::cout<<GridLogMessage<<"*********************************************************"<<std::endl;
  std::cout<<GridLogMessage<<"* Solving with g5-VPGCR "<<std::endl;
  std::cout<<GridLogMessage<<"*********************************************************"<<std::endl;
  Gamma5R5HermitianLinearOperator<DomainWallFermionR,LatticeFermion> g5HermOp(Ddwf);
  result=zero;
  PGCR(g5HermOp,src,result);

  std::cout<<GridLogMessage<<"*********************************************************"<<std::endl;
  std::cout<<GridLogMessage<<"* Solving with MdagM-CR "<<std::endl;
  std::cout<<GridLogMessage<<"*********************************************************"<<std::endl;
  result=zero;
  CR(HermOp,src,result);

  std::cout<<GridLogMessage<<"*********************************************************"<<std::endl;
  std::cout<<GridLogMessage<<"* Solving with g5-CR "<<std::endl;
  std::cout<<GridLogMessage<<"*********************************************************"<<std::endl;
  result=zero;
  CR(g5HermOp,src,result);

  std::cout<<GridLogMessage<<"*********************************************************"<<std::endl;
  std::cout<<GridLogMessage<<"* Solving with MdagM-CG "<<std::endl;
  std::cout<<GridLogMessage<<"*********************************************************"<<std::endl;
  result=zero;
  CG(HermOp,src,result);

  Grid_finalize();
}
