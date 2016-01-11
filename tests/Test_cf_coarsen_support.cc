    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_cf_coarsen_support.cc

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
#include <Grid.h>

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

  const int Ls=9;

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
  LatticeFermion    ref(FGrid); ref=zero;
  LatticeFermion    tmp(FGrid);
  LatticeFermion    err(FGrid);
  LatticeGaugeField Umu(UGrid); random(RNG4,Umu);

  std::vector<LatticeColourMatrix> U(4,UGrid);
  for(int mu=0;mu<Nd;mu++){
    U[mu] = PeekIndex<LorentzIndex>(Umu,mu);
  }
  
  RealD mass=0.1;
  RealD M5=1.8;

  {
    OverlapWilsonContFracTanhFermionR Dcf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,1.0);
    HermitianLinearOperator<OverlapWilsonContFracTanhFermionR,LatticeFermion> HermIndefOp(Dcf);

    HermIndefOp.Op(src,ref);
    HermIndefOp.OpDiag(src,result);
    
    for(int d=0;d<4;d++){
      HermIndefOp.OpDir(src,tmp,d,+1); result=result+tmp; 
      std::cout<<GridLogMessage<<"dir "<<d<<" tmp "<<norm2(tmp)<<std::endl;
      HermIndefOp.OpDir(src,tmp,d,-1); result=result+tmp;
      std::cout<<GridLogMessage<<"dir "<<d<<" tmp "<<norm2(tmp)<<std::endl;
    }
    err = result-ref;
    std::cout<<GridLogMessage<<"Error "<<norm2(err)<<std::endl;
  }

  {
    OverlapWilsonPartialFractionTanhFermionR Dpf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,1.0);
    HermitianLinearOperator<OverlapWilsonPartialFractionTanhFermionR,LatticeFermion> HermIndefOp(Dpf);
    
    HermIndefOp.Op(src,ref);
    HermIndefOp.OpDiag(src,result);
    
    for(int d=0;d<4;d++){
      HermIndefOp.OpDir(src,tmp,d,+1); result=result+tmp; 
      std::cout<<GridLogMessage<<"dir "<<d<<" tmp "<<norm2(tmp)<<std::endl;
      HermIndefOp.OpDir(src,tmp,d,-1); result=result+tmp;
      std::cout<<GridLogMessage<<"dir "<<d<<" tmp "<<norm2(tmp)<<std::endl;
    }

    err = result-ref;
    std::cout<<GridLogMessage<<"Error "<<norm2(err)<<std::endl;
  }


  Grid_finalize();
}
