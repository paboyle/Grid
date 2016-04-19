    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_dwf_even_odd.cc

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: paboyle <paboyle@ph.ed.ac.uk>

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

typedef WilsonFermion5D<DomainWallRedBlack5dImplR> WilsonFermion5DR;
typedef WilsonFermion5D<DomainWallRedBlack5dImplF> WilsonFermion5DF;
typedef WilsonFermion5D<DomainWallRedBlack5dImplD> WilsonFermion5DD;

typedef WilsonFermion5D<WilsonImplR> WilsonFermion5D_OKR;

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  int threads = GridThread::GetThreads();
  std::cout<<GridLogMessage << "Grid is setup to use "<<threads<<" threads"<<std::endl;


  const int Ls=32;

  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);

  GridCartesian         * sUGrid   = SpaceTimeGrid::makeFourDimDWFGrid(GridDefaultLatt(),GridDefaultMpi());
  GridRedBlackCartesian * sUrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(sUGrid);

  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimDWFGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimDWFRedBlackGrid(Ls,UGrid);

  std::vector<int> seeds4({1,2,3,4});
  std::vector<int> seeds5({5,6,7,8});

  GridParallelRNG          RNG4(UGrid);  RNG4.SeedFixedIntegers(seeds4);
  GridParallelRNG          RNG5(FGrid);  RNG5.SeedFixedIntegers(seeds5);

  LatticeFermion src   (FGrid); random(RNG5,src);
  LatticeFermion result(FGrid); result=zero;
  LatticeFermion    ref(FGrid);    ref=zero;
  LatticeFermion    tmp(FGrid);
  LatticeFermion    err(FGrid);

  LatticeGaugeField Umu(UGrid); random(RNG4,Umu);
  std::vector<LatticeColourMatrix> U(4,UGrid);

  // Only one non-zero (y)
  /*
  Umu=zero;
  for(int nn=0;nn<Nd;nn++){
    random(RNG4,U[nn]);
    PokeIndex<LorentzIndex>(Umu,U[nn],nn);
  }
  */

  RealD mass=0.1;
  RealD M5  =1.8;
  typename WilsonFermion5DR::ImplParams params;

  WilsonFermion5DR Dw(1,Umu,*FGrid,*FrbGrid,*sUGrid,*sUrbGrid,M5,params);

  Dw.Dhop(src,result,0);

  std::cout << "Norm src = "<<norm2(src)<<" Norm res = "<<norm2(result) << std::endl;



  GridCartesian         * FokGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FokrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);
  WilsonFermion5D_OKR Dok(Umu,*FokGrid,*FokrbGrid,*UGrid,*UrbGrid,M5,params);
  
  LatticeFermion src_ok   (FokGrid);
  LatticeFermion ref_ok(FokGrid); 
  LatticeFermion result_ok(FokGrid);
  
  
  for(int lidx=0;lidx<FGrid->lSites();lidx++){
    std::vector<int> lcoor;
    FGrid->LocalIndexToLocalCoor(lidx,lcoor);
    
    SpinColourVector siteSrc;

    peekLocalSite(siteSrc,src,lcoor);
    pokeLocalSite(siteSrc,src_ok,lcoor);

    peekLocalSite(siteSrc,result,lcoor);
    pokeLocalSite(siteSrc,result_ok,lcoor);
  }
  
  Dok.Dhop(src_ok,ref_ok,0);
  
  std::cout << "Reference = "<<norm2(src_ok)<<" res = "<<norm2(ref_ok) << std::endl;
  ref_ok = ref_ok - result_ok;
  std::cout << "Reference diff = "<<norm2(result_ok)<< std::endl;
  std::cout << "Reference diff = "<<norm2(ref_ok)<< std::endl;

  Grid_finalize();
}
