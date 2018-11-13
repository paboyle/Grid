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

  const int Ls=8;

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

  FermionField src(FGrid);
  FermionField tt(FGrid);
#if 1
  random(pRNG5,src);
#else
  src=zero;
  ComplexField coor(FGrid);
  LatticeCoordinate(coor,0);
  for(int ss=0;ss<FGrid->oSites();ss++){
    src._odata[ss]()()(0)=coor._odata[ss]()()();
  }
  LatticeCoordinate(coor,1);
  for(int ss=0;ss<FGrid->oSites();ss++){
    src._odata[ss]()()(0)+=coor._odata[ss]()()();
  }
#endif
  FermionField src_o(FrbGrid);   pickCheckerboard(Odd,src_o,src);
  FermionField result_o(FrbGrid); result_o=zero; 
  RealD nrm = norm2(src);

  LatticeGaugeField Umu(UGrid); SU3::HotConfiguration(pRNG,Umu);

  double volume=1;
  for(int mu=0;mu<Nd;mu++){
    volume=volume*latt_size[mu];
  }  

  RealD mass=0.003;
  RealD c1=9.0/8.0;
  RealD c2=-1.0/24.0;
  RealD u0=1.0;
  ImprovedStaggeredFermion5DR Ds(Umu,Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,c1,c2,u0); 
  SchurStaggeredOperator<ImprovedStaggeredFermion5DR,FermionField> HermOp(Ds);

  ConjugateGradient<FermionField> CG(1.0e-8,10000);
  int blockDim = 0;
  BlockConjugateGradient<FermionField>    BCGrQ(BlockCGrQ,blockDim,1.0e-8,10000);
  BlockConjugateGradient<FermionField>    BCG  (BlockCGrQ,blockDim,1.0e-8,10000);
  BlockConjugateGradient<FermionField>    BCGv (BlockCGrQVec,blockDim,1.0e-8,10000);
  BlockConjugateGradient<FermionField>    mCG  (CGmultiRHS,blockDim,1.0e-8,10000);

  std::cout << GridLogMessage << "****************************************************************** "<<std::endl;
  std::cout << GridLogMessage << " Calling 4d CG "<<std::endl;
  std::cout << GridLogMessage << "****************************************************************** "<<std::endl;
  ImprovedStaggeredFermionR Ds4d(Umu,Umu,*UGrid,*UrbGrid,mass,c1,c2,u0);
  SchurStaggeredOperator<ImprovedStaggeredFermionR,FermionField> HermOp4d(Ds4d);
  FermionField src4d(UGrid); random(pRNG,src4d);
  FermionField src4d_o(UrbGrid);   pickCheckerboard(Odd,src4d_o,src4d);
  FermionField result4d_o(UrbGrid); 

  double deodoe_flops=(16*(3*(6+8+8)) + 15*3*2)*volume; // == 66*16 +  == 1146
  result4d_o=zero;
  {
    double t1=usecond();
    CG(HermOp4d,src4d_o,result4d_o);
    double t2=usecond();
    double ncall=CG.IterationsToComplete;
    double flops = deodoe_flops * ncall;
    std::cout<<GridLogMessage << "usec    =   "<< (t2-t1)<<std::endl;
    std::cout<<GridLogMessage << "flops   =   "<< flops<<std::endl;
    std::cout<<GridLogMessage << "mflop/s =   "<< flops/(t2-t1)<<std::endl;
    HermOp4d.Report();
  }
  Ds4d.Report();
  std::cout << GridLogMessage << "************************************************************************ "<<std::endl;


  std::cout << GridLogMessage << "************************************************************************ "<<std::endl;
  std::cout << GridLogMessage << " Calling 5d CG for "<<Ls <<" right hand sides" <<std::endl;
  std::cout << GridLogMessage << "************************************************************************ "<<std::endl;
  Ds.ZeroCounters();
  result_o=zero;
  {
    double t1=usecond();
    CG(HermOp,src_o,result_o);
    double t2=usecond();
    double ncall=CG.IterationsToComplete*Ls;
    double flops = deodoe_flops * ncall;
    std::cout<<GridLogMessage << "usec    =   "<< (t2-t1)<<std::endl;
    std::cout<<GridLogMessage << "flops   =   "<< flops<<std::endl;
    std::cout<<GridLogMessage << "mflop/s =   "<< flops/(t2-t1)<<std::endl;
    HermOp.Report();
  }
  Ds.Report();
  std::cout << GridLogMessage << "************************************************************************ "<<std::endl;

  std::cout << GridLogMessage << "************************************************************************ "<<std::endl;
  std::cout << GridLogMessage << " Calling multiRHS CG for "<<Ls <<" right hand sides" <<std::endl;
  std::cout << GridLogMessage << "************************************************************************ "<<std::endl;
  Ds.ZeroCounters();
  result_o=zero;
  {
    double t1=usecond();
    mCG(HermOp,src_o,result_o);
    double t2=usecond();
    double ncall=mCG.IterationsToComplete*Ls;
    double flops = deodoe_flops * ncall;
    std::cout<<GridLogMessage << "usec    =   "<< (t2-t1)<<std::endl;
    std::cout<<GridLogMessage << "flops   =   "<< flops<<std::endl;
    std::cout<<GridLogMessage << "mflop/s =   "<< flops/(t2-t1)<<std::endl;
    HermOp.Report();
  }

  Ds.Report();
  std::cout << GridLogMessage << "************************************************************************ "<<std::endl;

  std::cout << GridLogMessage << "************************************************************************ "<<std::endl;
  std::cout << GridLogMessage << " Calling Block CGrQ for "<<Ls <<" right hand sides" <<std::endl;
  std::cout << GridLogMessage << "************************************************************************ "<<std::endl;
  Ds.ZeroCounters();
  result_o=zero;
  {
    double t1=usecond();
    BCGrQ(HermOp,src_o,result_o);
    double t2=usecond();
    double ncall=BCGrQ.IterationsToComplete*Ls;
    double flops = deodoe_flops * ncall;
    std::cout<<GridLogMessage << "usec    =   "<< (t2-t1)<<std::endl;
    std::cout<<GridLogMessage << "flops   =   "<< flops<<std::endl;
    std::cout<<GridLogMessage << "mflop/s =   "<< flops/(t2-t1)<<std::endl;
    HermOp.Report();
  }
  Ds.Report();
  std::cout << GridLogMessage << "************************************************************************ "<<std::endl;

  std::cout << GridLogMessage << "************************************************************************ "<<std::endl;
  std::cout << GridLogMessage << " Calling Block CG for "<<Ls <<" right hand sides" <<std::endl;
  std::cout << GridLogMessage << "************************************************************************ "<<std::endl;
  Ds.ZeroCounters();
  result_o=zero;
  {
    double t1=usecond();
    BCG(HermOp,src_o,result_o);
    double t2=usecond();
    double ncall=BCGrQ.IterationsToComplete*Ls;
    double flops = deodoe_flops * ncall;
    std::cout<<GridLogMessage << "usec    =   "<< (t2-t1)<<std::endl;
    std::cout<<GridLogMessage << "flops   =   "<< flops<<std::endl;
    std::cout<<GridLogMessage << "mflop/s =   "<< flops/(t2-t1)<<std::endl;
    HermOp.Report();
  }
  Ds.Report();
  std::cout << GridLogMessage << "************************************************************************ "<<std::endl;

  std::cout << GridLogMessage << "****************************************************************** "<<std::endl;
  std::cout << GridLogMessage << " Calling BCGvec "<<std::endl;
  std::cout << GridLogMessage << "****************************************************************** "<<std::endl;
  std::vector<FermionField> src_v   (Ls,UrbGrid);
  std::vector<FermionField> result_v(Ls,UrbGrid);
  for(int s=0;s<Ls;s++) result_v[s] = zero;
  for(int s=0;s<Ls;s++) {
    FermionField src4(UGrid);
    ExtractSlice(src4,src,s,0);
    pickCheckerboard(Odd,src_v[s],src4);  
  }

  {
    double t1=usecond();
    BCGv(HermOp4d,src_v,result_v);
    double t2=usecond();
    double ncall=BCGv.IterationsToComplete*Ls;
    double flops = deodoe_flops * ncall;
    std::cout<<GridLogMessage << "usec    =   "<< (t2-t1)<<std::endl;
    std::cout<<GridLogMessage << "flops   =   "<< flops<<std::endl;
    std::cout<<GridLogMessage << "mflop/s =   "<< flops/(t2-t1)<<std::endl;
    //    HermOp4d.Report();
  }


  Grid_finalize();
}
