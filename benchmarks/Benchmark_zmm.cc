    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_zmm.cc

    Copyright (C) 2015

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
#include <PerfCount.h>


using namespace Grid;
using namespace Grid::QCD;


int bench(std::ofstream &os, std::vector<int> &latt4,int Ls);

int main(int argc,char **argv)
{
  Grid_init(&argc,&argv);
  std::ofstream os("zmm.dat");

  os << "#V Ls Lxy Lzt C++ Asm OMP L1 " <<std::endl;
  for(int L=4;L<=32;L+=4){
    for(int m=1;m<=2;m++){
      for(int Ls=8;Ls<=16;Ls+=8){
	std::vector<int> grid({L,L,m*L,m*L});
	for(int i=0;i<4;i++) { 
	  std::cout << grid[i]<<"x";
	}
	std::cout << Ls<<std::endl;
	bench(os,grid,Ls);
      }
    }
  }
}

int bench(std::ofstream &os, std::vector<int> &latt4,int Ls)
{

  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(latt4, GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  std::vector<int> simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
  std::vector<int> mpi_layout  = GridDefaultMpi();
  int threads = GridThread::GetThreads();

  std::vector<int> seeds4({1,2,3,4});
  std::vector<int> seeds5({5,6,7,8});

  GridSerialRNG sRNG; sRNG.SeedFixedIntegers(seeds4);

  LatticeFermion src (FGrid);
  LatticeFermion tmp (FGrid);
  LatticeFermion srce(FrbGrid);

  LatticeFermion resulto(FrbGrid); resulto=zero;
  LatticeFermion resulta(FrbGrid); resulta=zero;
  LatticeFermion junk(FrbGrid); junk=zero;
  LatticeFermion diff(FrbGrid); 
  LatticeGaugeField Umu(UGrid);

  double mfc, mfa, mfo, mfl1;

  GridParallelRNG          RNG4(UGrid);  RNG4.SeedFixedIntegers(seeds4);
  GridParallelRNG          RNG5(FGrid);  RNG5.SeedFixedIntegers(seeds5);
  random(RNG5,src);
#if 1
  random(RNG4,Umu);
#else
  int mmu=2;
  std::vector<LatticeColourMatrix> U(4,UGrid);
  for(int mu=0;mu<Nd;mu++){
    U[mu] = PeekIndex<LorentzIndex>(Umu,mu);
    if ( mu!=mmu ) U[mu] = zero;
    if ( mu==mmu ) U[mu] = 1.0;
    PokeIndex<LorentzIndex>(Umu,U[mu],mu);
  }
#endif
 pickCheckerboard(Even,srce,src);

  RealD mass=0.1;
  RealD M5  =1.8;
  DomainWallFermionR Dw(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);

  std::cout<<GridLogMessage << "Calling Dw"<<std::endl;
  int ncall=50;
  double t0=usecond();
  for(int i=0;i<ncall;i++){
    Dw.DhopOE(srce,resulto,0);
  }
  double t1=usecond();

  double volume=Ls;  for(int mu=0;mu<Nd;mu++) volume=volume*latt4[mu];
  double flops=1344*volume/2;

  mfc = flops*ncall/(t1-t0);
  std::cout<<GridLogMessage << "Called C++ Dw"<< " mflop/s =   "<< mfc<<std::endl;

  QCD::WilsonFermion5DStatic::AsmOptDslash=1;
  t0=usecond();
  for(int i=0;i<ncall;i++){
    Dw.DhopOE(srce,resulta,0);
  }
  t1=usecond();
  mfa = flops*ncall/(t1-t0);
  std::cout<<GridLogMessage << "Called ASM Dw"<< " mflop/s =   "<< mfa<<std::endl;
  /*
  int dag=DaggerNo;
  t0=usecond();
  for(int i=0;i<1;i++){
    Dw.DhopInternalOMPbench(Dw.StencilEven,Dw.LebesgueEvenOdd,Dw.UmuOdd,srce,resulta,dag);
  }
  t1=usecond();
  mfo = flops*100/(t1-t0);
  std::cout<<GridLogMessage << "Called ASM-OMP Dw"<< " mflop/s =   "<< mfo<<std::endl;

  t0=usecond();
  for(int i=0;i<1;i++){
    Dw.DhopInternalL1bench(Dw.StencilEven,Dw.LebesgueEvenOdd,Dw.UmuOdd,srce,resulta,dag);
  }
  t1=usecond();
  mfl1= flops*100/(t1-t0);
  std::cout<<GridLogMessage << "Called ASM-L1 Dw"<< " mflop/s =   "<< mfl1<<std::endl;
  os << latt4[0]*latt4[1]*latt4[2]*latt4[3]<< " "<<Ls<<" "<< latt4[0] <<" " <<latt4[2]<< " "
     << mfc<<" "
     << mfa<<" "
     << mfo<<" "
     << mfl1<<std::endl;
  */

#if 0
  for(int i=0;i< PerformanceCounter::NumTypes(); i++ ){
    Dw.DhopOE(srce,resulta,0);
    PerformanceCounter Counter(i);
    Counter.Start();
    Dw.DhopOE(srce,resulta,0);
    Counter.Stop();
    Counter.Report();
  }
#endif
  //resulta = (-0.5) * resulta;

  diff = resulto-resulta;
  std::cout<<GridLogMessage << "diff "<< norm2(diff)<<std::endl;
  std::cout<<std::endl;
  return 0;
}


