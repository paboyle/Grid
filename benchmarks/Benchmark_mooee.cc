    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./benchmarks/Benchmark_dwf.cc

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
#include <Grid/Grid.h>

using namespace std;
using namespace Grid;



int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  int threads = GridThread::GetThreads();
  std::cout<<GridLogMessage << "Grid is setup to use "<<threads<<" threads"<<std::endl;

  Coordinate latt4 = GridDefaultLatt();
  const int Ls=16;
  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  std::vector<int> seeds4({1,2,3,4});
  std::vector<int> seeds5({5,6,7,8});

  GridParallelRNG          RNG4(UGrid);  RNG4.SeedFixedIntegers(seeds4);
  std::cout << GridLogMessage << "Seeded"<<std::endl;

  LatticeGaugeField Umu(UGrid); SU<Nc>::HotConfiguration(RNG4,Umu);

  std::cout << GridLogMessage << "made random gauge fields"<<std::endl;

  RealD mass=0.1;
  RealD M5  =1.8;
  //  RealD NP = UGrid->_Nprocessors;


  if (1)
  {
    const int ncall=1000;

    std::cout << GridLogMessage<< "*********************************************************" <<std::endl;
    std::cout << GridLogMessage<< "* Benchmarking DomainWallFermionD::Dhop "<<std::endl;
    std::cout << GridLogMessage<< "*********************************************************" <<std::endl;

    GridParallelRNG RNG5(FGrid); RNG5.SeedFixedIntegers(seeds5);
    LatticeFermion src(FGrid); random(RNG5,src);
    LatticeFermion result(FGrid);

    DomainWallFermionD Dw(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);
    double t0,t1;
    
    typedef typename DomainWallFermionD::Coeff_t Coeff_t;
    std::vector<Coeff_t> diag = Dw.bs;
    std::vector<Coeff_t> upper= Dw.cs;
    std::vector<Coeff_t> lower= Dw.cs;
    upper[Ls-1]=-Dw.mass_minus*upper[Ls-1];
    lower[0]   =-Dw.mass_plus*lower[0];
    
    LatticeFermion r_eo(FGrid);
    LatticeFermion src_e (FrbGrid);
    LatticeFermion src_o (FrbGrid);
    LatticeFermion r_e   (FrbGrid);
    LatticeFermion r_o   (FrbGrid);
    
    pickCheckerboard(Even,src_e,src);
    pickCheckerboard(Odd,src_o,src);
    
    setCheckerboard(r_eo,src_o);
    setCheckerboard(r_eo,src_e);
    
    r_e = Zero();
    r_o = Zero();


#define BENCH_DW(A,...)			\
    Dw. A (__VA_ARGS__);				\
    FGrid->Barrier();				\
    t0=usecond();				\
    for(int i=0;i<ncall;i++){			\
      Dw. A (__VA_ARGS__);				\
    }						\
    t1=usecond();				\
    FGrid->Barrier();				\
    std::cout<<GridLogMessage << "Called " #A " "<< (t1-t0)/ncall<<" us"<<std::endl;\
    std::cout<<GridLogMessage << "******************"<<std::endl;

#define BENCH_ZDW(A,in,out)			\
    zDw. A (in,out);				\
    FGrid->Barrier();				\
    t0=usecond();				\
    for(int i=0;i<ncall;i++){			\
      zDw. A (in,out);				\
    }						\
    t1=usecond();				\
    FGrid->Barrier();				\
    std::cout<<GridLogMessage << "Called ZDw " #A " "<< (t1-t0)/ncall<<" us"<<std::endl;\
    std::cout<<GridLogMessage << "******************"<<std::endl;

#define BENCH_DW_SSC(A,in,out)			\
    Dw. A (in,out);				\
    FGrid->Barrier();				\
    t0=usecond();				\
    for(int i=0;i<ncall;i++){			\
      __SSC_START ;				\
      Dw. A (in,out);				\
      __SSC_STOP ;				\
    }						\
    t1=usecond();				\
    FGrid->Barrier();				\
    std::cout<<GridLogMessage << "Called " #A " "<< (t1-t0)/ncall<<" us"<<std::endl;\
    std::cout<<GridLogMessage << "******************"<<std::endl;

    BENCH_DW(Dhop    ,src,result,0);
    BENCH_DW(DhopEO  ,src_o,r_e,0);
    BENCH_DW(Meooe   ,src_o,r_e);
    BENCH_DW(M5D     ,src_o,src_o,r_e,lower,diag,upper);
    BENCH_DW(Mooee   ,src_o,r_o);
    BENCH_DW(MooeeInv,src_o,r_o);

  }

  Grid_finalize();
}
