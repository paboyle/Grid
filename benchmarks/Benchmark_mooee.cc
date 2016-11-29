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
using namespace Grid::QCD;


int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  int threads = GridThread::GetThreads();
  std::cout<<GridLogMessage << "Grid is setup to use "<<threads<<" threads"<<std::endl;

  std::vector<int> latt4 = GridDefaultLatt();
  const int Ls=8;
  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  std::cout << GridLogMessage << "Making Vec5d innermost grids"<<std::endl;
  GridCartesian         * sUGrid   = SpaceTimeGrid::makeFourDimDWFGrid(GridDefaultLatt(),GridDefaultMpi());
  GridRedBlackCartesian * sUrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(sUGrid);
  GridCartesian         * sFGrid   = SpaceTimeGrid::makeFiveDimDWFGrid(Ls,UGrid);
  GridRedBlackCartesian * sFrbGrid = SpaceTimeGrid::makeFiveDimDWFRedBlackGrid(Ls,UGrid);

  std::vector<int> seeds4({1,2,3,4});
  std::vector<int> seeds5({5,6,7,8});

  GridParallelRNG          RNG4(UGrid);  RNG4.SeedFixedIntegers(seeds4);
  std::cout << GridLogMessage << "Seeded"<<std::endl;

  LatticeGaugeField Umu(UGrid); SU3::HotConfiguration(RNG4,Umu);

  std::cout << GridLogMessage << "made random gauge fields"<<std::endl;

  RealD mass=0.1;
  RealD M5  =1.8;
  RealD NP = UGrid->_Nprocessors;


  if (1)
  {
    const int ncall=100;

    std::cout << GridLogMessage<< "*********************************************************" <<std::endl;
    std::cout << GridLogMessage<< "* Benchmarking DomainWallFermionR::Dhop "<<std::endl;
    std::cout << GridLogMessage<< "*********************************************************" <<std::endl;

    GridParallelRNG RNG5(FGrid);
    LatticeFermion src(FGrid); random(RNG5,src);
    LatticeFermion result(FGrid);

    DomainWallFermionR Dw(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);

    FGrid->Barrier();

    double t0,t1;
    t0=usecond();
    for(int i=0;i<ncall;i++){
      Dw.Dhop(src,result,0);
    }
    t1=usecond();
    FGrid->Barrier();

    std::cout<<GridLogMessage << "Called Dhop "<< (t1-t0)/ncall<<" us"<<std::endl;

    LatticeFermion r_eo(FGrid);
    LatticeFermion src_e (FrbGrid);
    LatticeFermion src_o (FrbGrid);
    LatticeFermion r_e   (FrbGrid);
    LatticeFermion r_o   (FrbGrid);
    
    pickCheckerboard(Even,src_e,src);
    pickCheckerboard(Odd,src_o,src);
    
    setCheckerboard(r_eo,src_o);
    setCheckerboard(r_eo,src_e);
    
    r_e = zero;
    r_o = zero;

    FGrid->Barrier();
    t0=usecond();
    for (int i = 0; i < ncall; i++) {
      Dw.DhopEO(src_o, r_e, DaggerNo);
    }
    t1=usecond();
    FGrid->Barrier();
    std::cout<<GridLogMessage << "Called DhopEO "<< (t1-t0)/ncall<<" us"<<std::endl;

    FGrid->Barrier();
    t0=usecond();
    for (int i = 0; i < ncall; i++) {
      Dw.Mooee(src_o, r_o);
    }
    t1=usecond();
    FGrid->Barrier();
    std::cout<<GridLogMessage << "Called Mooee "<< (t1-t0)/ncall<<" us"<<std::endl;

    FGrid->Barrier();
    t0=usecond();
    for (int i = 0; i < ncall; i++) {
      Dw.MooeeInv(src_o, r_o);
    }
    t1=usecond();
    FGrid->Barrier();
    std::cout<<GridLogMessage << "Called MooeeInv "<< (t1-t0)/ncall<<" us"<<std::endl;


    FGrid->Barrier();
    t0=usecond();
    for (int i = 0; i < ncall; i++) {
      Dw.Meooe(src_o, r_e);
    }
    t1=usecond();
    FGrid->Barrier();
    std::cout<<GridLogMessage << "Called Meooe "<< (t1-t0)/ncall<<" us"<<std::endl;

  }

  if (1)
  {
    const int ncall=100;

    std::cout << GridLogMessage<< "*********************************************************" <<std::endl;
    std::cout << GridLogMessage<< "* Benchmarking DomainWallFermionVec5dR::Dhop "<<std::endl;
    std::cout << GridLogMessage<< "*********************************************************" <<std::endl;

    GridParallelRNG RNG5(sFGrid);
    LatticeFermion src(sFGrid); random(RNG5,src);
    LatticeFermion sref(sFGrid);
    LatticeFermion result(sFGrid);

    std::cout<<GridLogMessage << "Constructing Vec5D Dw "<<std::endl;
    DomainWallFermionVec5dR Dw(Umu,*sFGrid,*sFrbGrid,*sUGrid,*sUrbGrid,mass,M5);

    std::cout<<GridLogMessage << "Calling Dhop "<<std::endl;
    FGrid->Barrier();

    double t0,t1;
    t0=usecond();
    for(int i=0;i<ncall;i++){
      Dw.Dhop(src,result,0);
    }
    t1=usecond();
    FGrid->Barrier();

    std::cout<<GridLogMessage << "Called Vec5D Dhop "<< (t1-t0)/ncall<<" us"<<std::endl;

    LatticeFermion r_eo(sFGrid);
    LatticeFermion src_e (sFrbGrid);
    LatticeFermion src_o (sFrbGrid);
    LatticeFermion r_e   (sFrbGrid);
    LatticeFermion r_o   (sFrbGrid);
    
    pickCheckerboard(Even,src_e,src);
    pickCheckerboard(Odd,src_o,src);
    
    setCheckerboard(r_eo,src_o);
    setCheckerboard(r_eo,src_e);
    
    r_e = zero;
    r_o = zero;

    FGrid->Barrier();
    t0=usecond();
    for (int i = 0; i < ncall; i++) {
      Dw.DhopEO(src_o, r_e, DaggerNo);
    }
    t1=usecond();
    FGrid->Barrier();
    std::cout<<GridLogMessage << "Called Vec5D DhopEO "<< (t1-t0)/ncall<<" us"<<std::endl;

    FGrid->Barrier();
    t0=usecond();
    for (int i = 0; i < ncall; i++) {
      Dw.Mooee(src_o, r_o);
    }
    t1=usecond();
    FGrid->Barrier();
    std::cout<<GridLogMessage << "Called Vec5D Mooee "<< (t1-t0)/ncall<<" us"<<std::endl;

    FGrid->Barrier();
    t0=usecond();
    for (int i = 0; i < ncall; i++) {
      Dw.MooeeInv(src_o, r_o);
    }
    t1=usecond();
    FGrid->Barrier();
    std::cout<<GridLogMessage << "Called Vec5D MooeeInv "<< (t1-t0)/ncall<<" us"<<std::endl;


    FGrid->Barrier();
    t0=usecond();
    for (int i = 0; i < ncall; i++) {
      Dw.Meooe(src_o, r_e);
    }
    t1=usecond();
    FGrid->Barrier();
    std::cout<<GridLogMessage << "Called Vec5D Meooe "<< (t1-t0)/ncall<<" us"<<std::endl;

  }



  Grid_finalize();
}
