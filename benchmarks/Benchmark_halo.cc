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
#ifdef GRID_CUDA
#define CUDA_PROFILE
#endif

#ifdef CUDA_PROFILE
#include <cuda_profiler_api.h>
#endif

using namespace std;
using namespace Grid;

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

  Coordinate latt4= GridDefaultLatt();
  Coordinate mpi  = GridDefaultMpi();
  Coordinate simd = GridDefaultSimd(Nd,vComplexF::Nsimd());

  GridLogLayout();

  int Ls=16;
  for(int i=0;i<argc;i++)
    if(std::string(argv[i]) == "-Ls"){
      std::stringstream ss(argv[i+1]); ss >> Ls;
    }

  
  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(latt4,simd ,mpi);
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  std::cout << GridLogMessage << "Making s innermost grids"<<std::endl;
  GridCartesian         * sUGrid   = SpaceTimeGrid::makeFourDimDWFGrid(GridDefaultLatt(),GridDefaultMpi());
  GridRedBlackCartesian * sUrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(sUGrid);
  GridCartesian         * sFGrid   = SpaceTimeGrid::makeFiveDimDWFGrid(Ls,UGrid);
  GridRedBlackCartesian * sFrbGrid = SpaceTimeGrid::makeFiveDimDWFRedBlackGrid(Ls,UGrid);

  std::vector<int> seeds4({1,2,3,4});
  std::vector<int> seeds5({5,6,7,8});

  std::cout << GridLogMessage << "Initialising 4d RNG" << std::endl;
  GridParallelRNG          RNG4(UGrid);  RNG4.SeedUniqueString(std::string("The 4D RNG"));
  std::cout << GridLogMessage << "Initialising 5d RNG" << std::endl;
  GridParallelRNG          RNG5(FGrid);  RNG5.SeedUniqueString(std::string("The 5D RNG"));
  std::cout << GridLogMessage << "Initialised RNGs" << std::endl;

  LatticeFermionF src   (FGrid); random(RNG5,src);
  RealD N2 = 1.0/::sqrt(norm2(src));
  src = src*N2;

  std::cout << GridLogMessage << "Drawing gauge field" << std::endl;
  LatticeGaugeFieldF Umu(UGrid);
  SU<Nc>::HotConfiguration(RNG4,Umu);
  std::cout << GridLogMessage << "Random gauge initialised " << std::endl;

  RealD mass=0.1;
  RealD M5  =1.8;

  RealD NP = UGrid->_Nprocessors;
  RealD NN = UGrid->NodeCount();

  DomainWallFermionF Dw(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);

  const int ncall = 500;
  std::cout << GridLogMessage<< "*********************************************************" <<std::endl;
  std::cout << GridLogMessage<< "* Benchmarking DomainWallFermionF::HaloGatherOpt         "<<std::endl;
  std::cout << GridLogMessage<< "*********************************************************" <<std::endl;
  {
    typename DomainWallFermionF::Compressor compressor(0);
    FGrid->Barrier();
    Dw.Stencil.HaloExchangeOptGather(src,compressor);
    double t0=usecond();
    for(int i=0;i<ncall;i++){
      Dw.Stencil.HaloExchangeOptGather(src,compressor);
    }
    double t1=usecond();
    FGrid->Barrier();

    double bytes=0.0;
    if(mpi[0]) bytes+=latt4[1]*latt4[2]*latt4[3];
    if(mpi[1]) bytes+=latt4[0]*latt4[2]*latt4[3];
    if(mpi[2]) bytes+=latt4[0]*latt4[1]*latt4[3];
    if(mpi[3]) bytes+=latt4[0]*latt4[1]*latt4[2];
    bytes = bytes * Ls * 8.* (24.+12.)* 2.0;

    std::cout<<GridLogMessage << "Gather us /call =   "<< (t1-t0)/ncall<<std::endl;
    std::cout<<GridLogMessage << "Gather MBs /call =   "<< bytes*ncall/(t1-t0)<<std::endl;

  }

  Grid_finalize();
  exit(0);
}
