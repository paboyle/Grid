 /*************************************************************************************
    Grid physics library, www.github.com/paboyle/Grid
    Source file: ./benchmarks/Benchmark_gpdwf_Xconj.cc
    Copyright (C) 2015

    Author: Christopher Kelly <ckelly@bnl.gov>
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
#include <sstream>
using namespace std;
using namespace Grid;

typedef typename XconjugateDomainWallFermionF::FermionField LatticeFermionF;
typedef typename XconjugateDomainWallFermionD::FermionField LatticeFermionD;

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);
#ifdef ENABLE_GPARITY
  int Ls=16;
  bool do_fp32 = true;
  bool do_fp64 = false;

  for(int i=0;i<argc;i++){
    if(std::string(argv[i]) == "-Ls"){
      std::stringstream ss(argv[i+1]); ss >> Ls;
    }else if(std::string(argv[i]) == "-precision"){
      std::string p = argv[i+1];
      if(p == "both"){ do_fp32 = true; do_fp64 = true; }
      else if(p == "single"){ do_fp32 = true; do_fp64 = false; } //default
      else if(p == "double"){ do_fp32 = false; do_fp64 = true; }
      else{
	assert(0 && "Invalid precision argument");
      }
    }
  }

  int threads = GridThread::GetThreads();
  std::cout<<GridLogMessage << "Grid is setup to use "<<threads<<" threads"<<std::endl;
  std::cout<<GridLogMessage << "Ls = " << Ls << std::endl;

  Coordinate latt4 = GridDefaultLatt();

  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplexF::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  std::vector<int> seeds4({1,2,3,4});
  std::vector<int> seeds5({5,6,7,8});
  
  std::cout << GridLogMessage << "Initialising 4d RNG" << std::endl;
  GridParallelRNG          RNG4(UGrid);  RNG4.SeedFixedIntegers(seeds4);
  std::cout << GridLogMessage << "Initialising 5d RNG" << std::endl;
  GridParallelRNG          RNG5(FGrid);  RNG5.SeedFixedIntegers(seeds5);
  std::cout << GridLogMessage << "Initialised RNGs" << std::endl;

  LatticeFermionF src   (FGrid); random(RNG5,src);
  RealD N2 = 1.0/::sqrt(norm2(src));
  src = src*N2;

  LatticeFermionF result(FGrid); result=Zero();
  LatticeFermionF    ref(FGrid);    ref=Zero();
  LatticeFermionF    tmp(FGrid);
  LatticeFermionF    err(FGrid);

  std::cout << GridLogMessage << "Drawing gauge field" << std::endl;
  LatticeGaugeFieldF Umu(UGrid); 
  SU<Nc>::HotConfiguration(RNG4,Umu); 
  std::cout << GridLogMessage << "Random gauge initialised " << std::endl;

  RealD mass=0.1;
  RealD M5  =1.8;

  RealD NP = UGrid->_Nprocessors;
  RealD NN = UGrid->NodeCount();

  std::cout << GridLogMessage<< "*****************************************************************" <<std::endl;
  std::cout << GridLogMessage<< "* Kernel options --dslash-generic, --dslash-unroll, --dslash-asm" <<std::endl;
  std::cout << GridLogMessage<< "*****************************************************************" <<std::endl;
  std::cout << GridLogMessage<< "*****************************************************************" <<std::endl;
  std::cout << GridLogMessage<< "* Benchmarking XconjugateDomainWallFermion::Dhop                  "<<std::endl;
  std::cout << GridLogMessage<< "* Vectorising space-time by "<<vComplexF::Nsimd()<<std::endl;
#ifdef GRID_OMP
  if ( WilsonKernelsStatic::Comms == WilsonKernelsStatic::CommsAndCompute ) std::cout << GridLogMessage<< "* Using Overlapped Comms/Compute" <<std::endl;
  if ( WilsonKernelsStatic::Comms == WilsonKernelsStatic::CommsThenCompute) std::cout << GridLogMessage<< "* Using sequential comms compute" <<std::endl;
#endif
  if ( WilsonKernelsStatic::Opt == WilsonKernelsStatic::OptGeneric   ) std::cout << GridLogMessage<< "* Using GENERIC Nc WilsonKernels" <<std::endl;
  if ( WilsonKernelsStatic::Opt == WilsonKernelsStatic::OptHandUnroll) std::cout << GridLogMessage<< "* Using Nc=3       WilsonKernels" <<std::endl;
  if ( WilsonKernelsStatic::Opt == WilsonKernelsStatic::OptInlineAsm ) std::cout << GridLogMessage<< "* Using Asm Nc=3   WilsonKernels" <<std::endl;
  std::cout << GridLogMessage<< "*****************************************************************" <<std::endl;

  std::vector<int> twists({1,1,1,0});
  XconjugateDomainWallFermionF::ImplParams xparams;
  xparams.twists = twists;
  xparams.boundary_phase = 1.0;

  int ncall =1000;
  
  if(do_fp32){
    std::cout << GridLogMessage<< "* SINGLE/SINGLE"<<std::endl;
    XconjugateDomainWallFermionF Dw(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,xparams);
    FGrid->Barrier();
    Dw.Dhop(src,result,0);
    std::cout<<GridLogMessage<<"Called warmup"<<std::endl;
    double t0=usecond();
    for(int i=0;i<ncall;i++){
      Dw.Dhop(src,result,0);
    }
    double t1=usecond();
    FGrid->Barrier();
    
    double volume=Ls;  for(int mu=0;mu<Nd;mu++) volume=volume*latt4[mu];
    double flops=1320*volume*ncall;

    std::cout<<GridLogMessage << "Called Dw "<<ncall<<" times in "<<t1-t0<<" us"<<std::endl;
    std::cout<<GridLogMessage << "mflop/s =   "<< flops/(t1-t0)<<std::endl;
    std::cout<<GridLogMessage << "mflop/s per rank =  "<< flops/(t1-t0)/NP<<std::endl;
    std::cout<<GridLogMessage << "mflop/s per node =  "<< flops/(t1-t0)/NN<<std::endl;
  }

  if(do_fp64){
    GridCartesian         * UGrid_d   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplexD::Nsimd()),GridDefaultMpi());
    GridRedBlackCartesian * UrbGrid_d = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid_d);
    GridCartesian         * FGrid_d   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid_d);
    GridRedBlackCartesian * FrbGrid_d = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid_d);

    std::cout << GridLogMessage<< "* DOUBLE/DOUBLE"<<std::endl;
    LatticeFermionD src_d(FGrid_d);
    precisionChange(src_d,src);
    
    LatticeGaugeFieldD Umu_d(UGrid_d); 
    precisionChange(Umu_d,Umu);
    
    LatticeFermionD result_d(FGrid_d);

    XconjugateDomainWallFermionD DwD(Umu_d,*FGrid_d,*FrbGrid_d,*UGrid_d,*UrbGrid_d,mass,M5,xparams);
    FGrid_d->Barrier();
    DwD.Dhop(src_d,result_d,0);
    std::cout<<GridLogMessage<<"Called warmup"<<std::endl;
    double t0=usecond();
    for(int i=0;i<ncall;i++){
      DwD.Dhop(src_d,result_d,0);
    }
    double t1=usecond();
    FGrid_d->Barrier();
      
    double volume=Ls;  for(int mu=0;mu<Nd;mu++) volume=volume*latt4[mu];
    double flops=1320*volume*ncall;
      
    std::cout<<GridLogMessage << "Called Dw "<<ncall<<" times in "<<t1-t0<<" us"<<std::endl;
    std::cout<<GridLogMessage << "mflop/s =   "<< flops/(t1-t0)<<std::endl;
    std::cout<<GridLogMessage << "mflop/s per rank =  "<< flops/(t1-t0)/NP<<std::endl;
    std::cout<<GridLogMessage << "mflop/s per node =  "<< flops/(t1-t0)/NN<<std::endl;
  }

#endif
  Grid_finalize();
}
