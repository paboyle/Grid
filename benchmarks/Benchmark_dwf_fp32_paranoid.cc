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


  int threads = GridThread::GetThreads();

  Coordinate latt4 = GridDefaultLatt();
  int Ls=16;
  for(int i=0;i<argc;i++)
    if(std::string(argv[i]) == "-Ls"){
      std::stringstream ss(argv[i+1]); ss >> Ls;
    }

  GridLogLayout();

  long unsigned int single_site_flops = 8*Nc*(7+16*Nc);


  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplexF::Nsimd()),GridDefaultMpi());
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
  LatticeFermionF src1   (FGrid); random(RNG5,src1);
#if 0
  src = Zero();
  {
    Coordinate origin({0,0,0,latt4[2]-1,0});
    SpinColourVectorF tmp;
    tmp=Zero();
    tmp()(0)(0)=Complex(-2.0,0.0);
    std::cout << " source site 0 " << tmp<<std::endl;
    pokeSite(tmp,src,origin);
  }
#else
  RealD N2 = 1.0/::sqrt(norm2(src));
  src = src*N2;
#endif


  LatticeFermionF result(FGrid); result=Zero();
  LatticeFermionF    ref(FGrid);    ref=Zero();
  LatticeFermionF    tmp(FGrid);
  LatticeFermionF    err(FGrid);

  std::cout << GridLogMessage << "Drawing gauge field" << std::endl;
  LatticeGaugeFieldF Umu(UGrid);
  SU<Nc>::HotConfiguration(RNG4,Umu);
  std::cout << GridLogMessage << "Random gauge initialised " << std::endl;
#if 0
  Umu=1.0;
  for(int mu=0;mu<Nd;mu++){
    LatticeColourMatrixF ttmp(UGrid);
    ttmp = PeekIndex<LorentzIndex>(Umu,mu);
    //    if (mu !=2 ) ttmp = 0;
    //    ttmp = ttmp* pow(10.0,mu);
    PokeIndex<LorentzIndex>(Umu,ttmp,mu);
  }
  std::cout << GridLogMessage << "Forced to diagonal " << std::endl;
#endif

  ////////////////////////////////////
  // Naive wilson implementation
  ////////////////////////////////////
  // replicate across fifth dimension
  //  LatticeGaugeFieldF Umu5d(FGrid);
  std::vector<LatticeColourMatrixF> U(4,UGrid);
  for(int mu=0;mu<Nd;mu++){
    U[mu] = PeekIndex<LorentzIndex>(Umu,mu);
  }
  std::cout << GridLogMessage << "Setting up Cshift based reference " << std::endl;

  if (1)
  {
    ref = Zero();
    for(int mu=0;mu<Nd;mu++){

      tmp = Cshift(src,mu+1,1);
      {
	autoView( tmp_v  , tmp  , CpuWrite);
	autoView( U_v  , U[mu]  , CpuRead);
	for(int ss=0;ss<U[mu].Grid()->oSites();ss++){
	  for(int s=0;s<Ls;s++){
	    tmp_v[Ls*ss+s] = U_v[ss]*tmp_v[Ls*ss+s];
	  }
	}
      }
      ref=ref + tmp - Gamma(Gmu[mu])*tmp;

      {
	autoView( tmp_v  , tmp  , CpuWrite);
	autoView( U_v  , U[mu]  , CpuRead);
	autoView( src_v, src    , CpuRead);
	for(int ss=0;ss<U[mu].Grid()->oSites();ss++){
	  for(int s=0;s<Ls;s++){
	    tmp_v[Ls*ss+s] = adj(U_v[ss])*src_v[Ls*ss+s];
	  }
	}
      }
      tmp =Cshift(tmp,mu+1,-1);
      ref=ref + tmp + Gamma(Gmu[mu])*tmp;
    }
    ref = -0.5*ref;
  }

  RealD mass=0.1;
  RealD M5  =1.8;

  RealD NP = UGrid->_Nprocessors;
  RealD NN = UGrid->NodeCount();

  std::cout << GridLogMessage<< "*****************************************************************" <<std::endl;
  std::cout << GridLogMessage<< "* Kernel options --dslash-generic, --dslash-unroll, --dslash-asm" <<std::endl;
  std::cout << GridLogMessage<< "*****************************************************************" <<std::endl;
  std::cout << GridLogMessage<< "*****************************************************************" <<std::endl;
  std::cout << GridLogMessage<< "* Benchmarking DomainWallFermionR::Dhop                  "<<std::endl;
  std::cout << GridLogMessage<< "* Vectorising space-time by "<<vComplexF::Nsimd()<<std::endl;
  std::cout << GridLogMessage<< "* VComplexF size is "<<sizeof(vComplexF)<< " B"<<std::endl;
  if ( sizeof(RealF)==4 )   std::cout << GridLogMessage<< "* SINGLE precision "<<std::endl;
  if ( sizeof(RealF)==8 )   std::cout << GridLogMessage<< "* DOUBLE precision "<<std::endl;
#ifdef GRID_OMP
  if ( WilsonKernelsStatic::Comms == WilsonKernelsStatic::CommsAndCompute ) std::cout << GridLogMessage<< "* Using Overlapped Comms/Compute" <<std::endl;
  if ( WilsonKernelsStatic::Comms == WilsonKernelsStatic::CommsThenCompute) std::cout << GridLogMessage<< "* Using sequential comms compute" <<std::endl;
#endif
  if ( WilsonKernelsStatic::Opt == WilsonKernelsStatic::OptGeneric   ) std::cout << GridLogMessage<< "* Using GENERIC Nc WilsonKernels" <<std::endl;
  if ( WilsonKernelsStatic::Opt == WilsonKernelsStatic::OptHandUnroll) std::cout << GridLogMessage<< "* Using Nc=3       WilsonKernels" <<std::endl;
  if ( WilsonKernelsStatic::Opt == WilsonKernelsStatic::OptInlineAsm ) std::cout << GridLogMessage<< "* Using Asm Nc=3   WilsonKernels" <<std::endl;
  std::cout << GridLogMessage<< "*****************************************************************" <<std::endl;

  DomainWallFermionF Dw(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);
  int ncall =100;

  if (1) {
    FGrid->Barrier();
    Dw.Dhop(src,result,0);
    std::cout<<GridLogMessage<<"Called warmup"<<std::endl;
    double t0=usecond();
    for(int i=0;i<ncall;i++){
      Dw.Dhop(src1,result,0);
      Dw.Dhop(src,result,0);
      err = ref-result;
      std::cout<<GridLogMessage << "norm diff   "<< norm2(err)<<std::endl;
      assert (norm2(err)< 1.0e-4 );
    }
    double t1=usecond();
    FGrid->Barrier();

    double volume=Ls;  for(int mu=0;mu<Nd;mu++) volume=volume*latt4[mu];
    double flops=single_site_flops*volume*ncall;

    auto nsimd = vComplex::Nsimd();
    auto simdwidth = sizeof(vComplex);

    // RF: Nd Wilson * Ls, Nd gauge * Ls, Nc colors
    double data_rf = volume * ((2*Nd+1)*Nd*Nc + 2*Nd*Nc*Nc) * simdwidth / nsimd * ncall / (1024.*1024.*1024.);

    // mem: Nd Wilson * Ls, Nd gauge, Nc colors
    double data_mem = (volume * (2*Nd+1)*Nd*Nc + (volume/Ls) *2*Nd*Nc*Nc) * simdwidth / nsimd * ncall / (1024.*1024.*1024.);

    std::cout<<GridLogMessage << "Called Dw "<<ncall<<" times in "<<t1-t0<<" us"<<std::endl;
    //    std::cout<<GridLogMessage << "norm result "<< norm2(result)<<std::endl;
    //    std::cout<<GridLogMessage << "norm ref    "<< norm2(ref)<<std::endl;
    std::cout<<GridLogMessage << "mflop/s =   "<< flops/(t1-t0)<<std::endl;
    std::cout<<GridLogMessage << "mflop/s per rank =  "<< flops/(t1-t0)/NP<<std::endl;
    std::cout<<GridLogMessage << "mflop/s per node =  "<< flops/(t1-t0)/NN<<std::endl;
    std::cout<<GridLogMessage << "RF  GiB/s (base 2) =   "<< 1000000. * data_rf/((t1-t0))<<std::endl;
    std::cout<<GridLogMessage << "mem GiB/s (base 2) =   "<< 1000000. * data_mem/((t1-t0))<<std::endl;
    err = ref-result;
    std::cout<<GridLogMessage << "norm diff   "<< norm2(err)<<std::endl;
    //exit(0);

    if(( norm2(err)>1.0e-4) ) {

      /*
      std::cout << "RESULT\n " << result<<std::endl;
      std::cout << "REF   \n " << ref   <<std::endl;
      std::cout << "ERR   \n " << err   <<std::endl;
      */
      std::cout<<GridLogMessage << "WRONG RESULT" << std::endl;
      FGrid->Barrier();
      exit(-1);
    }
    assert (norm2(err)< 1.0e-4 );
  }

  if (1)
  { // Naive wilson dag implementation
    ref = Zero();
    for(int mu=0;mu<Nd;mu++){

      //    ref =  src - Gamma(Gamma::Algebra::GammaX)* src ; // 1+gamma_x
      tmp = Cshift(src,mu+1,1);
      {
	autoView( ref_v, ref, CpuWrite);
	autoView( tmp_v, tmp, CpuRead);
	autoView( U_v  , U[mu]  , CpuRead);
	for(int ss=0;ss<U[mu].Grid()->oSites();ss++){
	  for(int s=0;s<Ls;s++){
	    int i=s+Ls*ss;
	    ref_v[i]+= U_v[ss]*(tmp_v[i] + Gamma(Gmu[mu])*tmp_v[i]); ;
	  }
	}
      }
      
      {
	autoView( tmp_v  , tmp  , CpuWrite);
	autoView( U_v  , U[mu]  , CpuRead);
	autoView( src_v, src    , CpuRead);
	for(int ss=0;ss<U[mu].Grid()->oSites();ss++){
	  for(int s=0;s<Ls;s++){
	    tmp_v[Ls*ss+s] = adj(U_v[ss])*src_v[Ls*ss+s];
	  }
	}
      }
      //      tmp =adj(U[mu])*src;
      tmp =Cshift(tmp,mu+1,-1);
      {
	autoView( ref_v, ref, CpuWrite);
	autoView( tmp_v, tmp, CpuRead);
	for(int i=0;i<ref_v.size();i++){
	  ref_v[i]+= tmp_v[i] - Gamma(Gmu[mu])*tmp_v[i]; ;
	}
      }
    }
    ref = -0.5*ref;
  }
  //  dump=1;
  Dw.Dhop(src,result,1);
  std::cout << GridLogMessage << "Compare to naive wilson implementation Dag to verify correctness" << std::endl;
  std::cout<<GridLogMessage << "Called DwDag"<<std::endl;
  std::cout<<GridLogMessage << "norm dag result "<< norm2(result)<<std::endl;
  std::cout<<GridLogMessage << "norm dag ref    "<< norm2(ref)<<std::endl;
  err = ref-result;
  std::cout<<GridLogMessage << "norm dag diff   "<< norm2(err)<<std::endl;
  if((norm2(err)>1.0e-4)){
/*
	std::cout<< "DAG RESULT\n "  <<ref     << std::endl;
	std::cout<< "DAG sRESULT\n " <<result  << std::endl;
	std::cout<< "DAG ERR   \n "  << err    <<std::endl;
*/
  }
  LatticeFermionF src_e (FrbGrid);
  LatticeFermionF src_o (FrbGrid);
  LatticeFermionF r_e   (FrbGrid);
  LatticeFermionF r_o   (FrbGrid);
  LatticeFermionF r_eo  (FGrid);

  std::cout<<GridLogMessage << "Calling Deo and Doe and //assert Deo+Doe == Dunprec"<<std::endl;
  pickCheckerboard(Even,src_e,src);
  pickCheckerboard(Odd,src_o,src);

  std::cout<<GridLogMessage << "src_e"<<norm2(src_e)<<std::endl;
  std::cout<<GridLogMessage << "src_o"<<norm2(src_o)<<std::endl;


  // S-direction is INNERMOST and takes no part in the parity.
  std::cout << GridLogMessage<< "*********************************************************" <<std::endl;
  std::cout << GridLogMessage<< "* Benchmarking DomainWallFermionF::DhopEO                "<<std::endl;
  std::cout << GridLogMessage<< "* Vectorising space-time by "<<vComplexF::Nsimd()<<std::endl;
  if ( sizeof(RealF)==4 )   std::cout << GridLogMessage<< "* SINGLE precision "<<std::endl;
  if ( sizeof(RealF)==8 )   std::cout << GridLogMessage<< "* DOUBLE precision "<<std::endl;
#ifdef GRID_OMP
  if ( WilsonKernelsStatic::Comms == WilsonKernelsStatic::CommsAndCompute ) std::cout << GridLogMessage<< "* Using Overlapped Comms/Compute" <<std::endl;
  if ( WilsonKernelsStatic::Comms == WilsonKernelsStatic::CommsThenCompute) std::cout << GridLogMessage<< "* Using sequential comms compute" <<std::endl;
#endif
  if ( WilsonKernelsStatic::Opt == WilsonKernelsStatic::OptGeneric   ) std::cout << GridLogMessage<< "* Using GENERIC Nc WilsonKernels" <<std::endl;
  if ( WilsonKernelsStatic::Opt == WilsonKernelsStatic::OptHandUnroll) std::cout << GridLogMessage<< "* Using Nc=3       WilsonKernels" <<std::endl;
  if ( WilsonKernelsStatic::Opt == WilsonKernelsStatic::OptInlineAsm ) std::cout << GridLogMessage<< "* Using Asm Nc=3   WilsonKernels" <<std::endl;
  std::cout << GridLogMessage<< "*********************************************************" <<std::endl;
  {
    FGrid->Barrier();
    Dw.DhopEO(src_o,r_e,DaggerNo);
    double t0=usecond();
    for(int i=0;i<ncall;i++){
#ifdef CUDA_PROFILE
      if(i==10) cudaProfilerStart();
#endif
      Dw.DhopEO(src_o,r_e,DaggerNo);
#ifdef CUDA_PROFILE
      if(i==20) cudaProfilerStop();
#endif
    }
    double t1=usecond();
    FGrid->Barrier();

    double volume=Ls;  for(int mu=0;mu<Nd;mu++) volume=volume*latt4[mu];
    double flops=(single_site_flops*volume*ncall)/2.0;

    std::cout<<GridLogMessage << "Deo mflop/s =   "<< flops/(t1-t0)<<std::endl;
    std::cout<<GridLogMessage << "Deo mflop/s per rank   "<< flops/(t1-t0)/NP<<std::endl;
    std::cout<<GridLogMessage << "Deo mflop/s per node   "<< flops/(t1-t0)/NN<<std::endl;
  }
  Dw.DhopEO(src_o,r_e,DaggerNo);
  Dw.DhopOE(src_e,r_o,DaggerNo);
  Dw.Dhop  (src  ,result,DaggerNo);

  std::cout<<GridLogMessage << "r_e"<<norm2(r_e)<<std::endl;
  std::cout<<GridLogMessage << "r_o"<<norm2(r_o)<<std::endl;
  std::cout<<GridLogMessage << "res"<<norm2(result)<<std::endl;

  setCheckerboard(r_eo,r_o);
  setCheckerboard(r_eo,r_e);

  err = r_eo-result;
  std::cout<<GridLogMessage << "norm diff   "<< norm2(err)<<std::endl;
  if((norm2(err)>1.0e-4)){
    /*
	std::cout<< "Deo RESULT\n " <<r_eo << std::endl;
	std::cout<< "Deo REF\n " <<result  << std::endl;
	std::cout<< "Deo ERR   \n " << err <<std::endl;
    */
  }

  pickCheckerboard(Even,src_e,err);
  pickCheckerboard(Odd,src_o,err);
  std::cout<<GridLogMessage << "norm diff even  "<< norm2(src_e)<<std::endl;
  std::cout<<GridLogMessage << "norm diff odd   "<< norm2(src_o)<<std::endl;

  assert(norm2(src_e)<1.0e-4);
  assert(norm2(src_o)<1.0e-4);
  Grid_finalize();
  exit(0);
}
