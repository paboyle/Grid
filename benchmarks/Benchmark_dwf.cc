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

typedef WilsonFermion5D<DomainWallVec5dImplR> WilsonFermion5DR;
typedef WilsonFermion5D<DomainWallVec5dImplF> WilsonFermion5DF;
typedef WilsonFermion5D<DomainWallVec5dImplD> WilsonFermion5DD;


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

  std::cout << GridLogMessage << "Making s innermost grids"<<std::endl;
  GridCartesian         * sUGrid   = SpaceTimeGrid::makeFourDimDWFGrid(GridDefaultLatt(),GridDefaultMpi());
  GridRedBlackCartesian * sUrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(sUGrid);
  GridCartesian         * sFGrid   = SpaceTimeGrid::makeFiveDimDWFGrid(Ls,UGrid);
  GridRedBlackCartesian * sFrbGrid = SpaceTimeGrid::makeFiveDimDWFRedBlackGrid(Ls,UGrid);

  std::vector<int> seeds4({1,2,3,4});
  std::vector<int> seeds5({5,6,7,8});

  GridParallelRNG          RNG4(UGrid);  RNG4.SeedFixedIntegers(seeds4);
  GridParallelRNG          RNG5(FGrid);  RNG5.SeedFixedIntegers(seeds5);

  LatticeFermion src   (FGrid); random(RNG5,src);
  LatticeFermion result(FGrid); result=zero;
  LatticeFermion    ref(FGrid);    ref=zero;
  LatticeFermion    tmp(FGrid);
  LatticeFermion    err(FGrid);

  LatticeGaugeField Umu(UGrid); 
  random(RNG4,Umu);

  LatticeGaugeField Umu5d(FGrid); 

  // replicate across fifth dimension
  for(int ss=0;ss<Umu._grid->oSites();ss++){
    for(int s=0;s<Ls;s++){
      Umu5d._odata[Ls*ss+s] = Umu._odata[ss];
    }
  }

  ////////////////////////////////////
  // Naive wilson implementation
  ////////////////////////////////////
  std::vector<LatticeColourMatrix> U(4,FGrid);
  for(int mu=0;mu<Nd;mu++){
    U[mu] = PeekIndex<LorentzIndex>(Umu5d,mu);
  }

  if (1)
  {
    ref = zero;
    for(int mu=0;mu<Nd;mu++){

      tmp = U[mu]*Cshift(src,mu+1,1);
      ref=ref + tmp - Gamma(Gmu[mu])*tmp;

      tmp =adj(U[mu])*src;
      tmp =Cshift(tmp,mu+1,-1);
      ref=ref + tmp + Gamma(Gmu[mu])*tmp;
    }
    ref = -0.5*ref;
  }

  RealD mass=0.1;
  RealD M5  =1.8;

  RealD NP = UGrid->_Nprocessors;

  DomainWallFermionR Dw(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);

  std::cout << GridLogMessage<< "*****************************************************************" <<std::endl;
  std::cout << GridLogMessage<< "* Kernel options --dslash-generic, --dslash-unroll, --dslash-asm" <<std::endl;
  std::cout << GridLogMessage<< "*****************************************************************" <<std::endl;
  std::cout << GridLogMessage<< "*****************************************************************" <<std::endl;
  std::cout << GridLogMessage<< "* Benchmarking DomainWallFermionR::Dhop                  "<<std::endl;
  std::cout << GridLogMessage<< "* Vectorising space-time by "<<vComplex::Nsimd()<<std::endl;
  if ( sizeof(Real)==4 )   std::cout << GridLogMessage<< "* SINGLE precision "<<std::endl;
  if ( sizeof(Real)==8 )   std::cout << GridLogMessage<< "* DOUBLE precision "<<std::endl;
  if ( WilsonKernelsStatic::Opt == WilsonKernelsStatic::OptGeneric   ) std::cout << GridLogMessage<< "* Using GENERIC Nc WilsonKernels" <<std::endl;
  if ( WilsonKernelsStatic::Opt == WilsonKernelsStatic::OptHandUnroll) std::cout << GridLogMessage<< "* Using Nc=3       WilsonKernels" <<std::endl;
  if ( WilsonKernelsStatic::Opt == WilsonKernelsStatic::OptInlineAsm ) std::cout << GridLogMessage<< "* Using Asm Nc=3   WilsonKernels" <<std::endl;
  std::cout << GridLogMessage<< "*****************************************************************" <<std::endl;

  int ncall =100;
  if (1) {
    FGrid->Barrier();
    Dw.ZeroCounters();
    double t0=usecond();
    for(int i=0;i<ncall;i++){
      __SSC_START;
      Dw.Dhop(src,result,0);
      __SSC_STOP;
    }
    double t1=usecond();
    FGrid->Barrier();
    
    double volume=Ls;  for(int mu=0;mu<Nd;mu++) volume=volume*latt4[mu];
    double flops=1344*volume*ncall;

    std::cout<<GridLogMessage << "Called Dw "<<ncall<<" times in "<<t1-t0<<" us"<<std::endl;
    std::cout<<GridLogMessage << "norm result "<< norm2(result)<<std::endl;
    std::cout<<GridLogMessage << "norm ref    "<< norm2(ref)<<std::endl;
    std::cout<<GridLogMessage << "mflop/s =   "<< flops/(t1-t0)<<std::endl;
    std::cout<<GridLogMessage << "mflop/s per rank =  "<< flops/(t1-t0)/NP<<std::endl;
    err = ref-result; 
    std::cout<<GridLogMessage << "norm diff   "<< norm2(err)<<std::endl;
    assert (norm2(err)< 1.0e-4 );
    Dw.Report();
  }

  if (1)
  {

    std::cout << GridLogMessage<< "*********************************************************" <<std::endl;
    std::cout << GridLogMessage<< "* Benchmarking WilsonFermion5D<DomainWallVec5dImplR>::Dhop "<<std::endl;
    std::cout << GridLogMessage<< "* Vectorising fifth dimension by "<<vComplex::Nsimd()<<std::endl;
    if ( sizeof(Real)==4 )   std::cout << GridLogMessage<< "* SINGLE precision "<<std::endl;
    if ( sizeof(Real)==8 )   std::cout << GridLogMessage<< "* DOUBLE precision "<<std::endl;
    if ( WilsonKernelsStatic::Opt == WilsonKernelsStatic::OptGeneric   ) std::cout << GridLogMessage<< "* Using GENERIC Nc WilsonKernels" <<std::endl;
    if ( WilsonKernelsStatic::Opt == WilsonKernelsStatic::OptHandUnroll) std::cout << GridLogMessage<< "* Using Nc=3       WilsonKernels" <<std::endl;
    if ( WilsonKernelsStatic::Opt == WilsonKernelsStatic::OptInlineAsm ) std::cout << GridLogMessage<< "* Using Asm Nc=3   WilsonKernels" <<std::endl;
    std::cout << GridLogMessage<< "*********************************************************" <<std::endl;

    typedef WilsonFermion5D<DomainWallVec5dImplR> WilsonFermion5DR;
    LatticeFermion ssrc(sFGrid);
    LatticeFermion sref(sFGrid);
    LatticeFermion sresult(sFGrid);

    WilsonFermion5DR sDw(Umu,*sFGrid,*sFrbGrid,*sUGrid,*sUrbGrid,M5);
  
    for(int x=0;x<latt4[0];x++){
    for(int y=0;y<latt4[1];y++){
    for(int z=0;z<latt4[2];z++){
    for(int t=0;t<latt4[3];t++){
    for(int s=0;s<Ls;s++){
      std::vector<int> site({s,x,y,z,t});
      SpinColourVector tmp;
      peekSite(tmp,src,site);
      pokeSite(tmp,ssrc,site);
    }}}}}
    std::cout<<GridLogMessage<< "src norms "<< norm2(src)<<" " <<norm2(ssrc)<<std::endl;
    FGrid->Barrier();
    double t0=usecond();
    sDw.ZeroCounters();
    for(int i=0;i<ncall;i++){
      __SSC_START;
      sDw.Dhop(ssrc,sresult,0);
      __SSC_STOP;
    }
    double t1=usecond();
    FGrid->Barrier();
    double volume=Ls;  for(int mu=0;mu<Nd;mu++) volume=volume*latt4[mu];
    double flops=1344*volume*ncall;

    std::cout<<GridLogMessage << "Called Dw s_inner "<<ncall<<" times in "<<t1-t0<<" us"<<std::endl;
    std::cout<<GridLogMessage << "mflop/s =   "<< flops/(t1-t0)<<std::endl;
    std::cout<<GridLogMessage << "mflop/s per rank =  "<< flops/(t1-t0)/NP<<std::endl;
    sDw.Report();
  
    if(0){
      for(int i=0;i< PerformanceCounter::NumTypes(); i++ ){
	sDw.Dhop(ssrc,sresult,0);
	PerformanceCounter Counter(i);
	Counter.Start();
	sDw.Dhop(ssrc,sresult,0);
	Counter.Stop();
	Counter.Report();
      }
    }

    std::cout<<GridLogMessage<< "res norms "<< norm2(result)<<" " <<norm2(sresult)<<std::endl;

    RealD sum=0;
    for(int x=0;x<latt4[0];x++){
    for(int y=0;y<latt4[1];y++){
    for(int z=0;z<latt4[2];z++){
    for(int t=0;t<latt4[3];t++){
    for(int s=0;s<Ls;s++){
      std::vector<int> site({s,x,y,z,t});
      SpinColourVector normal, simd;
      peekSite(normal,result,site);
      peekSite(simd,sresult,site);
      sum=sum+norm2(normal-simd);
      if (norm2(normal-simd) > 1.0e-6 ) {
	std::cout << "site "<<x<<","<<y<<","<<z<<","<<t<<","<<s<<" "<<norm2(normal-simd)<<std::endl;
	std::cout << "site "<<x<<","<<y<<","<<z<<","<<t<<","<<s<<" normal "<<normal<<std::endl;
	std::cout << "site "<<x<<","<<y<<","<<z<<","<<t<<","<<s<<" simd   "<<simd<<std::endl;
      }
    }}}}}
    std::cout<<GridLogMessage<<" difference between normal and simd is "<<sum<<std::endl;
    assert (sum< 1.0e-4 );


    if (1) {

      LatticeFermion sr_eo(sFGrid);

      LatticeFermion ssrc_e (sFrbGrid);
      LatticeFermion ssrc_o (sFrbGrid);
      LatticeFermion sr_e   (sFrbGrid);
      LatticeFermion sr_o   (sFrbGrid);

      pickCheckerboard(Even,ssrc_e,ssrc);
      pickCheckerboard(Odd,ssrc_o,ssrc);

      setCheckerboard(sr_eo,ssrc_o);
      setCheckerboard(sr_eo,ssrc_e);

      sr_e = zero;
      sr_o = zero;

      std::cout << GridLogMessage<< "*********************************************************" <<std::endl;
      std::cout << GridLogMessage<< "* Benchmarking WilsonFermion5D<DomainWallVec5dImplR>::DhopEO "<<std::endl;
      std::cout << GridLogMessage<< "* Vectorising fifth dimension by "<<vComplex::Nsimd()<<std::endl;
      if ( sizeof(Real)==4 )   std::cout << GridLogMessage<< "* SINGLE precision "<<std::endl;
      if ( sizeof(Real)==8 )   std::cout << GridLogMessage<< "* DOUBLE precision "<<std::endl;
      if ( WilsonKernelsStatic::Opt == WilsonKernelsStatic::OptGeneric   ) std::cout << GridLogMessage<< "* Using GENERIC Nc WilsonKernels" <<std::endl;
      if ( WilsonKernelsStatic::Opt == WilsonKernelsStatic::OptHandUnroll) std::cout << GridLogMessage<< "* Using Nc=3       WilsonKernels" <<std::endl;
      if ( WilsonKernelsStatic::Opt == WilsonKernelsStatic::OptInlineAsm ) std::cout << GridLogMessage<< "* Using Asm Nc=3   WilsonKernels" <<std::endl;
      std::cout << GridLogMessage<< "*********************************************************" <<std::endl;

      FGrid->Barrier();
      sDw.ZeroCounters();
      sDw.stat.init("DhopEO");
      double t0=usecond();
      for (int i = 0; i < ncall; i++) {
        sDw.DhopEO(ssrc_o, sr_e, DaggerNo);
      }
      double t1=usecond();
      FGrid->Barrier();
      sDw.stat.print();

      double volume=Ls;  for(int mu=0;mu<Nd;mu++) volume=volume*latt4[mu];
      double flops=(1344.0*volume*ncall)/2;

      std::cout<<GridLogMessage << "sDeo mflop/s =   "<< flops/(t1-t0)<<std::endl;
      std::cout<<GridLogMessage << "sDeo mflop/s per rank   "<< flops/(t1-t0)/NP<<std::endl;
      sDw.Report();

      sDw.DhopEO(ssrc_o,sr_e,DaggerNo);
      sDw.DhopOE(ssrc_e,sr_o,DaggerNo);
      sDw.Dhop  (ssrc  ,sresult,DaggerNo);

      pickCheckerboard(Even,ssrc_e,sresult);
      pickCheckerboard(Odd ,ssrc_o,sresult);
      ssrc_e = ssrc_e - sr_e;
      RealD error = norm2(ssrc_e);

      std::cout<<GridLogMessage << "sE norm diff   "<< norm2(ssrc_e)<< "  vec nrm"<<norm2(sr_e) <<std::endl;
      ssrc_o = ssrc_o - sr_o;

      error+= norm2(ssrc_o);
      std::cout<<GridLogMessage << "sO norm diff   "<< norm2(ssrc_o)<< "  vec nrm"<<norm2(sr_o) <<std::endl;
      if(error>1.0e-4) { 
	setCheckerboard(ssrc,ssrc_o);
	setCheckerboard(ssrc,ssrc_e);
	std::cout<< ssrc << std::endl;
      }
    }


  }

  if (1)
  { // Naive wilson dag implementation
    ref = zero;
    for(int mu=0;mu<Nd;mu++){

      //    ref =  src - Gamma(Gamma::Algebra::GammaX)* src ; // 1+gamma_x
      tmp = U[mu]*Cshift(src,mu+1,1);
      for(int i=0;i<ref._odata.size();i++){
  ref._odata[i]+= tmp._odata[i] + Gamma(Gmu[mu])*tmp._odata[i]; ;
      }

      tmp =adj(U[mu])*src;
      tmp =Cshift(tmp,mu+1,-1);
      for(int i=0;i<ref._odata.size();i++){
  ref._odata[i]+= tmp._odata[i] - Gamma(Gmu[mu])*tmp._odata[i]; ;
      }
    }
    ref = -0.5*ref;
  }
  Dw.Dhop(src,result,1);
  std::cout << GridLogMessage << "Compare to naive wilson implementation Dag to verify correctness" << std::endl;
  std::cout<<GridLogMessage << "Called DwDag"<<std::endl;
  std::cout<<GridLogMessage << "norm result "<< norm2(result)<<std::endl;
  std::cout<<GridLogMessage << "norm ref    "<< norm2(ref)<<std::endl;
  err = ref-result; 
  std::cout<<GridLogMessage << "norm diff   "<< norm2(err)<<std::endl;
  assert(norm2(err)<1.0e-4);
  LatticeFermion src_e (FrbGrid);
  LatticeFermion src_o (FrbGrid);
  LatticeFermion r_e   (FrbGrid);
  LatticeFermion r_o   (FrbGrid);
  LatticeFermion r_eo  (FGrid);


  std::cout<<GridLogMessage << "Calling Deo and Doe and assert Deo+Doe == Dunprec"<<std::endl;
  pickCheckerboard(Even,src_e,src);
  pickCheckerboard(Odd,src_o,src);

  std::cout<<GridLogMessage << "src_e"<<norm2(src_e)<<std::endl;
  std::cout<<GridLogMessage << "src_o"<<norm2(src_o)<<std::endl;

  std::cout << GridLogMessage<< "*********************************************************" <<std::endl;
  std::cout << GridLogMessage<< "* Benchmarking DomainWallFermionR::DhopEO                "<<std::endl;
  std::cout << GridLogMessage<< "* Vectorising space-time by "<<vComplex::Nsimd()<<std::endl;
  if ( sizeof(Real)==4 )   std::cout << GridLogMessage<< "* SINGLE precision "<<std::endl;
  if ( sizeof(Real)==8 )   std::cout << GridLogMessage<< "* DOUBLE precision "<<std::endl;
  if ( WilsonKernelsStatic::Opt == WilsonKernelsStatic::OptGeneric   ) std::cout << GridLogMessage<< "* Using GENERIC Nc WilsonKernels" <<std::endl;
  if ( WilsonKernelsStatic::Opt == WilsonKernelsStatic::OptHandUnroll) std::cout << GridLogMessage<< "* Using Nc=3       WilsonKernels" <<std::endl;
  if ( WilsonKernelsStatic::Opt == WilsonKernelsStatic::OptInlineAsm ) std::cout << GridLogMessage<< "* Using Asm Nc=3   WilsonKernels" <<std::endl;
  std::cout << GridLogMessage<< "*********************************************************" <<std::endl;
  {
    Dw.ZeroCounters();
    FGrid->Barrier();
    double t0=usecond();
    for(int i=0;i<ncall;i++){
      Dw.DhopEO(src_o,r_e,DaggerNo);
    }
    double t1=usecond();
    FGrid->Barrier();
    
    double volume=Ls;  for(int mu=0;mu<Nd;mu++) volume=volume*latt4[mu];
    double flops=(1344.0*volume*ncall)/2;

    std::cout<<GridLogMessage << "Deo mflop/s =   "<< flops/(t1-t0)<<std::endl;
    std::cout<<GridLogMessage << "Deo mflop/s per rank   "<< flops/(t1-t0)/NP<<std::endl;
    Dw.Report();
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
  assert(norm2(err)<1.0e-4);

  pickCheckerboard(Even,src_e,err);
  pickCheckerboard(Odd,src_o,err);
  std::cout<<GridLogMessage << "norm diff even  "<< norm2(src_e)<<std::endl;
  std::cout<<GridLogMessage << "norm diff odd   "<< norm2(src_o)<<std::endl;
  assert(norm2(src_e)<1.0e-4);
  assert(norm2(src_o)<1.0e-4);

  Grid_finalize();
}
