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
#include <Grid.h>
#include <PerfCount.h>

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

bool overlapComms = false;
typedef WilsonFermion5D<DomainWallRedBlack5dImplR> WilsonFermion5DR;
typedef WilsonFermion5D<DomainWallRedBlack5dImplF> WilsonFermion5DF;
typedef WilsonFermion5D<DomainWallRedBlack5dImplD> WilsonFermion5DD;


int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  if( GridCmdOptionExists(argv,argv+argc,"--asynch") ){
    overlapComms = true;
  }

  int threads = GridThread::GetThreads();
  std::cout<<GridLogMessage << "Grid is setup to use "<<threads<<" threads"<<std::endl;

  std::vector<int> latt4 = GridDefaultLatt();
  const int Ls=16;
  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  std::cout << GridLogMessage << "Making s innermost grids"<<std::endl;
  GridCartesian         * sUGrid   = SpaceTimeGrid::makeFourDimDWFGrid(GridDefaultLatt(),GridDefaultMpi());
  GridCartesian         * sFGrid   = SpaceTimeGrid::makeFiveDimDWFGrid(Ls,UGrid);
  std::cout << GridLogMessage << "Making s innermost rb grids"<<std::endl;
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

  ColourMatrix cm = Complex(1.0,0.0);

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

  typename DomainWallFermionR::ImplParams params; 
  params.overlapCommsCompute = overlapComms;
  
  RealD NP = UGrid->_Nprocessors;

  for(int doasm=0;doasm<1;doasm++){

    QCD::WilsonKernelsStatic::AsmOpt=doasm;

  DomainWallFermionR Dw(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,params);
  
  std::cout<<GridLogMessage << "Calling Dw"<<std::endl;
  int ncall =50;
  if (1) {

    double t0=usecond();
    for(int i=0;i<ncall;i++){
      Dw.Dhop(src,result,0);
    }
    double t1=usecond();
    
    double volume=Ls;  for(int mu=0;mu<Nd;mu++) volume=volume*latt4[mu];
    double flops=1344*volume*ncall;

    std::cout<<GridLogMessage << "Called Dw "<<ncall<<" times in "<<t1-t0<<" us"<<std::endl;
    std::cout<<GridLogMessage << "norm result "<< norm2(result)<<std::endl;
    std::cout<<GridLogMessage << "norm ref    "<< norm2(ref)<<std::endl;
    std::cout<<GridLogMessage << "mflop/s =   "<< flops/(t1-t0)<<std::endl;
    std::cout<<GridLogMessage << "mflop/s per node =  "<< flops/(t1-t0)/NP<<std::endl;
    err = ref-result; 
    std::cout<<GridLogMessage << "norm diff   "<< norm2(err)<<std::endl;
    //    Dw.Report();
  }

  if (1)
  {
    typedef WilsonFermion5D<DomainWallRedBlack5dImplF> WilsonFermion5DF;
    LatticeFermionF ssrc(sFGrid);
    LatticeFermionF sref(sFGrid);
    LatticeFermionF sresult(sFGrid);
    WilsonFermion5DF sDw(1,Umu,*sFGrid,*sFrbGrid,*sUGrid,M5,params);
  
    for(int x=0;x<latt4[0];x++){
    for(int y=0;y<latt4[1];y++){
    for(int z=0;z<latt4[2];z++){
    for(int t=0;t<latt4[3];t++){
    for(int s=0;s<Ls;s++){
      std::vector<int> site({s,x,y,z,t});
      SpinColourVectorF tmp;
      peekSite(tmp,src,site);
      pokeSite(tmp,ssrc,site);
    }}}}}

    double t0=usecond();
    for(int i=0;i<ncall;i++){
      __SSC_START;
      sDw.Dhop(ssrc,sresult,0);
      __SSC_STOP;
    }
    double t1=usecond();
    double volume=Ls;  for(int mu=0;mu<Nd;mu++) volume=volume*latt4[mu];
    double flops=1344*volume*ncall;

    std::cout<<GridLogMessage << "Called Dw sinner "<<ncall<<" times in "<<t1-t0<<" us"<<std::endl;
    std::cout<<GridLogMessage << "mflop/s =   "<< flops/(t1-t0)<<std::endl;
    std::cout<<GridLogMessage << "mflop/s per node =  "<< flops/(t1-t0)/NP<<std::endl;
    //  sDw.Report();
  
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



    RealF sum=0;
    for(int x=0;x<latt4[0];x++){
    for(int y=0;y<latt4[1];y++){
    for(int z=0;z<latt4[2];z++){
    for(int t=0;t<latt4[3];t++){
    for(int s=0;s<Ls;s++){
      std::vector<int> site({s,x,y,z,t});
      SpinColourVectorF normal, simd;
      peekSite(normal,result,site);
      peekSite(simd,sresult,site);
      sum=sum+norm2(normal-simd);
      //      std::cout << "site "<<x<<","<<y<<","<<z<<","<<t<<","<<s<<" "<<norm2(normal-simd)<<std::endl;
      //      std::cout << "site "<<x<<","<<y<<","<<z<<","<<t<<","<<s<<" "<<normal<<std::endl;
      //      std::cout << "site "<<x<<","<<y<<","<<z<<","<<t<<","<<s<<" "<<simd<<std::endl;
    }}}}}
    std::cout<<" difference between normal and simd is "<<sum<<std::endl;


    if (1) {

      LatticeFermionF sr_eo(sFGrid);
      LatticeFermionF serr(sFGrid);

      LatticeFermion ssrc_e (sFrbGrid);
      LatticeFermion ssrc_o (sFrbGrid);
      LatticeFermion sr_e   (sFrbGrid);
      LatticeFermion sr_o   (sFrbGrid);

      pickCheckerboard(Even,ssrc_e,ssrc);
      pickCheckerboard(Odd,ssrc_o,ssrc);

      setCheckerboard(sr_eo,ssrc_o);
      setCheckerboard(sr_eo,ssrc_e);
      serr = sr_eo-ssrc; 
      std::cout<<GridLogMessage << "EO src norm diff   "<< norm2(serr)<<std::endl;

      sr_e = zero;
      sr_o = zero;

      double t0=usecond();
      for(int i=0;i<ncall;i++){
	sDw.DhopEO(ssrc_o,sr_e,DaggerNo);
      }
      double t1=usecond();

      double volume=Ls;  for(int mu=0;mu<Nd;mu++) volume=volume*latt4[mu];
      double flops=(1344.0*volume*ncall)/2;

      std::cout<<GridLogMessage << "sDeo mflop/s =   "<< flops/(t1-t0)<<std::endl;
      std::cout<<GridLogMessage << "sDeo mflop/s per node   "<< flops/(t1-t0)/NP<<std::endl;

      sDw.DhopEO(ssrc_o,sr_e,DaggerNo);
      sDw.DhopOE(ssrc_e,sr_o,DaggerNo);
      sDw.Dhop  (ssrc  ,sresult,DaggerNo);

      pickCheckerboard(Even,ssrc_e,sresult);
      pickCheckerboard(Odd ,ssrc_o,sresult);
      ssrc_e = ssrc_e - sr_e;
      std::cout<<GridLogMessage << "sE norm diff   "<< norm2(ssrc_e)<<std::endl;
      ssrc_o = ssrc_o - sr_o;
      std::cout<<GridLogMessage << "sO norm diff   "<< norm2(ssrc_o)<<std::endl;
    }


  }

  if (1)
  { // Naive wilson dag implementation
    ref = zero;
    for(int mu=0;mu<Nd;mu++){

      //    ref =  src - Gamma(Gamma::GammaX)* src ; // 1+gamma_x
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
  std::cout<<GridLogMessage << "Called DwDag"<<std::endl;
  std::cout<<GridLogMessage << "norm result "<< norm2(result)<<std::endl;
  std::cout<<GridLogMessage << "norm ref    "<< norm2(ref)<<std::endl;
  err = ref-result; 
  std::cout<<GridLogMessage << "norm diff   "<< norm2(err)<<std::endl;

  LatticeFermion src_e (FrbGrid);
  LatticeFermion src_o (FrbGrid);
  LatticeFermion r_e   (FrbGrid);
  LatticeFermion r_o   (FrbGrid);
  LatticeFermion r_eo  (FGrid);


  std::cout<<GridLogMessage << "Calling Deo and Doe"<<std::endl;
  pickCheckerboard(Even,src_e,src);
  pickCheckerboard(Odd,src_o,src);

  std::cout<<GridLogMessage << "src_e"<<norm2(src_e)<<std::endl;
  std::cout<<GridLogMessage << "src_o"<<norm2(src_o)<<std::endl;

  {
    double t0=usecond();
    for(int i=0;i<ncall;i++){
      Dw.DhopEO(src_o,r_e,DaggerNo);
    }
    double t1=usecond();
    
    double volume=Ls;  for(int mu=0;mu<Nd;mu++) volume=volume*latt4[mu];
    double flops=(1344.0*volume*ncall)/2;

    std::cout<<GridLogMessage << "Deo mflop/s =   "<< flops/(t1-t0)<<std::endl;
    std::cout<<GridLogMessage << "Deo mflop/s per node   "<< flops/(t1-t0)/NP<<std::endl;
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

  pickCheckerboard(Even,src_e,err);
  pickCheckerboard(Odd,src_o,err);
  std::cout<<GridLogMessage << "norm diff even  "<< norm2(src_e)<<std::endl;
  std::cout<<GridLogMessage << "norm diff odd   "<< norm2(src_o)<<std::endl;


  }

  Grid_finalize();
}
