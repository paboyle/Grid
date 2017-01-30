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

void benchDw(std::vector<int> & L, int Ls, int threads, int report =0 );
void benchsDw(std::vector<int> & L, int Ls, int threads, int report=0 );

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  std::cout << GridLogMessage<< "*****************************************************************" <<std::endl;
  std::cout << GridLogMessage<< "* Kernel options --dslash-generic, --dslash-unroll, --dslash-asm" <<std::endl;
  std::cout << GridLogMessage<< "*****************************************************************" <<std::endl;

  if ( WilsonKernelsStatic::Opt == WilsonKernelsStatic::OptGeneric   ) std::cout << GridLogMessage<< "* Using GENERIC Nc WilsonKernels" <<std::endl;
  if ( WilsonKernelsStatic::Opt == WilsonKernelsStatic::OptHandUnroll) std::cout << GridLogMessage<< "* Using Nc=3       WilsonKernels" <<std::endl;
  if ( WilsonKernelsStatic::Opt == WilsonKernelsStatic::OptInlineAsm ) std::cout << GridLogMessage<< "* Using Asm Nc=3   WilsonKernels" <<std::endl;
  std::cout << GridLogMessage<< "*****************************************************************" <<std::endl;

  const int Ls=8;
  int threads = GridThread::GetThreads();
  std::cout<<GridLogMessage << "Grid is setup to use "<<threads<<" threads"<<std::endl;
  std::cout<<GridLogMessage << "=========================================================================="<<std::endl;
  std::cout<<GridLogMessage << "= Benchmarking DWF"<<std::endl;
  std::cout<<GridLogMessage << "=========================================================================="<<std::endl;
  std::cout<<GridLogMessage << "Volume \t\t\tProcs \t Dw \t eoDw \t sDw \t eosDw (Mflop/s)  "<<std::endl;
  std::cout<<GridLogMessage << "=========================================================================="<<std::endl;

  int Lmax=16;
  int dmin=2;
  if ( getenv("LMAX") ) Lmax=atoi(getenv("LMAX"));
  if ( getenv("DMIN") ) dmin=atoi(getenv("DMIN"));
  for (int L=8;L<=Lmax;L*=2){
    std::vector<int> latt4(4,L);
    for(int d=4;d>dmin;d--){
      if ( d<=3 ) latt4[d]*=2;
      std::cout << GridLogMessage <<"\t";
      for(int d=0;d<Nd;d++){
	std::cout<<latt4[d]<<"x";
      }
      std::cout <<Ls<<"\t" ;
      benchDw (latt4,Ls,threads,0);
      benchsDw(latt4,Ls,threads,0);
      std::cout<<std::endl;
    }
  }
  std::cout<<GridLogMessage << "=========================================================================="<<std::endl;
  {
    std::vector<int> latt4(4,16);
    std::cout<<GridLogMessage << "16^4 Dw miss rate"<<std::endl;
    benchDw (latt4,Ls,threads,1);
    std::cout<<GridLogMessage << "16^4 sDw miss rate"<<std::endl;
    benchsDw(latt4,Ls,threads,1);
  }

  Grid_finalize();
}

#undef CHECK

void benchDw(std::vector<int> & latt4, int Ls, int threads,int report )
{
  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(latt4, GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  std::vector<int> seeds4({1,2,3,4});
  std::vector<int> seeds5({5,6,7,8});

#ifdef CHECK 
  GridParallelRNG          RNG4(UGrid);  RNG4.SeedFixedIntegers(seeds4);
  GridParallelRNG          RNG5(FGrid);  RNG5.SeedFixedIntegers(seeds5);
  LatticeFermion src   (FGrid); random(RNG5,src);
  LatticeGaugeField Umu(UGrid); 
  random(RNG4,Umu);
#else 
  LatticeFermion src   (FGrid); src=zero;
  LatticeGaugeField Umu(UGrid); Umu=zero;
#endif

  LatticeFermion result(FGrid); result=zero;
  LatticeFermion    ref(FGrid);    ref=zero;
  LatticeFermion    tmp(FGrid);
  LatticeFermion    err(FGrid);

  ColourMatrix cm = Complex(1.0,0.0);

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

#ifdef CHECK
  if (1) {

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
#endif

  RealD mass=0.1;
  RealD M5  =1.8;
  RealD NP = UGrid->_Nprocessors;

  DomainWallFermionR Dw(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);
  
  double t0=usecond();
  Dw.Dhop(src,result,0);
  double t1=usecond();

#ifdef TIMERS_OFF
    int ncall =10;
#else
  int ncall =1+(int) ((5.0*1000*1000)/(t1-t0));
#endif

  if (ncall < 5 ) exit(0);

  Dw.Dhop(src,result,0);

  PerformanceCounter Counter(8);
  Counter.Start();
  t0=usecond();
  for(int i=0;i<ncall;i++){
    Dw.Dhop(src,result,0);
  }
  t1=usecond();
  Counter.Stop();
  if ( report ) {
    Counter.Report();
  }
  
  if ( ! report ) {
    double volume=Ls;  for(int mu=0;mu<Nd;mu++) volume=volume*latt4[mu];
    double flops=1344*volume*ncall;
    std::cout <<"\t"<<NP<< "\t"<<flops/(t1-t0)<< "\t";
  }
  
#ifdef CHECK
  err = ref-result; 
  RealD errd = norm2(err);
  if ( errd> 1.0e-4 ) {
    std::cout<<GridLogMessage << "oops !!! norm diff   "<< norm2(err)<<std::endl;
    exit(-1);
  }
#endif
    
  LatticeFermion src_e (FrbGrid);
  LatticeFermion src_o (FrbGrid);
  LatticeFermion r_e   (FrbGrid);
  LatticeFermion r_o   (FrbGrid);
  LatticeFermion r_eo  (FGrid);
  
  pickCheckerboard(Even,src_e,src);
  pickCheckerboard(Odd,src_o,src);
  
  {
    Dw.DhopEO(src_o,r_e,DaggerNo);
    double t0=usecond();
    for(int i=0;i<ncall;i++){
      Dw.DhopEO(src_o,r_e,DaggerNo);
    }
    double t1=usecond();
    
    if(!report){
      double volume=Ls;  for(int mu=0;mu<Nd;mu++) volume=volume*latt4[mu];
      double flops=(1344.0*volume*ncall)/2;
      std::cout<< flops/(t1-t0);
    }
  }
}

#define CHECK_SDW
void benchsDw(std::vector<int> & latt4, int Ls, int threads, int report )
{

  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(latt4, GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  GridCartesian         * sUGrid   = SpaceTimeGrid::makeFourDimDWFGrid(latt4,GridDefaultMpi());
  GridRedBlackCartesian * sUrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(sUGrid);
  GridCartesian         * sFGrid   = SpaceTimeGrid::makeFiveDimDWFGrid(Ls,UGrid);
  GridRedBlackCartesian * sFrbGrid = SpaceTimeGrid::makeFiveDimDWFRedBlackGrid(Ls,UGrid);

  std::vector<int> seeds4({1,2,3,4});
  std::vector<int> seeds5({5,6,7,8});

#ifdef CHECK_SDW
  GridParallelRNG          RNG4(UGrid);  RNG4.SeedFixedIntegers(seeds4);
  GridParallelRNG          RNG5(FGrid);  RNG5.SeedFixedIntegers(seeds5);
  LatticeFermion src   (FGrid); random(RNG5,src);
  LatticeGaugeField Umu(UGrid); 
  random(RNG4,Umu);
#else 
  LatticeFermion src   (FGrid); src=zero;
  LatticeGaugeField Umu(UGrid); Umu=zero;
#endif

  LatticeFermion result(FGrid); result=zero;
  LatticeFermion    ref(FGrid);    ref=zero;
  LatticeFermion    tmp(FGrid);
  LatticeFermion    err(FGrid);

  ColourMatrix cm = Complex(1.0,0.0);

  LatticeGaugeField Umu5d(FGrid); 

  // replicate across fifth dimension
  for(int ss=0;ss<Umu._grid->oSites();ss++){
    for(int s=0;s<Ls;s++){
      Umu5d._odata[Ls*ss+s] = Umu._odata[ss];
    }
  }

  RealD mass=0.1;
  RealD M5  =1.8;

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

  double t0=usecond();
  sDw.Dhop(ssrc,sresult,0);
  double t1=usecond();

#ifdef TIMERS_OFF
  int ncall =10;
#else 
  int ncall =1+(int) ((5.0*1000*1000)/(t1-t0));
#endif

  PerformanceCounter Counter(8);
  Counter.Start();
  t0=usecond();
  for(int i=0;i<ncall;i++){
    sDw.Dhop(ssrc,sresult,0);
  }
  t1=usecond();
  Counter.Stop();
  
  if ( report ) {
    Counter.Report();
  } else { 
    double volume=Ls;  for(int mu=0;mu<Nd;mu++) volume=volume*latt4[mu];
    double flops=1344*volume*ncall;
    std::cout<<"\t"<< flops/(t1-t0);
  }

  LatticeFermion sr_eo(sFGrid);
  LatticeFermion serr(sFGrid);
  
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
  
  sDw.DhopEO(ssrc_o,sr_e,DaggerNo);
  PerformanceCounter CounterSdw(8);
  CounterSdw.Start();
  t0=usecond();
  for(int i=0;i<ncall;i++){
    __SSC_START;
    sDw.DhopEO(ssrc_o,sr_e,DaggerNo);
    __SSC_STOP;
  }
  t1=usecond();
  CounterSdw.Stop();

  if ( report ) { 
    CounterSdw.Report();
  } else {
    double volume=Ls;  for(int mu=0;mu<Nd;mu++) volume=volume*latt4[mu];
    double flops=(1344.0*volume*ncall)/2;
    std::cout<<"\t"<< flops/(t1-t0);
  }
}


