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

  const int Ls=12;
  int threads = GridThread::GetThreads();
  std::cout<<GridLogMessage << "Grid is setup to use "<<threads<<" threads"<<std::endl;
  std::cout<<GridLogMessage << "=========================================================================="<<std::endl;
  std::cout<<GridLogMessage << "= Benchmarking DWF"<<std::endl;
  std::cout<<GridLogMessage << "=========================================================================="<<std::endl;
  std::cout<<GridLogMessage << "Volume \t\t\tProcs \t Dw \t eoDw   "<<std::endl;
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
      std::cout<<std::endl;
    }
  }
  std::cout<<GridLogMessage << "=========================================================================="<<std::endl;
#ifdef PERFCOUNT
  {
    std::vector<int> latt4(4,16);
    std::cout<<GridLogMessage << "16^4 Dw miss rate"<<std::endl;
    benchDw (latt4,Ls,threads,1);
  }
#endif
  Grid_finalize();
}

#undef CHECK

void benchDw(std::vector<int> & latt4, int Ls, int threads,int report )
{
  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(latt4, GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);
  long unsigned int single_site_flops = 8*Nc*(7+16*Nc);

  std::vector<int> seeds4({1,2,3,4});
  std::vector<int> seeds5({5,6,7,8});

#ifdef CHECK 
  GridParallelRNG          RNG4(UGrid);  RNG4.SeedFixedIntegers(seeds4);
  GridParallelRNG          RNG5(FGrid);  RNG5.SeedFixedIntegers(seeds5);
  LatticeFermion src   (FGrid); random(RNG5,src);
  LatticeGaugeField Umu(UGrid); 
  random(RNG4,Umu);
#else 
  LatticeFermion src   (FGrid); src=Zero();
  LatticeGaugeField Umu(UGrid); Umu=Zero();
#endif

  LatticeFermion result(FGrid); result=Zero();
  LatticeFermion    ref(FGrid);    ref=Zero();
  LatticeFermion    tmp(FGrid);
  LatticeFermion    err(FGrid);

  ColourMatrix cm = Complex(1.0,0.0);

  LatticeGaugeField Umu5d(FGrid); 

  // replicate across fifth dimension
  auto Umu5d_v = Umu5d.View();
  auto Umu_v   = Umu.View();
  for(int ss=0;ss<Umu.Grid()->oSites();ss++){
    for(int s=0;s<Ls;s++){
      Umu5d_v[Ls*ss+s] = Umu_v[ss];
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

    ref = Zero();
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
#ifdef PERFCOUNT
  PerformanceCounter Counter(8);
  Counter.Start();
#endif
  t0=usecond();
  for(int i=0;i<ncall;i++){
    Dw.Dhop(src,result,0);
  }
  t1=usecond();
#ifdef PERFCOUNT
  Counter.Stop();
  if ( report ) {
    Counter.Report();
  }
#endif  
  if ( ! report ) {
    double volume=Ls;  for(int mu=0;mu<Nd;mu++) volume=volume*latt4[mu];
    double flops=single_site_flops*volume*ncall;
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
      double flops=(single_site_flops*volume*ncall)/2.0;
      std::cout<< flops/(t1-t0);
    }
  }
}



