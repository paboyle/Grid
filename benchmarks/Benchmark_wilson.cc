    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./benchmarks/Benchmark_wilson.cc

    Copyright (C) 2018

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
 ;


#include "Grid/util/Profiling.h"

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

bool overlapComms = false;
bool perfProfiling = false;

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  if( GridCmdOptionExists(argv,argv+argc,"--asynch") ){
    overlapComms = true;
  }
  if( GridCmdOptionExists(argv,argv+argc,"--perf") ){
    perfProfiling = true;
  }

  long unsigned int single_site_flops = 8*Nc*(7+16*Nc);


  auto latt_size   = GridDefaultLatt();
  auto simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
  auto mpi_layout  = GridDefaultMpi();

  GridCartesian               Grid(latt_size,simd_layout,mpi_layout);
  GridRedBlackCartesian     RBGrid(&Grid);

  int threads = GridThread::GetThreads();

  GridLogLayout();

  std::cout<<GridLogMessage << "Grid floating point word size is REALF"<< sizeof(RealF)<<std::endl;
  std::cout<<GridLogMessage << "Grid floating point word size is REALD"<< sizeof(RealD)<<std::endl;
  std::cout<<GridLogMessage << "Grid floating point word size is REAL"<< sizeof(Real)<<std::endl;
  std::cout<<GridLogMessage << "Grid number of colours : "<< Nc <<std::endl;
  std::cout<<GridLogMessage << "Benchmarking Wilson operator in the fundamental representation" << std::endl;


  std::vector<int> seeds({1,2,3,4});
  GridParallelRNG          pRNG(&Grid);
  pRNG.SeedFixedIntegers(seeds);
  //  pRNG.SeedFixedIntegers(std::vector<int>({45,12,81,9});

  LatticeFermion src   (&Grid); random(pRNG,src);
  LatticeFermion result(&Grid); result=Zero();
  LatticeFermion    ref(&Grid);    ref=Zero();
  LatticeFermion    tmp(&Grid);    tmp=Zero();
  LatticeFermion    err(&Grid);    tmp=Zero();
  LatticeGaugeField Umu(&Grid); random(pRNG,Umu);
  std::vector<LatticeColourMatrix> U(4,&Grid);

  double volume=1;
  for(int mu=0;mu<Nd;mu++){
    volume=volume*latt_size[mu];
  }

  // Only one non-zero (y)
#if 0
  Umu=Zero();
  Complex cone(1.0,0.0);
  for(int nn=0;nn<Nd;nn++){
    random(pRNG,U[nn]);
    if(1) {
      if (nn!=2) { U[nn]=Zero(); std::cout<<GridLogMessage << "zeroing gauge field in dir "<<nn<<std::endl; }
      //      else       { U[nn]= cone;std::cout<<GridLogMessage << "unit gauge field in dir "<<nn<<std::endl; }
      else       { std::cout<<GridLogMessage << "random gauge field in dir "<<nn<<std::endl; }
    }
    PokeIndex<LorentzIndex>(Umu,U[nn],nn);
  }
#endif

  for(int mu=0;mu<Nd;mu++){
    U[mu] = PeekIndex<LorentzIndex>(Umu,mu);
  }

  { // Naive wilson implementation
    ref = Zero();
    for(int mu=0;mu<Nd;mu++){
      //    ref =  src + Gamma(Gamma::Algebra::GammaX)* src ; // 1-gamma_x
      tmp = U[mu]*Cshift(src,mu,1);
      {
	autoView( ref_v, ref, CpuWrite);
	autoView( tmp_v, tmp, CpuWrite);
	for(int i=0;i<ref_v.size();i++){
	  ref_v[i]+= tmp_v[i] - Gamma(Gmu[mu])*tmp_v[i]; ;
	}
      }

      tmp =adj(U[mu])*src;
      tmp =Cshift(tmp,mu,-1);
      {
	autoView( ref_v, ref, CpuWrite);
	autoView( tmp_v, tmp, CpuWrite);
	for(int i=0;i<ref_v.size();i++){
	  ref_v[i]+= tmp_v[i] + Gamma(Gmu[mu])*tmp_v[i]; ;
	}
      }
    }
  }
  ref = -0.5*ref;
  RealD mass=0.1;

  typename WilsonFermionR::ImplParams params;

  WilsonFermionR Dw(Umu,Grid,RBGrid,mass,params);

  std::cout<<GridLogMessage << "Calling Dw"<<std::endl;
  int ncall=1000;
  //int ncall=1;

  // Counters
  Dw.ZeroCounters();
  Grid.Barrier();

  double t0=usecond();
  for(int i=0;i<ncall;i++){
    Dw.Dhop(src,result,0);
  }

  // Counters
  Grid.Barrier();

  double t1=usecond();
  double flops=single_site_flops*volume*ncall;

  if (perfProfiling){
  std::cout<<GridLogMessage << "Profiling Dw with perf"<<std::endl;

  System::profile("kernel", [&]() {
    for(int i=0;i<ncall;i++){
      Dw.Dhop(src,result,0);
    }
  });

  std::cout<<GridLogMessage << "Generated kernel.data"<<std::endl;
  std::cout<<GridLogMessage << "Use with: perf report -i kernel.data"<<std::endl;

  }

  auto nsimd = vComplex::Nsimd();
  auto simdwidth = sizeof(vComplex);

  std::cout<<GridLogMessage << "Nsimd "<< nsimd << std::endl;
  std::cout<<GridLogMessage << "Simd width "<< simdwidth << std::endl;

  // RF: Nd Wilson, Nd gauge, Nc colors
  double data = volume * ((2*Nd+1)*Nd*Nc + 2*Nd*Nc*Nc) * simdwidth / nsimd * ncall / (1024.*1024.*1024.);

  std::cout<<GridLogMessage << "Called Dw"<<std::endl;
  std::cout<<GridLogMessage << "flops per site " << single_site_flops << std::endl;
  std::cout<<GridLogMessage << "norm result "<< norm2(result)<<std::endl;
  std::cout<<GridLogMessage << "norm ref    "<< norm2(ref)<<std::endl;
  std::cout<<GridLogMessage << "mflop/s =   "<< flops/(t1-t0)<<std::endl;
  std::cout<<GridLogMessage << "RF  GiB/s (base 2) =   "<< 1000000. * data/(t1-t0)<<std::endl;
  err = ref-result;
  std::cout<<GridLogMessage << "norm diff   "<< norm2(err)<<std::endl;

  Dw.Report();
  
  // guard
  double err0 = norm2(err);

  //  for(int ss=0;ss<10;ss++ ){
  for(int ss=0;ss<0;ss++ ){
    for(int i=0;i<Ns;i++){
      for(int j=0;j<Nc;j++){
	autoView( ref_v, ref, CpuWrite);
	autoView( result_v, result, CpuWrite);
	ComplexF * ref_p = (ComplexF *)&ref_v[ss]()(i)(j);
	ComplexF * res_p = (ComplexF *)&result_v[ss]()(i)(j);
	std::cout<<GridLogMessage << ss<< " "<<i<<" "<<j<<" "<< (*ref_p)<<" " <<(*res_p)<<std::endl;
      }
    }
  }

  { // Naive wilson dag implementation
    ref = Zero();
    for(int mu=0;mu<Nd;mu++){


      //    ref =  src - Gamma(Gamma::Algebra::GammaX)* src ; // 1+gamma_x
      tmp = U[mu]*Cshift(src,mu,1);
      {
	autoView( ref_v, ref, CpuWrite);
	autoView( tmp_v, tmp, CpuWrite);
	for(int i=0;i<ref_v.size();i++){
	  ref_v[i]+= tmp_v[i] + Gamma(Gmu[mu])*tmp_v[i]; ;
	}
      }

      tmp =adj(U[mu])*src;
      tmp =Cshift(tmp,mu,-1);
      {
	autoView( ref_v, ref, CpuWrite);
	autoView( tmp_v, tmp, CpuWrite);
	for(int i=0;i<ref_v.size();i++){
	  ref_v[i]+= tmp_v[i] - Gamma(Gmu[mu])*tmp_v[i]; ;
	}
      }
    }
  }
  ref = -0.5*ref;
  Dw.Dhop(src,result,1);
  std::cout<<GridLogMessage << "Called DwDag"<<std::endl;
  std::cout<<GridLogMessage << "norm result "<< norm2(result)<<std::endl;
  std::cout<<GridLogMessage << "norm ref    "<< norm2(ref)<<std::endl;
  err = ref-result;
  std::cout<<GridLogMessage << "norm diff   "<< norm2(err)<<std::endl;

  // guard
  double err1 = norm2(err);
  assert(fabs(err0) < 1.0e-3);
  assert(fabs(err1) < 1.0e-3);

  Grid_finalize();
}
