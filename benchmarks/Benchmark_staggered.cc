    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./benchmarks/Benchmark_staggered.cc

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

  std::vector<int> latt_size   = GridDefaultLatt();
  std::vector<int> simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
  std::vector<int> mpi_layout  = GridDefaultMpi();
  GridCartesian               Grid(latt_size,simd_layout,mpi_layout);
  GridRedBlackCartesian     RBGrid(&Grid);

  int threads = GridThread::GetThreads();
  std::cout<<GridLogMessage << "Grid is setup to use "<<threads<<" threads"<<std::endl;
  std::cout<<GridLogMessage << "Grid floating point word size is REALF"<< sizeof(RealF)<<std::endl;
  std::cout<<GridLogMessage << "Grid floating point word size is REALD"<< sizeof(RealD)<<std::endl;
  std::cout<<GridLogMessage << "Grid floating point word size is REAL"<< sizeof(Real)<<std::endl;

  std::vector<int> seeds({1,2,3,4});
  GridParallelRNG          pRNG(&Grid);
  pRNG.SeedFixedIntegers(seeds);
  //  pRNG.SeedFixedIntegers(std::vector<int>({45,12,81,9});

  typedef typename ImprovedStaggeredFermionR::FermionField FermionField; 
  typename ImprovedStaggeredFermionR::ImplParams params; 

  FermionField src   (&Grid); random(pRNG,src);
  FermionField result(&Grid); result=zero;
  FermionField    ref(&Grid);    ref=zero;
  FermionField    tmp(&Grid);    tmp=zero;
  FermionField    err(&Grid);    tmp=zero;
  LatticeGaugeField Umu(&Grid); random(pRNG,Umu);
  std::vector<LatticeColourMatrix> U(4,&Grid);

  double volume=1;
  for(int mu=0;mu<Nd;mu++){
    volume=volume*latt_size[mu];
  }  

  // Only one non-zero (y)
#if 0
  Umu=zero;
  Complex cone(1.0,0.0);
  for(int nn=0;nn<Nd;nn++){
    random(pRNG,U[nn]);
    if(1) {
      if (nn!=2) { U[nn]=zero; std::cout<<GridLogMessage << "zeroing gauge field in dir "<<nn<<std::endl; }
      //      else       { U[nn]= cone;std::cout<<GridLogMessage << "unit gauge field in dir "<<nn<<std::endl; }
      else       { std::cout<<GridLogMessage << "random gauge field in dir "<<nn<<std::endl; }
    }
    PokeIndex<LorentzIndex>(Umu,U[nn],nn);
  }
#endif

  for(int mu=0;mu<Nd;mu++){
    U[mu] = PeekIndex<LorentzIndex>(Umu,mu);
  }
  ref = zero;
  /*  
  { // Naive wilson implementation
    ref = zero;
    for(int mu=0;mu<Nd;mu++){
      //    ref =  src + Gamma(Gamma::GammaX)* src ; // 1-gamma_x
      tmp = U[mu]*Cshift(src,mu,1);
      for(int i=0;i<ref._odata.size();i++){
	ref._odata[i]+= tmp._odata[i] - Gamma(Gmu[mu])*tmp._odata[i]; ;
      }

      tmp =adj(U[mu])*src;
      tmp =Cshift(tmp,mu,-1);
      for(int i=0;i<ref._odata.size();i++){
	ref._odata[i]+= tmp._odata[i] + Gamma(Gmu[mu])*tmp._odata[i]; ;
      }
    }
  }
  ref = -0.5*ref;
  */

  RealD mass=0.1;
  RealD c1=9.0/8.0;
  RealD c2=-1.0/24.0;
  RealD u0=1.0;
  ImprovedStaggeredFermionR Ds(Umu,Umu,Grid,RBGrid,mass,c1,c2,u0,params);
  
  std::cout<<GridLogMessage << "Calling Ds"<<std::endl;
  int ncall=1000;
  double t0=usecond();
  for(int i=0;i<ncall;i++){
    Ds.Dhop(src,result,0);
  }
  double t1=usecond();
  double flops=(16*(3*(6+8+8)) + 15*3*2)*volume*ncall; // == 66*16 +  == 1146
  
  std::cout<<GridLogMessage << "Called Ds"<<std::endl;
  std::cout<<GridLogMessage << "norm result "<< norm2(result)<<std::endl;
  std::cout<<GridLogMessage << "norm ref    "<< norm2(ref)<<std::endl;
  std::cout<<GridLogMessage << "mflop/s =   "<< flops/(t1-t0)<<std::endl;
  err = ref-result; 
  std::cout<<GridLogMessage << "norm diff   "<< norm2(err)<<std::endl;

  Grid_finalize();
}
