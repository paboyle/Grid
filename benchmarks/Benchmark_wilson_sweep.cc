/*************************************************************************************
    Grid physics library, www.github.com/paboyle/Grid 
    Source file: ./benchmarks/Benchmark_wilson.cc
    Copyright (C) 2015
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: paboyle <paboyle@ph.ed.ac.uk>
Author: Richard Rollins <rprollins@users.noreply.github.com>
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

void bench_wilson (
		   LatticeFermion &    src,
		   LatticeFermion & result,
		   WilsonFermionR &     Dw,
		   double const     volume,
		   int const           dag );

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);
  typename WilsonFermionR::ImplParams params;

  Coordinate simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
  Coordinate mpi_layout  = GridDefaultMpi();
  std::vector<int> seeds({1,2,3,4});
  RealD mass = 0.1;

  std::cout << GridLogMessage<< "*****************************************************************" <<std::endl;
  std::cout << GridLogMessage<< "* Kernel options --dslash-generic, --dslash-unroll, --dslash-asm" <<std::endl;
  std::cout << GridLogMessage<< "*****************************************************************" <<std::endl;
  std::cout << GridLogMessage<< "*****************************************************************" <<std::endl;
  std::cout << GridLogMessage<< "* Number of colours "<< Nc <<std::endl;
  std::cout << GridLogMessage<< "* Benchmarking WilsonFermionR::Dhop                  "<<std::endl;
  std::cout << GridLogMessage<< "* Vectorising space-time by "<<vComplex::Nsimd()<<std::endl;
  if ( sizeof(Real)==4 )   std::cout << GridLogMessage<< "* SINGLE precision "<<std::endl;
  if ( sizeof(Real)==8 )   std::cout << GridLogMessage<< "* DOUBLE precision "<<std::endl;
  if ( WilsonKernelsStatic::Opt == WilsonKernelsStatic::OptGeneric   ) std::cout << GridLogMessage<< "* Using GENERIC Nc WilsonKernels" <<std::endl;
  if ( WilsonKernelsStatic::Opt == WilsonKernelsStatic::OptHandUnroll) std::cout << GridLogMessage<< "* Using Nc=3       WilsonKernels" <<std::endl;
  if ( WilsonKernelsStatic::Opt == WilsonKernelsStatic::OptInlineAsm ) std::cout << GridLogMessage<< "* Using Asm Nc=3   WilsonKernels" <<std::endl;
  std::cout << GridLogMessage << "* OpenMP threads       : "<< GridThread::GetThreads() <<std::endl;
  std::cout << GridLogMessage << "* MPI tasks            : "<< GridCmdVectorIntToString(mpi_layout) << std::endl;
  std::cout << GridLogMessage<< "*****************************************************************" <<std::endl;

  std::cout<<GridLogMessage << "================================================================================================="<< std::endl;
  std::cout<<GridLogMessage << "= Benchmarking Wilson operator in the fundamental representation" << std::endl;
  std::cout<<GridLogMessage << "================================================================================================="<< std::endl;
  std::cout<<GridLogMessage << "Volume\t\t\tWilson/MFLOPs\tWilsonDag/MFLOPs\tWilsonEO/MFLOPs\tWilsonDagEO/MFLOPs" << std::endl;
  std::cout<<GridLogMessage << "================================================================================================="<< std::endl;

  int Lmax = 32;
  int dmin = 0;
  if ( getenv("LMAX") ) Lmax=atoi(getenv("LMAX"));
  if ( getenv("DMIN") ) dmin=atoi(getenv("DMIN"));
  for (int L=8; L<=Lmax; L*=2)
    {
      Coordinate latt_size = Coordinate(4,L);
      for(int d=4; d>dmin; d--)
	{
	  if ( d<=3 ) { latt_size[d] *= 2; }

	  std::cout << GridLogMessage;
	  std::cout << latt_size;

	  GridCartesian           Grid(latt_size,simd_layout,mpi_layout);
	  GridRedBlackCartesian RBGrid(&Grid);

	  GridParallelRNG  pRNG(&Grid); pRNG.SeedFixedIntegers(seeds);
	  LatticeGaugeField Umu(&Grid); random(pRNG,Umu);
	  LatticeFermion    src(&Grid); random(pRNG,src);
	  LatticeFermion    src_o(&RBGrid); pickCheckerboard(Odd,src_o,src);
	  LatticeFermion     result(&Grid); result=Zero();
	  LatticeFermion result_e(&RBGrid); result_e=Zero();

	  double volume = std::accumulate(latt_size.begin(),latt_size.end(),1,std::multiplies<int>());

	  WilsonFermionR Dw(Umu,Grid,RBGrid,mass,params);
      
    // Full operator      
	  bench_wilson(src,result,Dw,volume,DaggerNo);
	  bench_wilson(src,result,Dw,volume,DaggerYes);
	  std::cout << "\t";
    // EO
	  bench_wilson(src,result,Dw,volume,DaggerNo);
	  bench_wilson(src,result,Dw,volume,DaggerYes);
	  std::cout << std::endl;
	}
    }

  std::cout<<GridLogMessage << "============================================================================="<< std::endl;
  Grid_finalize();
}

void bench_wilson (
		   LatticeFermion &    src,
		   LatticeFermion & result,
		   WilsonFermionR &     Dw,
		   double const     volume,
		   int const           dag )
{
  int ncall    = 1000;
  long unsigned int single_site_flops = 8*Nc*(7+16*Nc);
  double t0    = usecond();
  for(int i=0; i<ncall; i++) { Dw.Dhop(src,result,dag); }
  double t1    = usecond();
  double flops = single_site_flops * volume * ncall;
  std::cout << flops/(t1-t0) << "\t\t";
}

void bench_wilson_eo (
		   LatticeFermion &    src,
		   LatticeFermion & result,
		   WilsonFermionR &     Dw,
		   double const     volume,
		   int const           dag )
{
  int ncall    = 1000;
  long unsigned int single_site_flops = 8*Nc*(7+16*Nc);
  double t0    = usecond();
  for(int i=0; i<ncall; i++) { Dw.DhopEO(src,result,dag); }
  double t1    = usecond();
  double flops = (single_site_flops * volume * ncall)/2.0;
  std::cout << flops/(t1-t0) << "\t\t";
}
