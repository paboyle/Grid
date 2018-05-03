    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./benchmarks/Benchmark_su3.cc

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: Peter Boyle <peterboyle@Peters-MacBook-Pro-2.local>

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
#define LMAX (32)
#define LMIN (16)
#define LINC (4)

  int64_t Nloop=2000;

  std::vector<int> simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
  std::vector<int> mpi_layout  = GridDefaultMpi();

  int64_t threads = GridThread::GetThreads();
  std::cout<<GridLogMessage << "Grid is setup to use "<<threads<<" threads"<<std::endl;

  std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
  std::cout<<GridLogMessage << "= Benchmarking SU3xSU3  x= x*y"<<std::endl;
  std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
  std::cout<<GridLogMessage << "  L  "<<"\t\t"<<"bytes"<<"\t\t\t"<<"GB/s\t\t GFlop/s"<<std::endl;
  std::cout<<GridLogMessage << "----------------------------------------------------------"<<std::endl;

  for(int lat=LMIN;lat<=LMAX;lat+=LINC){

      std::vector<int> latt_size  ({lat*mpi_layout[0],lat*mpi_layout[1],lat*mpi_layout[2],lat*mpi_layout[3]});
      int64_t vol = latt_size[0]*latt_size[1]*latt_size[2]*latt_size[3];
      GridCartesian     Grid(latt_size,simd_layout,mpi_layout);
      GridParallelRNG          pRNG(&Grid);      pRNG.SeedFixedIntegers(std::vector<int>({45,12,81,9}));

      LatticeColourMatrix z(&Grid); random(pRNG,z);
      LatticeColourMatrix x(&Grid); random(pRNG,x);
      LatticeColourMatrix y(&Grid); random(pRNG,y);

      double start=usecond();
      for(int64_t i=0;i<Nloop;i++){
	x=x*y;
      }
      double stop=usecond();
      double time = (stop-start)/Nloop*1000.0;
      
      double bytes=3.0*vol*Nc*Nc*sizeof(Complex);
      double footprint=2.0*vol*Nc*Nc*sizeof(Complex);
      double flops=Nc*Nc*(6.0+8.0+8.0)*vol;
      std::cout<<GridLogMessage<<std::setprecision(3) << lat<<"\t\t"<<footprint<<"    \t\t"<<bytes/time<<"\t\t" << flops/time<<std::endl;

    }


  std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
  std::cout<<GridLogMessage << "= Benchmarking SU3xSU3  z= x*y"<<std::endl;
  std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
  std::cout<<GridLogMessage << "  L  "<<"\t\t"<<"bytes"<<"\t\t\t"<<"GB/s\t\t GFlop/s"<<std::endl;
  std::cout<<GridLogMessage << "----------------------------------------------------------"<<std::endl;

  for(int lat=LMIN;lat<=LMAX;lat+=LINC){

      std::vector<int> latt_size  ({lat*mpi_layout[0],lat*mpi_layout[1],lat*mpi_layout[2],lat*mpi_layout[3]});
      int64_t vol = latt_size[0]*latt_size[1]*latt_size[2]*latt_size[3];

      GridCartesian     Grid(latt_size,simd_layout,mpi_layout);
      GridParallelRNG          pRNG(&Grid);      pRNG.SeedFixedIntegers(std::vector<int>({45,12,81,9}));

      LatticeColourMatrix z(&Grid); random(pRNG,z);
      LatticeColourMatrix x(&Grid); random(pRNG,x);
      LatticeColourMatrix y(&Grid); random(pRNG,y);

      double start=usecond();
      for(int64_t i=0;i<Nloop;i++){
	z=x*y;
      }
      double stop=usecond();
      double time = (stop-start)/Nloop*1000.0;
      
      double bytes=3*vol*Nc*Nc*sizeof(Complex);
      double flops=Nc*Nc*(6+8+8)*vol;
      std::cout<<GridLogMessage<<std::setprecision(3) << lat<<"\t\t"<<bytes<<"    \t\t"<<bytes/time<<"\t\t" << flops/time<<std::endl;

    }

  std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
  std::cout<<GridLogMessage << "= Benchmarking SU3xSU3  mult(z,x,y)"<<std::endl;
  std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
  std::cout<<GridLogMessage << "  L  "<<"\t\t"<<"bytes"<<"\t\t\t"<<"GB/s\t\t GFlop/s"<<std::endl;
  std::cout<<GridLogMessage << "----------------------------------------------------------"<<std::endl;

  for(int lat=LMIN;lat<=LMAX;lat+=LINC){

      std::vector<int> latt_size  ({lat*mpi_layout[0],lat*mpi_layout[1],lat*mpi_layout[2],lat*mpi_layout[3]});
      int64_t vol = latt_size[0]*latt_size[1]*latt_size[2]*latt_size[3];

      GridCartesian     Grid(latt_size,simd_layout,mpi_layout);
      GridParallelRNG          pRNG(&Grid);      pRNG.SeedFixedIntegers(std::vector<int>({45,12,81,9}));

      LatticeColourMatrix z(&Grid); random(pRNG,z);
      LatticeColourMatrix x(&Grid); random(pRNG,x);
      LatticeColourMatrix y(&Grid); random(pRNG,y);

      double start=usecond();
      for(int64_t i=0;i<Nloop;i++){
	mult(z,x,y);
      }
      double stop=usecond();
      double time = (stop-start)/Nloop*1000.0;
      
      double bytes=3*vol*Nc*Nc*sizeof(Complex);
      double flops=Nc*Nc*(6+8+8)*vol;
      std::cout<<GridLogMessage<<std::setprecision(3) << lat<<"\t\t"<<bytes<<"    \t\t"<<bytes/time<<"\t\t" << flops/time<<std::endl;

    }

  std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
  std::cout<<GridLogMessage << "= Benchmarking SU3xSU3  mac(z,x,y)"<<std::endl;
  std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
  std::cout<<GridLogMessage << "  L  "<<"\t\t"<<"bytes"<<"\t\t\t"<<"GB/s\t\t GFlop/s"<<std::endl;
  std::cout<<GridLogMessage << "----------------------------------------------------------"<<std::endl;

  for(int lat=LMIN;lat<=LMAX;lat+=LINC){
    
    std::vector<int> latt_size  ({lat*mpi_layout[0],lat*mpi_layout[1],lat*mpi_layout[2],lat*mpi_layout[3]});
    int64_t vol = latt_size[0]*latt_size[1]*latt_size[2]*latt_size[3];
    
    GridCartesian     Grid(latt_size,simd_layout,mpi_layout);
    GridParallelRNG          pRNG(&Grid);      pRNG.SeedFixedIntegers(std::vector<int>({45,12,81,9}));
    
    LatticeColourMatrix z(&Grid); random(pRNG,z);
    LatticeColourMatrix x(&Grid); random(pRNG,x);
    LatticeColourMatrix y(&Grid); random(pRNG,y);
    
    double start=usecond();
    for(int64_t i=0;i<Nloop;i++){
      mac(z,x,y);
    }
    double stop=usecond();
    double time = (stop-start)/Nloop*1000.0;
    
    double bytes=3*vol*Nc*Nc*sizeof(Complex);
    double flops=Nc*Nc*(6+8+8)*vol;
    std::cout<<GridLogMessage<<std::setprecision(3) << lat<<"\t\t"<<bytes<<"   \t\t"<<bytes/time<<"\t\t" << flops/time<<std::endl;
    
  }

  std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
  std::cout<<GridLogMessage << "= Benchmarking SU3xSU3  CovShiftForward(z,x,y)"<<std::endl;
  std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
  std::cout<<GridLogMessage << "  L  "<<"\t\t"<<"bytes"<<"\t\t\t"<<"GB/s\t\t GFlop/s"<<std::endl;
  std::cout<<GridLogMessage << "----------------------------------------------------------"<<std::endl;

  for(int lat=LMIN;lat<=LMAX;lat+=LINC){

      std::vector<int> latt_size  ({lat*mpi_layout[0],lat*mpi_layout[1],lat*mpi_layout[2],lat*mpi_layout[3]});
      int64_t vol = latt_size[0]*latt_size[1]*latt_size[2]*latt_size[3];

      GridCartesian     Grid(latt_size,simd_layout,mpi_layout);
      GridParallelRNG          pRNG(&Grid);      pRNG.SeedFixedIntegers(std::vector<int>({45,12,81,9}));

      LatticeColourMatrix z(&Grid); random(pRNG,z);
      LatticeColourMatrix x(&Grid); random(pRNG,x);
      LatticeColourMatrix y(&Grid); random(pRNG,y);

      for(int mu=0;mu<4;mu++){
	      double start=usecond();
	      for(int64_t i=0;i<Nloop;i++){
	        z = PeriodicBC::CovShiftForward(x,mu,y);
	    }
	    double stop=usecond();
	    double time = (stop-start)/Nloop*1000.0;
	
	
	    double bytes=3*vol*Nc*Nc*sizeof(Complex);
	    double flops=Nc*Nc*(6+8+8)*vol;
	    std::cout<<GridLogMessage<<std::setprecision(3) << lat<<"\t\t"<<bytes<<"   \t\t"<<bytes/time<<"\t\t" << flops/time<<std::endl;
      }
  }
#if 1
  std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
  std::cout<<GridLogMessage << "= Benchmarking SU3xSU3  z= x * Cshift(y)"<<std::endl;
  std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
  std::cout<<GridLogMessage << "  L  "<<"\t\t"<<"bytes"<<"\t\t\t"<<"GB/s\t\t GFlop/s"<<std::endl;
  std::cout<<GridLogMessage << "----------------------------------------------------------"<<std::endl;

  for(int lat=LMIN;lat<=LMAX;lat+=LINC){
      std::vector<int> latt_size  ({lat*mpi_layout[0],lat*mpi_layout[1],lat*mpi_layout[2],lat*mpi_layout[3]});
      int64_t vol = latt_size[0]*latt_size[1]*latt_size[2]*latt_size[3];

      GridCartesian     Grid(latt_size,simd_layout,mpi_layout);
      GridParallelRNG          pRNG(&Grid);      pRNG.SeedFixedIntegers(std::vector<int>({45,12,81,9}));

      LatticeColourMatrix z(&Grid); random(pRNG,z);
      LatticeColourMatrix x(&Grid); random(pRNG,x);
      LatticeColourMatrix y(&Grid); random(pRNG,y);
      LatticeColourMatrix tmp(&Grid);

      for(int mu=0;mu<4;mu++){
	double tshift=0;
	double tmult =0;

	double start=usecond();
	for(int64_t i=0;i<Nloop;i++){
	  tshift-=usecond();
	  tmp = Cshift(y,mu,-1);
	  tshift+=usecond();
	  tmult-=usecond();
	  z   = x*tmp;
	  tmult+=usecond();
	}
	double stop=usecond();
	double time = (stop-start)/Nloop;
	tshift = tshift/Nloop;
	tmult  = tmult /Nloop;
	
	double bytes=3*vol*Nc*Nc*sizeof(Complex);
	double flops=Nc*Nc*(6+8+8)*vol;
	std::cout<<GridLogMessage<<std::setprecision(3) << "total us "<<time<<" shift "<<tshift <<" mult "<<tmult<<std::endl;
	time = time * 1000; // convert to NS for GB/s
	std::cout<<GridLogMessage<<std::setprecision(3) << lat<<"\t\t"<<bytes<<"   \t\t"<<bytes/time<<"\t\t" << flops/time<<std::endl;
      }
    }
#endif
  Grid_finalize();
}
