    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./benchmarks/Benchmark_memory_bandwidth.cc

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

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  const int Nvec=8;
  typedef Lattice< iVector< vReal,Nvec> > LatticeVec;
  typedef iVector<vReal,Nvec> Vec;


  Vec rn = zero;

  std::vector<int> simd_layout = GridDefaultSimd(Nd,vReal::Nsimd());
  std::vector<int> mpi_layout  = GridDefaultMpi();

  int threads = GridThread::GetThreads();
  std::cout<<GridLogMessage << "Grid is setup to use "<<threads<<" threads"<<std::endl;

  
  std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
  std::cout<<GridLogMessage << "= Benchmarking fused AXPY bandwidth ; sizeof(Real) "<<sizeof(Real)<<std::endl;
  std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
  std::cout<<GridLogMessage << "  L  "<<"\t\t"<<"bytes"<<"\t\t\t"<<"GB/s"<<"\t\t"<<"Gflop/s"<<"\t\t seconds"<<std::endl;
  std::cout<<GridLogMessage << "----------------------------------------------------------"<<std::endl;
  uint64_t lmax=44;
#define NLOOP (1*lmax*lmax*lmax*lmax/vol)
  for(int lat=4;lat<=lmax;lat+=4){

      std::vector<int> latt_size  ({lat*mpi_layout[0],lat*mpi_layout[1],lat*mpi_layout[2],lat*mpi_layout[3]});
      int vol = latt_size[0]*latt_size[1]*latt_size[2]*latt_size[3];
      GridCartesian     Grid(latt_size,simd_layout,mpi_layout);

      uint64_t Nloop=NLOOP;

      //      GridParallelRNG          pRNG(&Grid);      pRNG.SeedRandomDevice();

      LatticeVec z(&Grid); //random(pRNG,z);
      LatticeVec x(&Grid); //random(pRNG,x);
      LatticeVec y(&Grid); //random(pRNG,y);
      double a=2.0;


      double start=usecond();
      for(int i=0;i<Nloop;i++){
	axpy(z,a,x,y);
        x._odata[0]=z._odata[0]; // serial loop dependence to prevent optimise
        y._odata[4]=z._odata[4];
      }
      double stop=usecond();
      double time = (stop-start)/Nloop*1000;
      
      double flops=vol*Nvec*2;// mul,add
      double bytes=3*vol*Nvec*sizeof(Real);
      std::cout<<GridLogMessage<<std::setprecision(3) << lat<<"\t\t"<<bytes<<"   \t\t"<<bytes/time<<"\t\t"<<flops/time<<"\t\t"<<(stop-start)/1000./1000.<<std::endl;

    }

  std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
  std::cout<<GridLogMessage << "= Benchmarking a*x + y bandwidth"<<std::endl;
  std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
  std::cout<<GridLogMessage << "  L  "<<"\t\t"<<"bytes"<<"\t\t\t"<<"GB/s"<<"\t\t"<<"Gflop/s"<<"\t\t seconds"<<std::endl;
  std::cout<<GridLogMessage << "----------------------------------------------------------"<<std::endl;
  
  for(int lat=4;lat<=lmax;lat+=4){

      std::vector<int> latt_size  ({lat*mpi_layout[0],lat*mpi_layout[1],lat*mpi_layout[2],lat*mpi_layout[3]});
      int vol = latt_size[0]*latt_size[1]*latt_size[2]*latt_size[3];
      GridCartesian     Grid(latt_size,simd_layout,mpi_layout);

      //      GridParallelRNG          pRNG(&Grid);      pRNG.SeedRandomDevice();

      LatticeVec z(&Grid); //random(pRNG,z);
      LatticeVec x(&Grid); //random(pRNG,x);
      LatticeVec y(&Grid); //random(pRNG,y);
      double a=2.0;

      uint64_t Nloop=NLOOP;

      double start=usecond();
      for(int i=0;i<Nloop;i++){
	z=a*x-y;
        x._odata[0]=z._odata[0]; // force serial dependency to prevent optimise away
        y._odata[4]=z._odata[4];
      }
      double stop=usecond();
      double time = (stop-start)/Nloop*1000;
     
      double flops=vol*Nvec*2;// mul,add
      double bytes=3*vol*Nvec*sizeof(Real);
      std::cout<<GridLogMessage<<std::setprecision(3) << lat<<"\t\t"<<bytes<<"   \t\t"<<bytes/time<<"\t\t"<<flops/time<<"\t\t"<<(stop-start)/1000./1000.<<std::endl;

    }

  std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
  std::cout<<GridLogMessage << "= Benchmarking SCALE bandwidth"<<std::endl;
  std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
  std::cout<<GridLogMessage << "  L  "<<"\t\t"<<"bytes"<<"\t\t\t"<<"GB/s"<<"\t\t"<<"Gflop/s"<<"\t\t seconds"<<std::endl;

  for(int lat=4;lat<=lmax;lat+=4){


      std::vector<int> latt_size  ({lat*mpi_layout[0],lat*mpi_layout[1],lat*mpi_layout[2],lat*mpi_layout[3]});
      int vol = latt_size[0]*latt_size[1]*latt_size[2]*latt_size[3];
      uint64_t Nloop=NLOOP;

      GridCartesian     Grid(latt_size,simd_layout,mpi_layout);

      //      GridParallelRNG          pRNG(&Grid);      pRNG.SeedRandomDevice();

      LatticeVec z(&Grid); //random(pRNG,z);
      LatticeVec x(&Grid); //random(pRNG,x);
      LatticeVec y(&Grid); //random(pRNG,y);
      RealD a=2.0;


      double start=usecond();
      for(int i=0;i<Nloop;i++){
	z=a*x;
        x._odata[0]=z._odata[0]*2.0;
      }
      double stop=usecond();
      double time = (stop-start)/Nloop*1000;
      
      double bytes=2*vol*Nvec*sizeof(Real);
      double flops=vol*Nvec*1;// mul
      std::cout<<GridLogMessage <<std::setprecision(3) << lat<<"\t\t"<<bytes<<"   \t\t"<<bytes/time<<"\t\t"<<flops/time<<"\t\t"<<(stop-start)/1000./1000.<<std::endl;

  }

  std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
  std::cout<<GridLogMessage << "= Benchmarking READ bandwidth"<<std::endl;
  std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
  std::cout<<GridLogMessage << "  L  "<<"\t\t"<<"bytes"<<"\t\t\t"<<"GB/s"<<"\t\t"<<"Gflop/s"<<"\t\t seconds"<<std::endl;
  std::cout<<GridLogMessage << "----------------------------------------------------------"<<std::endl;

  for(int lat=4;lat<=lmax;lat+=4){

      std::vector<int> latt_size  ({lat*mpi_layout[0],lat*mpi_layout[1],lat*mpi_layout[2],lat*mpi_layout[3]});
      int vol = latt_size[0]*latt_size[1]*latt_size[2]*latt_size[3];
      uint64_t Nloop=NLOOP;
      GridCartesian     Grid(latt_size,simd_layout,mpi_layout);

      //      GridParallelRNG          pRNG(&Grid);      pRNG.SeedRandomDevice();
      LatticeVec z(&Grid); //random(pRNG,z);
      LatticeVec x(&Grid); //random(pRNG,x);
      LatticeVec y(&Grid); //random(pRNG,y);
      RealD a=2.0;
      Real nn;      
      double start=usecond();
      for(int i=0;i<Nloop;i++){
	nn=norm2(x);
	vsplat(x._odata[0]._internal[0],nn);
      }
      double stop=usecond();
      double time = (stop-start)/Nloop*1000;
      
      double bytes=vol*Nvec*sizeof(Real);
      double flops=vol*Nvec*2;// mul,add
      std::cout<<GridLogMessage<<std::setprecision(3) << lat<<"\t\t"<<bytes<<"  \t\t"<<bytes/time<<"\t\t"<<flops/time<< "\t\t"<<(stop-start)/1000./1000.<< "\t\t " <<std::endl;

  }    

  Grid_finalize();
}
