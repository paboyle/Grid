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


int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

#define LMIN (16)
#define LMAX (48)
#define LADD (8)

  int64_t Nwarm=20;
  int64_t Nloop=500;

  Coordinate simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
  Coordinate mpi_layout  = GridDefaultMpi();

  int64_t threads = GridThread::GetThreads();
  std::cout<<GridLogMessage << "Grid is setup to use "<<threads<<" threads"<<std::endl;

  std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
  std::cout<<GridLogMessage << "= Benchmarking SU3xSU3  ext/ins z = x*y"<<std::endl;
  std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
  std::cout<<GridLogMessage << "  L  "<<"\t\t"<<"bytes"<<"\t\t\t"<<"GB/s\t\t GFlop/s"<<std::endl;
  std::cout<<GridLogMessage << "----------------------------------------------------------"<<std::endl;

  for(int lat=LMIN;lat<=LMAX;lat+=LADD){

      Coordinate latt_size  ({lat*mpi_layout[0],lat*mpi_layout[1],lat*mpi_layout[2],lat*mpi_layout[3]});
      int64_t vol = latt_size[0]*latt_size[1]*latt_size[2]*latt_size[3];
      GridCartesian     Grid(latt_size,simd_layout,mpi_layout);
      GridParallelRNG          pRNG(&Grid);      pRNG.SeedFixedIntegers(std::vector<int>({45,12,81,9}));

      LatticeColourMatrix z(&Grid);// random(pRNG,z);
      LatticeColourMatrix x(&Grid);// random(pRNG,x);
      LatticeColourMatrix y(&Grid);// random(pRNG,y);

      auto x_v = x.View();
      auto y_v = y.View();
      auto z_v = z.View();
      const uint64_t Nsite = x_v.size();
      const uint64_t nsimd = vComplex::Nsimd();
      const uint64_t    NN = Nsite*nsimd;
      for(int64_t i=0;i<Nwarm;i++){
	accelerator_loopN( sss, NN, {
	  uint64_t lane = sss % nsimd; 
	  uint64_t ss   = sss / nsimd;
	  auto xx = extractLane(lane,x_v[ss]);
	  auto yy = extractLane(lane,y_v[ss]);
	  auto zz = xx*yy;
	  insertLane(lane,z_v[ss],zz);
        });
      }
      double start=usecond();
      for(int64_t i=0;i<Nloop;i++){
	accelerator_loopN( sss, NN, {
	  uint64_t lane = sss % nsimd; 
	  uint64_t ss   = sss / nsimd;
	  auto xx = extractLane(lane,x_v[ss]);
	  auto yy = extractLane(lane,y_v[ss]);
	  auto zz = xx*yy;
	  insertLane(lane,z_v[ss],zz);
        });
      }
      double stop=usecond();
      double time = (stop-start)/Nloop*1000.0;
      
      double bytes=3.0*vol*Nc*Nc*sizeof(Complex);
      double footprint=3.0*vol*Nc*Nc*sizeof(Complex);
      double flops=Nc*Nc*(6.0+8.0+8.0)*vol;
      std::cout<<GridLogMessage<<std::setprecision(3) << lat<<"\t\t"<<footprint<<"    \t\t"<<bytes/time<<"\t\t" << flops/time<<std::endl;

    }


  std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
  std::cout<<GridLogMessage << "= Benchmarking SU3xSU3  ext/ins x = x*y"<<std::endl;
  std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
  std::cout<<GridLogMessage << "  L  "<<"\t\t"<<"bytes"<<"\t\t\t"<<"GB/s\t\t GFlop/s"<<std::endl;
  std::cout<<GridLogMessage << "----------------------------------------------------------"<<std::endl;

  for(int lat=LMIN;lat<=LMAX;lat+=LADD){

      Coordinate latt_size  ({lat*mpi_layout[0],lat*mpi_layout[1],lat*mpi_layout[2],lat*mpi_layout[3]});
      int64_t vol = latt_size[0]*latt_size[1]*latt_size[2]*latt_size[3];
      GridCartesian     Grid(latt_size,simd_layout,mpi_layout);
      GridParallelRNG          pRNG(&Grid);      pRNG.SeedFixedIntegers(std::vector<int>({45,12,81,9}));

      LatticeColourMatrix z(&Grid);// random(pRNG,z);
      LatticeColourMatrix x(&Grid);// random(pRNG,x);
      LatticeColourMatrix y(&Grid);// random(pRNG,y);

      auto x_v = x.View();
      auto y_v = y.View();
      auto z_v = z.View();
      const uint64_t Nsite = x_v.size();
      const uint64_t nsimd = vComplex::Nsimd();
      const uint64_t    NN = Nsite*nsimd;
      for(int64_t i=0;i<Nwarm;i++){
	accelerator_loopN( sss, NN, {
	  uint64_t lane = sss % nsimd; 
	  uint64_t ss   = sss / nsimd;
	  auto xx = extractLane(lane,x_v[ss]);
	  auto yy = extractLane(lane,y_v[ss]);
	  auto zz = xx*yy;
	  insertLane(lane,x_v[ss],zz);
        });
      }
      double start=usecond();
      for(int64_t i=0;i<Nloop;i++){
	accelerator_loopN( sss, NN, {
	  uint64_t lane = sss % nsimd; 
	  uint64_t ss   = sss / nsimd;
	  auto xx = extractLane(lane,x_v[ss]);
	  auto yy = extractLane(lane,y_v[ss]);
	  auto zz = xx*yy;
	  insertLane(lane,x_v[ss],zz);
        });
      }
      double stop=usecond();
      double time = (stop-start)/Nloop*1000.0;
      
      double bytes=3.0*vol*Nc*Nc*sizeof(Complex);
      double footprint=2.0*vol*Nc*Nc*sizeof(Complex);
      double flops=Nc*Nc*(6.0+8.0+8.0)*vol;
      std::cout<<GridLogMessage<<std::setprecision(3) << lat<<"\t\t"<<footprint<<"    \t\t"<<bytes/time<<"\t\t" << flops/time<<std::endl;

    }



  std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
  std::cout<<GridLogMessage << "= Benchmarking SU3xSU3  ext/ins z = z + x*y"<<std::endl;
  std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
  std::cout<<GridLogMessage << "  L  "<<"\t\t"<<"bytes"<<"\t\t\t"<<"GB/s\t\t GFlop/s"<<std::endl;
  std::cout<<GridLogMessage << "----------------------------------------------------------"<<std::endl;

  for(int lat=LMIN;lat<=LMAX;lat+=LADD){

      Coordinate latt_size  ({lat*mpi_layout[0],lat*mpi_layout[1],lat*mpi_layout[2],lat*mpi_layout[3]});
      int64_t vol = latt_size[0]*latt_size[1]*latt_size[2]*latt_size[3];
      GridCartesian     Grid(latt_size,simd_layout,mpi_layout);
      GridParallelRNG          pRNG(&Grid);      pRNG.SeedFixedIntegers(std::vector<int>({45,12,81,9}));

      LatticeColourMatrix z(&Grid);// random(pRNG,z);
      LatticeColourMatrix x(&Grid);// random(pRNG,x);
      LatticeColourMatrix y(&Grid);// random(pRNG,y);

      auto x_v = x.View();
      auto y_v = y.View();
      auto z_v = z.View();
      const uint64_t Nsite = x_v.size();
      const uint64_t nsimd = vComplex::Nsimd();
      const uint64_t    NN = Nsite*nsimd;
      for(int64_t i=0;i<Nwarm;i++){
	accelerator_loopN( sss, NN, {
	  uint64_t lane = sss % nsimd; 
	  uint64_t ss   = sss / nsimd;
	  auto xx = extractLane(lane,x_v[ss]);
	  auto yy = extractLane(lane,y_v[ss]);
	  auto zz = extractLane(lane,z_v[ss]);
	  zz = zz + xx * yy;
	  insertLane(lane,z_v[ss],zz);
        });
      }
      double start=usecond();
      for(int64_t i=0;i<Nloop;i++){
	accelerator_loopN( sss, NN, {
	  uint64_t lane = sss % nsimd; 
	  uint64_t ss   = sss / nsimd;
	  auto xx = extractLane(lane,x_v[ss]);
	  auto yy = extractLane(lane,y_v[ss]);
	  auto zz = extractLane(lane,z_v[ss]);
	  zz = zz + xx * yy;
	  insertLane(lane,x_v[ss],zz);
        });
      }
      double stop=usecond();
      double time = (stop-start)/Nloop*1000.0;
      
      double bytes=4.0*vol*Nc*Nc*sizeof(Complex);
      double footprint=3.0*vol*Nc*Nc*sizeof(Complex);
      double flops=Nc*Nc*(6.0+8.0+8.0)*vol;
      std::cout<<GridLogMessage<<std::setprecision(3) << lat<<"\t\t"<<footprint<<"    \t\t"<<bytes/time<<"\t\t" << flops/time<<std::endl;

    }


  std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
  std::cout<<GridLogMessage << "= Benchmarking SU3xSU3  ext/ins w = z + x*y"<<std::endl;
  std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
  std::cout<<GridLogMessage << "  L  "<<"\t\t"<<"bytes"<<"\t\t\t"<<"GB/s\t\t GFlop/s"<<std::endl;
  std::cout<<GridLogMessage << "----------------------------------------------------------"<<std::endl;

  for(int lat=LMIN;lat<=LMAX;lat+=LADD){

      Coordinate latt_size  ({lat*mpi_layout[0],lat*mpi_layout[1],lat*mpi_layout[2],lat*mpi_layout[3]});
      int64_t vol = latt_size[0]*latt_size[1]*latt_size[2]*latt_size[3];
      GridCartesian     Grid(latt_size,simd_layout,mpi_layout);
      GridParallelRNG          pRNG(&Grid);      pRNG.SeedFixedIntegers(std::vector<int>({45,12,81,9}));

      LatticeColourMatrix z(&Grid);// random(pRNG,z);
      LatticeColourMatrix x(&Grid);// random(pRNG,x);
      LatticeColourMatrix y(&Grid);// random(pRNG,y);
      LatticeColourMatrix w(&Grid);// random(pRNG,y);

      auto x_v = x.View();
      auto y_v = y.View();
      auto z_v = z.View();
      auto w_v = z.View();
      const uint64_t Nsite = x_v.size();
      const uint64_t nsimd = vComplex::Nsimd();
      const uint64_t    NN = Nsite*nsimd;
      for(int64_t i=0;i<Nwarm;i++){
	accelerator_loopN( sss, NN, {
	  uint64_t lane = sss % nsimd; 
	  uint64_t ss   = sss / nsimd;
	  auto xx = extractLane(lane,x_v[ss]);
	  auto yy = extractLane(lane,y_v[ss]);
	  auto zz = extractLane(lane,z_v[ss]);
	  auto ww = zz + xx * yy;
	  insertLane(lane,w_v[ss],ww);
        });
      }
      double start=usecond();
      for(int64_t i=0;i<Nloop;i++){
	accelerator_loopN( sss, NN, {
	  uint64_t lane = sss % nsimd; 
	  uint64_t ss   = sss / nsimd;
	  auto xx = extractLane(lane,x_v[ss]);
	  auto yy = extractLane(lane,y_v[ss]);
	  auto zz = extractLane(lane,z_v[ss]);
	  auto ww = zz + xx * yy;
	  insertLane(lane,w_v[ss],ww);
        });
      }
      double stop=usecond();
      double time = (stop-start)/Nloop*1000.0;
      
      double bytes=4.0*vol*Nc*Nc*sizeof(Complex);
      double footprint=3.0*vol*Nc*Nc*sizeof(Complex);
      double flops=Nc*Nc*(6.0+8.0+8.0)*vol;
      std::cout<<GridLogMessage<<std::setprecision(3) << lat<<"\t\t"<<footprint<<"    \t\t"<<bytes/time<<"\t\t" << flops/time<<std::endl;

    }



#if 1
  std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
  std::cout<<GridLogMessage << "= Benchmarking SU3xSU3  x= x*y"<<std::endl;
  std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
  std::cout<<GridLogMessage << "  L  "<<"\t\t"<<"bytes"<<"\t\t\t"<<"GB/s\t\t GFlop/s"<<std::endl;
  std::cout<<GridLogMessage << "----------------------------------------------------------"<<std::endl;

  for(int lat=LMIN;lat<=LMAX;lat+=LADD){

      Coordinate latt_size  ({lat*mpi_layout[0],lat*mpi_layout[1],lat*mpi_layout[2],lat*mpi_layout[3]});
      int64_t vol = latt_size[0]*latt_size[1]*latt_size[2]*latt_size[3];
      GridCartesian     Grid(latt_size,simd_layout,mpi_layout);
      GridParallelRNG          pRNG(&Grid);      pRNG.SeedFixedIntegers(std::vector<int>({45,12,81,9}));

      LatticeColourMatrix z(&Grid); random(pRNG,z);
      LatticeColourMatrix x(&Grid); random(pRNG,x);
      LatticeColourMatrix y(&Grid); random(pRNG,y);

      for(int64_t i=0;i<Nwarm;i++){
	x=x*y;
      }
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

  for(int lat=LMIN;lat<=LMAX;lat+=LADD){

      Coordinate latt_size  ({lat*mpi_layout[0],lat*mpi_layout[1],lat*mpi_layout[2],lat*mpi_layout[3]});
      int64_t vol = latt_size[0]*latt_size[1]*latt_size[2]*latt_size[3];

      GridCartesian     Grid(latt_size,simd_layout,mpi_layout);
      GridParallelRNG          pRNG(&Grid);      pRNG.SeedFixedIntegers(std::vector<int>({45,12,81,9}));

      LatticeColourMatrix z(&Grid); random(pRNG,z);
      LatticeColourMatrix x(&Grid); random(pRNG,x);
      LatticeColourMatrix y(&Grid); random(pRNG,y);

      for(int64_t i=0;i<Nwarm;i++){
	z=x*y;
      }
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

  for(int lat=LMIN;lat<=LMAX;lat+=LADD){

      Coordinate latt_size  ({lat*mpi_layout[0],lat*mpi_layout[1],lat*mpi_layout[2],lat*mpi_layout[3]});
      int64_t vol = latt_size[0]*latt_size[1]*latt_size[2]*latt_size[3];

      GridCartesian     Grid(latt_size,simd_layout,mpi_layout);
      GridParallelRNG          pRNG(&Grid);      pRNG.SeedFixedIntegers(std::vector<int>({45,12,81,9}));

      LatticeColourMatrix z(&Grid); random(pRNG,z);
      LatticeColourMatrix x(&Grid); random(pRNG,x);
      LatticeColourMatrix y(&Grid); random(pRNG,y);

      for(int64_t i=0;i<Nwarm;i++){
	mult(z,x,y);
      }
      double start=usecond();
      for(int64_t i=0;i<Nloop;i++){
	mult(z,x,y);
	//	mac(z,x,y);
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

  for(int lat=LMIN;lat<=LMAX;lat+=LADD){

      Coordinate latt_size  ({lat*mpi_layout[0],lat*mpi_layout[1],lat*mpi_layout[2],lat*mpi_layout[3]});
      int64_t vol = latt_size[0]*latt_size[1]*latt_size[2]*latt_size[3];

      GridCartesian     Grid(latt_size,simd_layout,mpi_layout);
      GridParallelRNG          pRNG(&Grid);      pRNG.SeedFixedIntegers(std::vector<int>({45,12,81,9}));

      LatticeColourMatrix z(&Grid); random(pRNG,z);
      LatticeColourMatrix x(&Grid); random(pRNG,x);
      LatticeColourMatrix y(&Grid); random(pRNG,y);

      for(int64_t i=0;i<Nwarm;i++){
	mac(z,x,y);
      }
      double start=usecond();
      for(int64_t i=0;i<Nloop;i++){
	mac(z,x,y);
      }
      double stop=usecond();
      double time = (stop-start)/Nloop*1000.0;
      
      double bytes=4*vol*Nc*Nc*sizeof(Complex);
      double flops=Nc*Nc*(8+8+8)*vol;
      std::cout<<GridLogMessage<<std::setprecision(3) << lat<<"\t\t"<<bytes<<"   \t\t"<<bytes/time<<"\t\t" << flops/time<<std::endl;

    }
#endif
  Grid_finalize();
}
