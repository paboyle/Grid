    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./benchmarks/Benchmark_memory_asynch.cc

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


  std::vector<int> simd_layout = GridDefaultSimd(Nd,vReal::Nsimd());
  std::vector<int> mpi_layout  = GridDefaultMpi();

  int threads = GridThread::GetThreads();
  std::cout<<GridLogMessage << "Grid is setup to use "<<threads<<" threads"<<std::endl;
  
  std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
  std::cout<<GridLogMessage << "= Benchmarking READ bandwidth"<<std::endl;
  std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
  std::cout<<GridLogMessage << "  L  "<<"\t\t"<<"bytes"<<"\t\t"<<"bytes/thread"<<"\t\t\t"<<"GB/s"<<"\t\t\t"<<"GB/s per thread"<<std::endl;
  std::cout<<GridLogMessage << "----------------------------------------------------------"<<std::endl;

  const int lmax = 16536*16;
  for(int lat=4;lat<=lmax;lat*=2){

    int Nloop=lmax*4/lat;

    std::vector<int> latt_size  ({2*mpi_layout[0],2*mpi_layout[1],4*mpi_layout[2],lat*mpi_layout[3]});

    GridCartesian     Grid(latt_size,simd_layout,mpi_layout);

    int vol = latt_size[0]*latt_size[1]*latt_size[2]*latt_size[3]*threads;

    Vec tsum; tsum = zero;

    GridParallelRNG          pRNG(&Grid);      pRNG.SeedRandomDevice();

    std::vector<double> stop(threads);
    Vector<Vec> sum(threads);

    std::vector<LatticeVec> x(threads,&Grid);
    for(int t=0;t<threads;t++){
      //      random(pRNG,x[t]);
    }

    double start=usecond();
PARALLEL_FOR_LOOP
    for(int t=0;t<threads;t++){

      sum[t] = x[t]._odata[0];
      for(int i=0;i<Nloop;i++){
	for(auto ss=x[t].begin();ss<x[t].end();ss++){
	  sum[t] = sum[t]+x[t]._odata[ss];
	}
      }
      stop[t]=usecond();
    }

    double max_stop=stop[0];
    double min_stop=stop[0];
    
    for(int t=0;t<threads;t++){
      tsum+=sum[t];
      if ( stop[t]<min_stop ) min_stop=stop[t];
      if ( stop[t]>max_stop ) max_stop=stop[t];
    }

    

    double max_time = (max_stop-start)/Nloop*1000;
    double min_time = (min_stop-start)/Nloop*1000;
      
    double bytes=vol*Nvec*sizeof(Real);
    std::cout<<GridLogMessage<<std::setprecision(3) << lat<<"\t\t"<<bytes<<"\t\t"<<bytes/threads<<"\t\t"<<bytes/max_time<<" - "<< bytes/min_time<<"\t\t"<<bytes/min_time/threads <<std::endl;
    
  }    

  Grid_finalize();
}
