    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./benchmarks/Benchmark_comms.cc

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>

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

struct time_statistics{
  double mean;
  double err;
  double min;
  double max;

  void statistics(std::vector<double> v){
      double sum = std::accumulate(v.begin(), v.end(), 0.0);
      mean = sum / v.size();

      std::vector<double> diff(v.size());
      std::transform(v.begin(), v.end(), diff.begin(), [=](double x) { return x - mean; });
      double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
      err = std::sqrt(sq_sum / (v.size()*(v.size() - 1)));

      auto result = std::minmax_element(v.begin(), v.end());
      min = *result.first;
      max = *result.second;
}
};

void header(){
  std::cout <<GridLogMessage << " L  "<<"\t"<<" Ls  "<<"\t"
            <<std::setw(11)<<"bytes\t\t"<<"MB/s uni"<<"\t"<<"MB/s bidi"<<std::endl;
};

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  Coordinate simd_layout = GridDefaultSimd(Nd,vComplexD::Nsimd());
  Coordinate mpi_layout  = GridDefaultMpi();
  int threads = GridThread::GetThreads();
  std::cout<<GridLogMessage << "Grid is setup to use "<<threads<<" threads"<<std::endl;

  int Nloop=250;
  int nmu=0;
  int maxlat=32;
  for(int mu=0;mu<Nd;mu++) if (mpi_layout[mu]>1) nmu++;

  std::cout << GridLogMessage << "Number of iterations to average: "<< Nloop << std::endl;
  std::vector<double> t_time(Nloop);
  //  time_statistics timestat;

  std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
  std::cout<<GridLogMessage << "= Benchmarking sequential halo exchange from host memory "<<std::endl;
  std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
  header();

  for(int lat=8;lat<=maxlat;lat+=4){
    for(int Ls=8;Ls<=8;Ls*=2){

      Coordinate latt_size  ({lat*mpi_layout[0],
	                      lat*mpi_layout[1],
      			      lat*mpi_layout[2],
      			      lat*mpi_layout[3]});

      GridCartesian     Grid(latt_size,simd_layout,mpi_layout);
      RealD Nrank = Grid._Nprocessors;
      RealD Nnode = Grid.NodeCount();
      RealD ppn = Nrank/Nnode;

      std::vector<std::vector<HalfSpinColourVectorD> > xbuf(8);
      std::vector<std::vector<HalfSpinColourVectorD> > rbuf(8);

      for(int mu=0;mu<8;mu++){
	xbuf[mu].resize(lat*lat*lat*Ls);
	rbuf[mu].resize(lat*lat*lat*Ls);
      }
      uint64_t bytes=lat*lat*lat*Ls*sizeof(HalfSpinColourVectorD);

      int ncomm;

      for(int mu=0;mu<4;mu++){
	if (mpi_layout[mu]>1 ) {
	double start=usecond();
	for(int i=0;i<Nloop;i++){

	  ncomm=0;
	
	  
	    ncomm++;
	    int comm_proc=1;
	    int xmit_to_rank;
	    int recv_from_rank;
	    
	    {
	      std::vector<CommsRequest_t> requests;
	      Grid.ShiftedRanks(mu,comm_proc,xmit_to_rank,recv_from_rank);
	      Grid.SendToRecvFrom((void *)&xbuf[mu][0],
				  xmit_to_rank,
				  (void *)&rbuf[mu][0],
				  recv_from_rank,
				  bytes);
	    }

	    comm_proc = mpi_layout[mu]-1;
	    {
	      std::vector<CommsRequest_t> requests;
	      Grid.ShiftedRanks(mu,comm_proc,xmit_to_rank,recv_from_rank);
	      Grid.SendToRecvFrom((void *)&xbuf[mu+4][0],
				  xmit_to_rank,
				  (void *)&rbuf[mu+4][0],
				  recv_from_rank,
				  bytes);
	    }
	}
	Grid.Barrier();
	double stop=usecond();
        double mean=(stop-start)/Nloop;      
      double dbytes    = bytes*ppn;
      double xbytes    = dbytes*2.0*ncomm;
      double rbytes    = xbytes;
      double bidibytes = xbytes+rbytes;

      std::cout<<GridLogMessage << std::setw(4) << lat<<"\t"<<Ls<<"\t"
               <<std::setw(11) << bytes<< std::fixed << std::setprecision(1) << std::setw(7)<<" "
               <<std::right<< xbytes/mean<<"  "
               << "\t\t"<<std::setw(7)<< bidibytes/mean<< std::endl;


	
	}
      }


      
    }
  }

  std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
  std::cout<<GridLogMessage << "= Benchmarking sequential halo exchange from GPU memory "<<std::endl;
  std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
  header();

  for(int lat=8;lat<=maxlat;lat+=4){
    for(int Ls=8;Ls<=8;Ls*=2){

      Coordinate latt_size  ({lat*mpi_layout[0],
	                      lat*mpi_layout[1],
      			      lat*mpi_layout[2],
      			      lat*mpi_layout[3]});

      GridCartesian     Grid(latt_size,simd_layout,mpi_layout);
      RealD Nrank = Grid._Nprocessors;
      RealD Nnode = Grid.NodeCount();
      RealD ppn = Nrank/Nnode;


      std::vector<HalfSpinColourVectorD *> xbuf(8);
      std::vector<HalfSpinColourVectorD *> rbuf(8);

      uint64_t bytes = lat*lat*lat*Ls*sizeof(HalfSpinColourVectorD);
      for(int d=0;d<8;d++){
	xbuf[d] = (HalfSpinColourVectorD *)acceleratorAllocDevice(bytes);
	rbuf[d] = (HalfSpinColourVectorD *)acceleratorAllocDevice(bytes);
      }

      int ncomm;

      for(int mu=0;mu<4;mu++){
	if (mpi_layout[mu]>1 ) {
	double start=usecond();
	for(int i=0;i<Nloop;i++){

	  ncomm=0;
	
	  
	    ncomm++;
	    int comm_proc=1;
	    int xmit_to_rank;
	    int recv_from_rank;
	    
	    {
	      std::vector<CommsRequest_t> requests;
	      Grid.ShiftedRanks(mu,comm_proc,xmit_to_rank,recv_from_rank);
	      Grid.SendToRecvFrom((void *)&xbuf[mu][0],
				  xmit_to_rank,
				  (void *)&rbuf[mu][0],
				  recv_from_rank,
				  bytes);
	    }

	    comm_proc = mpi_layout[mu]-1;
	    {
	      std::vector<CommsRequest_t> requests;
	      Grid.ShiftedRanks(mu,comm_proc,xmit_to_rank,recv_from_rank);
	      Grid.SendToRecvFrom((void *)&xbuf[mu+4][0],
				  xmit_to_rank,
				  (void *)&rbuf[mu+4][0],
				  recv_from_rank,
				  bytes);
	    }
	}
	Grid.Barrier();
	double stop=usecond();
        double mean=(stop-start)/Nloop;      
      double dbytes    = bytes*ppn;
      double xbytes    = dbytes*2.0*ncomm;
      double rbytes    = xbytes;
      double bidibytes = xbytes+rbytes;

      std::cout<<GridLogMessage << std::setw(4) << lat<<"\t"<<Ls<<"\t"
               <<std::setw(11) << bytes<< std::fixed << std::setprecision(1) << std::setw(7)<<" "
               <<std::right<< xbytes/mean<<"  "
               << "\t\t"<<std::setw(7)<< bidibytes/mean<< std::endl;


	
	}
      }

      for(int d=0;d<8;d++){
	acceleratorFreeDevice(xbuf[d]);
	acceleratorFreeDevice(rbuf[d]);
      }

      
    }
  }


  std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
  std::cout<<GridLogMessage << "= All done; Bye Bye"<<std::endl;
  std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;

  Grid_finalize();
}
