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
            <<std::setw(11)<<"bytes"<<"MB/s uni (err/min/max)"<<"\t\t"<<"MB/s bidi (err/min/max)"<<std::endl;
};

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  Coordinate simd_layout = GridDefaultSimd(Nd,vComplexD::Nsimd());
  Coordinate mpi_layout  = GridDefaultMpi();
  int threads = GridThread::GetThreads();
  std::cout<<GridLogMessage << "Grid is setup to use "<<threads<<" threads"<<std::endl;

  int Nloop=100;
  int nmu=0;
  int maxlat=32;
  for(int mu=0;mu<Nd;mu++) if (mpi_layout[mu]>1) nmu++;

  std::cout << GridLogMessage << "Number of iterations to average: "<< Nloop << std::endl;
  std::vector<double> t_time(Nloop);
  time_statistics timestat;



  std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
  std::cout<<GridLogMessage << "= Benchmarking sequential halo exchange in "<<nmu<<" dimensions"<<std::endl;
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
	//	std::cout << " buffers " << std::hex << (uint64_t)&xbuf[mu][0] <<" " << (uint64_t)&rbuf[mu][0] <<std::endl;
      }
      uint64_t bytes=lat*lat*lat*Ls*sizeof(HalfSpinColourVectorD);

      int ncomm;

      for(int i=0;i<Nloop;i++){
      double start=usecond();
    
	ncomm=0;
	for(int mu=0;mu<4;mu++){
	
	  if (mpi_layout[mu]>1 ) {
	  
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
	}
	Grid.Barrier();
	double stop=usecond();
	t_time[i] = stop-start; // microseconds

      }

      timestat.statistics(t_time);
      
      double dbytes    = bytes*ppn;
      double xbytes    = dbytes*2.0*ncomm;
      double rbytes    = xbytes;
      double bidibytes = xbytes+rbytes;

    std::cout<<GridLogMessage << std::setw(4) << lat<<"\t"<<Ls<<"\t"
               <<std::setw(11) << bytes<< std::fixed << std::setprecision(1) << std::setw(7)
               <<std::right<< xbytes/timestat.mean<<"  "<< xbytes*timestat.err/(timestat.mean*timestat.mean)<< " "
               <<xbytes/timestat.max <<" "<< xbytes/timestat.min  
               << "\t\t"<<std::setw(7)<< bidibytes/timestat.mean<< "  " << bidibytes*timestat.err/(timestat.mean*timestat.mean) << " "
               << bidibytes/timestat.max << " " << bidibytes/timestat.min << std::endl;

      
    }
  }  


  std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
  std::cout<<GridLogMessage << "= Benchmarking concurrent STENCIL halo exchange in "<<nmu<<" dimensions"<<std::endl;
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
      Grid.ShmBufferFreeAll();
      uint64_t bytes = lat*lat*lat*Ls*sizeof(HalfSpinColourVectorD);
      for(int d=0;d<8;d++){
	xbuf[d] = (HalfSpinColourVectorD *)Grid.ShmBufferMalloc(bytes);
	rbuf[d] = (HalfSpinColourVectorD *)Grid.ShmBufferMalloc(bytes);
      }

      int ncomm;

      double dbytes;
      for(int i=0;i<Nloop;i++){
	double start=usecond();

	dbytes=0;
	ncomm=0;

	std::vector<CommsRequest_t> requests;

	for(int mu=0;mu<4;mu++){
	

	  if (mpi_layout[mu]>1 ) {
	  
	    ncomm++;
	    int comm_proc=1;
	    int xmit_to_rank;
	    int recv_from_rank;
	    Grid.ShiftedRanks(mu,comm_proc,xmit_to_rank,recv_from_rank);
	    dbytes+=
	      Grid.StencilSendToRecvFromBegin(requests,
					      (void *)&xbuf[mu][0],
					      xmit_to_rank,1,
					      (void *)&rbuf[mu][0],
					      recv_from_rank,1,
					      bytes,bytes,mu,sizeof(double));
	
	    comm_proc = mpi_layout[mu]-1;
	  
	    Grid.ShiftedRanks(mu,comm_proc,xmit_to_rank,recv_from_rank);
	    dbytes+=
	      Grid.StencilSendToRecvFromBegin(requests,
					      (void *)&xbuf[mu+4][0],
					      xmit_to_rank,1,
					      (void *)&rbuf[mu+4][0],
					      recv_from_rank,1,
					      bytes,bytes,mu+4,sizeof(double));
	  
	  }
	}
	Grid.StencilSendToRecvFromComplete(requests,0,sizeof(double));
	Grid.Barrier();
	double stop=usecond();
	t_time[i] = stop-start; // microseconds
	
      }

      timestat.statistics(t_time);

      dbytes=dbytes*ppn;
      double xbytes    = dbytes*0.5;
      //      double rbytes    = dbytes*0.5;
      double bidibytes = dbytes;

      std::cout<<GridLogMessage << std::setw(4) << lat<<"\t"<<Ls<<"\t"
               <<std::setw(11) << bytes<< std::fixed << std::setprecision(1) << std::setw(7)
               <<std::right<< xbytes/timestat.mean<<"  "<< xbytes*timestat.err/(timestat.mean*timestat.mean)<< " "
               <<xbytes/timestat.max <<" "<< xbytes/timestat.min  
               << "\t\t"<<std::setw(7)<< bidibytes/timestat.mean<< "  " << bidibytes*timestat.err/(timestat.mean*timestat.mean) << " "
               << bidibytes/timestat.max << " " << bidibytes/timestat.min << std::endl;


    }
  }    


  std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
  std::cout<<GridLogMessage << "= Benchmarking sequential STENCIL halo exchange in "<<nmu<<" dimensions"<<std::endl;
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
      Grid.ShmBufferFreeAll();
      uint64_t bytes=lat*lat*lat*Ls*sizeof(HalfSpinColourVectorD);
      for(int d=0;d<8;d++){
	xbuf[d] = (HalfSpinColourVectorD *)Grid.ShmBufferMalloc(bytes);
	rbuf[d] = (HalfSpinColourVectorD *)Grid.ShmBufferMalloc(bytes);
      }

      int ncomm;
      double dbytes;
      for(int i=0;i<Nloop;i++){
	double start=usecond();

	std::vector<CommsRequest_t> requests;
	dbytes=0;
	ncomm=0;
	for(int mu=0;mu<4;mu++){
	
	  if (mpi_layout[mu]>1 ) {
	  
	    ncomm++;
	    int comm_proc=1;
	    int xmit_to_rank;
	    int recv_from_rank;
	    
	    Grid.ShiftedRanks(mu,comm_proc,xmit_to_rank,recv_from_rank);
	    dbytes+=
	      Grid.StencilSendToRecvFromBegin(requests,
					      (void *)&xbuf[mu][0],
					      xmit_to_rank,1,
					      (void *)&rbuf[mu][0],
					      recv_from_rank,1,
					      bytes,bytes,mu,sizeof(double));
	    Grid.StencilSendToRecvFromComplete(requests,mu,sizeof(double));
	    requests.resize(0);

	    comm_proc = mpi_layout[mu]-1;
	  
	    Grid.ShiftedRanks(mu,comm_proc,xmit_to_rank,recv_from_rank);
	    dbytes+=
	      Grid.StencilSendToRecvFromBegin(requests,
					      (void *)&xbuf[mu+4][0],
					      xmit_to_rank,1,
					      (void *)&rbuf[mu+4][0],
					      recv_from_rank,1,
					      bytes,bytes,mu+4,sizeof(double));
	    Grid.StencilSendToRecvFromComplete(requests,mu+4,sizeof(double));
	    requests.resize(0);
	  
	  }
	}
	Grid.Barrier();
	double stop=usecond();
	t_time[i] = stop-start; // microseconds
	
      }

      timestat.statistics(t_time);

      dbytes=dbytes*ppn;
      double xbytes    = dbytes*0.5;
      //      double rbytes    = dbytes*0.5;
      double bidibytes = dbytes;


      std::cout<<GridLogMessage << std::setw(4) << lat<<"\t"<<Ls<<"\t"
               <<std::setw(11) << bytes<< std::fixed << std::setprecision(1) << std::setw(7)
               <<std::right<< xbytes/timestat.mean<<"  "<< xbytes*timestat.err/(timestat.mean*timestat.mean)<< " "
               <<xbytes/timestat.max <<" "<< xbytes/timestat.min  
               << "\t\t"<<std::setw(7)<< bidibytes/timestat.mean<< "  " << bidibytes*timestat.err/(timestat.mean*timestat.mean) << " "
               << bidibytes/timestat.max << " " << bidibytes/timestat.min << std::endl;
 
    }
  }    


#ifdef GRID_OMP
  std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
  std::cout<<GridLogMessage << "= Benchmarking threaded STENCIL halo exchange in "<<nmu<<" dimensions"<<std::endl;
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
      Grid.ShmBufferFreeAll();
      uint64_t bytes = lat*lat*lat*Ls*sizeof(HalfSpinColourVectorD);
      for(int d=0;d<8;d++){
	xbuf[d] = (HalfSpinColourVectorD *)Grid.ShmBufferMalloc(bytes);
	rbuf[d] = (HalfSpinColourVectorD *)Grid.ShmBufferMalloc(bytes);
      }

      int ncomm;
      double dbytes;
      for(int i=0;i<Nloop;i++){
	double start=usecond();

	std::vector<CommsRequest_t> requests;
	dbytes=0;
	ncomm=0;

	thread_for(dir,8,{

	  double tbytes;
	  int mu =dir % 4;

	  if (mpi_layout[mu]>1 ) {
	  
	    ncomm++;
	    int xmit_to_rank;
	    int recv_from_rank;
	    if ( dir == mu ) { 
	      int comm_proc=1;
	      Grid.ShiftedRanks(mu,comm_proc,xmit_to_rank,recv_from_rank);
	    } else { 
	      int comm_proc = mpi_layout[mu]-1;
	      Grid.ShiftedRanks(mu,comm_proc,xmit_to_rank,recv_from_rank);
	    }
            int tid = omp_get_thread_num();
	    tbytes= Grid.StencilSendToRecvFrom((void *)&xbuf[dir][0], xmit_to_rank,1,
					       (void *)&rbuf[dir][0], recv_from_rank,1, bytes,tid,sizeof(double));

	    thread_critical { dbytes+=tbytes; }
	  }
        });
	Grid.Barrier();
	double stop=usecond();
	t_time[i] = stop-start; // microseconds
      }

      timestat.statistics(t_time);

      dbytes=dbytes*ppn;
      double xbytes    = dbytes*0.5;
      //      double rbytes    = dbytes*0.5;
      double bidibytes = dbytes;


      std::cout<<GridLogMessage << std::setw(4) << lat<<"\t"<<Ls<<"\t"
               <<std::setw(11) << bytes<< std::fixed << std::setprecision(1) << std::setw(7)
               <<std::right<< xbytes/timestat.mean<<"  "<< xbytes*timestat.err/(timestat.mean*timestat.mean)<< " "
               <<xbytes/timestat.max <<" "<< xbytes/timestat.min  
               << "\t\t"<<std::setw(7)<< bidibytes/timestat.mean<< "  " << bidibytes*timestat.err/(timestat.mean*timestat.mean) << " "
               << bidibytes/timestat.max << " " << bidibytes/timestat.min << std::endl;
 
    }
  }    
#endif
  std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
  std::cout<<GridLogMessage << "= All done; Bye Bye"<<std::endl;
  std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;

  Grid_finalize();
}
