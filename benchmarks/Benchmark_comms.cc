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
using namespace Grid::QCD;

void statistics(std::vector<double> v, double &mean, double &std_err){
      double sum = std::accumulate(v.begin(), v.end(), 0.0);
      mean = sum / v.size();

      std::vector<double> diff(v.size());
      std::transform(v.begin(), v.end(), diff.begin(), [mean](double x) { return x - mean; });
      double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
      std_err = std::sqrt(sq_sum / (v.size()*(v.size() - 1)));
}


int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  std::vector<int> simd_layout = GridDefaultSimd(Nd,vComplexD::Nsimd());
  std::vector<int> mpi_layout  = GridDefaultMpi();
  int threads = GridThread::GetThreads();
  std::cout<<GridLogMessage << "Grid is setup to use "<<threads<<" threads"<<std::endl;

  int Nloop=500;
  int nmu=0;
  int maxlat=24;
  for(int mu=0;mu<Nd;mu++) if (mpi_layout[mu]>1) nmu++;

  std::cout << GridLogMessage << "Number of iterations to average: "<< Nloop << std::endl;
  std::vector<double> t_time(Nloop);
  double mean_time, std_err_time;

  std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
  std::cout<<GridLogMessage << "= Benchmarking concurrent halo exchange in "<<nmu<<" dimensions"<<std::endl;
  std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
  std::cout<<GridLogMessage << "  L  "<<"\t\t"<<" Ls  "<<"\t\t"<<"bytes"<<"\t\t"<<"MB/s uni"<<"\t\t"<<"MB/s bidi"<<std::endl;
  for(int lat=4;lat<=maxlat;lat+=4){
    for(int Ls=8;Ls<=32;Ls*=2){

      std::vector<int> latt_size  ({lat*mpi_layout[0],
      				    lat*mpi_layout[1],
      				    lat*mpi_layout[2],
      				    lat*mpi_layout[3]});

      GridCartesian     Grid(latt_size,simd_layout,mpi_layout);

      std::vector<std::vector<HalfSpinColourVectorD> > xbuf(8,std::vector<HalfSpinColourVectorD>(lat*lat*lat*Ls));
      std::vector<std::vector<HalfSpinColourVectorD> > rbuf(8,std::vector<HalfSpinColourVectorD>(lat*lat*lat*Ls));

      int ncomm;
      int bytes=lat*lat*lat*Ls*sizeof(HalfSpinColourVectorD);

      for(int i=0;i<Nloop;i++){
      double start=usecond();

	std::vector<CartesianCommunicator::CommsRequest_t> requests;

	ncomm=0;
	for(int mu=0;mu<4;mu++){
	
	  if (mpi_layout[mu]>1 ) {
	  
	    ncomm++;
	    int comm_proc=1;
	    int xmit_to_rank;
	    int recv_from_rank;
	    
	    Grid.ShiftedRanks(mu,comm_proc,xmit_to_rank,recv_from_rank);
	    Grid.SendToRecvFromBegin(requests,
				   (void *)&xbuf[mu][0],
				   xmit_to_rank,
				   (void *)&rbuf[mu][0],
				   recv_from_rank,
				   bytes);
	
	    comm_proc = mpi_layout[mu]-1;
	  
	    Grid.ShiftedRanks(mu,comm_proc,xmit_to_rank,recv_from_rank);
	    Grid.SendToRecvFromBegin(requests,
				     (void *)&xbuf[mu+4][0],
				     xmit_to_rank,
				     (void *)&rbuf[mu+4][0],
				     recv_from_rank,
				     bytes);
	  
	  }
	}
	Grid.SendToRecvFromComplete(requests);
	Grid.Barrier();
  double stop=usecond();
  t_time[i] = stop-start; // microseconds
      }

      statistics(t_time, mean_time, std_err_time);

      double dbytes    = bytes;
      double xbytes    = dbytes*2.0*ncomm;
      double rbytes    = xbytes;
      double bidibytes = xbytes+rbytes;

      std::cout<<GridLogMessage << lat<<"\t\t"<<Ls<<"\t\t"<< std::setw(15) <<bytes<<"\t"
               <<xbytes/mean_time<<" +- "<< xbytes*std_err_time/(mean_time*mean_time)
               << "\t\t"<<bidibytes/mean_time<< " +- " << bidibytes*std_err_time/(mean_time*mean_time) << std::endl;

    }
  }    


  std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
  std::cout<<GridLogMessage << "= Benchmarking sequential halo exchange in "<<nmu<<" dimensions"<<std::endl;
  std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
  std::cout<<GridLogMessage << "  L  "<<"\t\t"<<" Ls  "<<"\t\t"<<"bytes"<<"\t\t"<<"MB/s uni"<<"\t\t"<<"MB/s bidi"<<std::endl;


  for(int lat=4;lat<=maxlat;lat+=4){
    for(int Ls=8;Ls<=32;Ls*=2){

      std::vector<int> latt_size  ({lat,lat,lat,lat});

      GridCartesian     Grid(latt_size,simd_layout,mpi_layout);

      std::vector<std::vector<HalfSpinColourVectorD> > xbuf(8,std::vector<HalfSpinColourVectorD>(lat*lat*lat*Ls));
      std::vector<std::vector<HalfSpinColourVectorD> > rbuf(8,std::vector<HalfSpinColourVectorD>(lat*lat*lat*Ls));


      int ncomm;
      int bytes=lat*lat*lat*Ls*sizeof(HalfSpinColourVectorD);

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
	      std::vector<CartesianCommunicator::CommsRequest_t> requests;
	      Grid.ShiftedRanks(mu,comm_proc,xmit_to_rank,recv_from_rank);
	      Grid.SendToRecvFromBegin(requests,
				       (void *)&xbuf[mu][0],
				       xmit_to_rank,
				       (void *)&rbuf[mu][0],
				       recv_from_rank,
				       bytes);
	      Grid.SendToRecvFromComplete(requests);
	    }

	    comm_proc = mpi_layout[mu]-1;
	    {
	      std::vector<CartesianCommunicator::CommsRequest_t> requests;
	      Grid.ShiftedRanks(mu,comm_proc,xmit_to_rank,recv_from_rank);
	      Grid.SendToRecvFromBegin(requests,
				       (void *)&xbuf[mu+4][0],
				       xmit_to_rank,
				       (void *)&rbuf[mu+4][0],
				       recv_from_rank,
				       bytes);
	      Grid.SendToRecvFromComplete(requests);
	    }
	  }
	}
	Grid.Barrier();
      double stop=usecond();
    t_time[i] = stop-start; // microseconds

      }

      statistics(t_time, mean_time, std_err_time);
      
      double dbytes    = bytes;
      double xbytes    = Nloop*dbytes*2.0*ncomm;
      double rbytes    = xbytes;
      double bidibytes = xbytes+rbytes;

     std::cout<<GridLogMessage << lat<<"\t\t"<<Ls<<"\t\t"<< std::setw(15) <<bytes<<"\t"
               <<xbytes/mean_time<<" +- "<< xbytes*std_err_time/(mean_time*mean_time)
               << "\t\t"<<bidibytes/mean_time<< " +- " << bidibytes*std_err_time/(mean_time*mean_time) << std::endl;
    }
  }  


  std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
  std::cout<<GridLogMessage << "= Benchmarking concurrent STENCIL halo exchange in "<<nmu<<" dimensions"<<std::endl;
  std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
  std::cout<<GridLogMessage << "  L  "<<"\t\t"<<" Ls  "<<"\t\t"<<"bytes"<<"\t\t"<<"MB/s uni"<<"\t\t"<<"MB/s bidi"<<std::endl;

  for(int lat=4;lat<=maxlat;lat+=4){
    for(int Ls=8;Ls<=32;Ls*=2){

      std::vector<int> latt_size  ({lat*mpi_layout[0],
      				    lat*mpi_layout[1],
      				    lat*mpi_layout[2],
      				    lat*mpi_layout[3]});

      GridCartesian     Grid(latt_size,simd_layout,mpi_layout);

      std::vector<HalfSpinColourVectorD *> xbuf(8);
      std::vector<HalfSpinColourVectorD *> rbuf(8);
      Grid.ShmBufferFreeAll();
      for(int d=0;d<8;d++){
	xbuf[d] = (HalfSpinColourVectorD *)Grid.ShmBufferMalloc(lat*lat*lat*Ls*sizeof(HalfSpinColourVectorD));
	rbuf[d] = (HalfSpinColourVectorD *)Grid.ShmBufferMalloc(lat*lat*lat*Ls*sizeof(HalfSpinColourVectorD));
      }

      int ncomm;
      int bytes=lat*lat*lat*Ls*sizeof(HalfSpinColourVectorD);

      for(int i=0;i<Nloop;i++){
      double start=usecond();

	std::vector<CartesianCommunicator::CommsRequest_t> requests;

	ncomm=0;
	for(int mu=0;mu<4;mu++){
	
	  if (mpi_layout[mu]>1 ) {
	  
	    ncomm++;
	    int comm_proc=1;
	    int xmit_to_rank;
	    int recv_from_rank;
	    
	    Grid.ShiftedRanks(mu,comm_proc,xmit_to_rank,recv_from_rank);
	    Grid.StencilSendToRecvFromBegin(requests,
					    (void *)&xbuf[mu][0],
					    xmit_to_rank,
					    (void *)&rbuf[mu][0],
					    recv_from_rank,
					    bytes);
	
	    comm_proc = mpi_layout[mu]-1;
	  
	    Grid.ShiftedRanks(mu,comm_proc,xmit_to_rank,recv_from_rank);
	    Grid.StencilSendToRecvFromBegin(requests,
					    (void *)&xbuf[mu+4][0],
					    xmit_to_rank,
					    (void *)&rbuf[mu+4][0],
					    recv_from_rank,
					    bytes);
	  
	  }
	}
	Grid.StencilSendToRecvFromComplete(requests);
	Grid.Barrier();
      double stop=usecond();
    t_time[i] = stop-start; // microseconds

      }

      double dbytes    = bytes;
      double xbytes    = Nloop*dbytes*2.0*ncomm;
      double rbytes    = xbytes;
      double bidibytes = xbytes+rbytes;

   std::cout<<GridLogMessage << lat<<"\t\t"<<Ls<<"\t\t"<< std::setw(15) <<bytes<<"\t"
               <<xbytes/mean_time<<" +- "<< xbytes*std_err_time/(mean_time*mean_time)
               << "\t\t"<<bidibytes/mean_time<< " +- " << bidibytes*std_err_time/(mean_time*mean_time) << std::endl;
    }
  }    


  std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
  std::cout<<GridLogMessage << "= Benchmarking sequential STENCIL halo exchange in "<<nmu<<" dimensions"<<std::endl;
  std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
  std::cout<<GridLogMessage << "  L  "<<"\t\t"<<" Ls  "<<"\t\t"<<"bytes"<<"\t\t"<<"MB/s uni"<<"\t\t"<<"MB/s bidi"<<std::endl;

  for(int lat=4;lat<=maxlat;lat+=4){
    for(int Ls=8;Ls<=32;Ls*=2){

      std::vector<int> latt_size  ({lat*mpi_layout[0],
      				    lat*mpi_layout[1],
      				    lat*mpi_layout[2],
      				    lat*mpi_layout[3]});

      GridCartesian     Grid(latt_size,simd_layout,mpi_layout);

      std::vector<HalfSpinColourVectorD *> xbuf(8);
      std::vector<HalfSpinColourVectorD *> rbuf(8);
      Grid.ShmBufferFreeAll();
      for(int d=0;d<8;d++){
	xbuf[d] = (HalfSpinColourVectorD *)Grid.ShmBufferMalloc(lat*lat*lat*Ls*sizeof(HalfSpinColourVectorD));
	rbuf[d] = (HalfSpinColourVectorD *)Grid.ShmBufferMalloc(lat*lat*lat*Ls*sizeof(HalfSpinColourVectorD));
      }

      int ncomm;
      int bytes=lat*lat*lat*Ls*sizeof(HalfSpinColourVectorD);

      for(int i=0;i<Nloop;i++){
      double start=usecond();

	std::vector<CartesianCommunicator::CommsRequest_t> requests;

	ncomm=0;
	for(int mu=0;mu<4;mu++){
	
	  if (mpi_layout[mu]>1 ) {
	  
	    ncomm++;
	    int comm_proc=1;
	    int xmit_to_rank;
	    int recv_from_rank;
	    
	    Grid.ShiftedRanks(mu,comm_proc,xmit_to_rank,recv_from_rank);
	    Grid.StencilSendToRecvFromBegin(requests,
					    (void *)&xbuf[mu][0],
					    xmit_to_rank,
					    (void *)&rbuf[mu][0],
					    recv_from_rank,
					    bytes);
	    Grid.StencilSendToRecvFromComplete(requests);
	    requests.resize(0);

	    comm_proc = mpi_layout[mu]-1;
	  
	    Grid.ShiftedRanks(mu,comm_proc,xmit_to_rank,recv_from_rank);
	    Grid.StencilSendToRecvFromBegin(requests,
					    (void *)&xbuf[mu+4][0],
					    xmit_to_rank,
					    (void *)&rbuf[mu+4][0],
					    recv_from_rank,
					    bytes);
	    Grid.StencilSendToRecvFromComplete(requests);
	    requests.resize(0);
	  
	  }
	}
	    Grid.Barrier();
      double stop=usecond();
      t_time[i] = stop-start; // microseconds

      }

      double dbytes    = bytes;
      double xbytes    = Nloop*dbytes*2.0*ncomm;
      double rbytes    = xbytes;
      double bidibytes = xbytes+rbytes;


      std::cout<<GridLogMessage << lat<<"\t\t"<<Ls<<"\t\t"<< std::setw(15) <<bytes<<"\t"
               <<xbytes/mean_time<<" +- "<< xbytes*std_err_time/(mean_time*mean_time)
               << "\t\t"<<bidibytes/mean_time<< " +- " << bidibytes*std_err_time/(mean_time*mean_time) << std::endl;
 
    }
  }    

  Grid_finalize();
}
