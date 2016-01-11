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
#include <Grid.h>

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  std::vector<int> simd_layout = GridDefaultSimd(Nd,vComplexD::Nsimd());
  std::vector<int> mpi_layout  = GridDefaultMpi();
  int threads = GridThread::GetThreads();
  std::cout<<GridLogMessage << "Grid is setup to use "<<threads<<" threads"<<std::endl;

  int Nloop=10;
  int nmu=0;
  for(int mu=0;mu<4;mu++) if (mpi_layout[mu]>1) nmu++;
  std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
  std::cout<<GridLogMessage << "= Benchmarking concurrent halo exchange in "<<nmu<<" dimensions"<<std::endl;
  std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
  std::cout<<GridLogMessage << "  L  "<<"\t\t"<<" Ls  "<<"\t\t"<<"bytes"<<"\t\t"<<"MB/s uni"<<"\t\t"<<"MB/s bidi"<<std::endl;



  for(int lat=4;lat<=32;lat+=2){
    for(int Ls=1;Ls<=16;Ls*=2){

      std::vector<int> latt_size  ({lat*mpi_layout[0],
      				    lat*mpi_layout[1],
      				    lat*mpi_layout[2],
      				    lat*mpi_layout[3]});

      GridCartesian     Grid(latt_size,simd_layout,mpi_layout);

      std::vector<std::vector<HalfSpinColourVectorD> > xbuf(8,std::vector<HalfSpinColourVectorD>(lat*lat*lat*Ls));
      std::vector<std::vector<HalfSpinColourVectorD> > rbuf(8,std::vector<HalfSpinColourVectorD>(lat*lat*lat*Ls));

      int ncomm;
      int bytes=lat*lat*lat*Ls*sizeof(HalfSpinColourVectorD);

      double start=usecond();
      for(int i=0;i<Nloop;i++){

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

      }
      double stop=usecond();

      double dbytes    = bytes;
      double xbytes    = Nloop*dbytes*2.0*ncomm;
      double rbytes    = xbytes;
      double bidibytes = xbytes+rbytes;

      double time = stop-start; // microseconds

      std::cout<<GridLogMessage << lat<<"\t\t"<<Ls<<"\t\t"<<bytes<<"\t\t"<<xbytes/time<<"\t\t"<<bidibytes/time<<std::endl;
    }
  }    


  std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
  std::cout<<GridLogMessage << "= Benchmarking sequential halo exchange in "<<nmu<<" dimensions"<<std::endl;
  std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
  std::cout<<GridLogMessage << "  L  "<<"\t\t"<<" Ls  "<<"\t\t"<<"bytes"<<"\t\t"<<"MB/s uni"<<"\t\t"<<"MB/s bidi"<<std::endl;


  for(int lat=4;lat<=32;lat+=2){
    for(int Ls=1;Ls<=16;Ls*=2){

      std::vector<int> latt_size  ({lat,lat,lat,lat});

      GridCartesian     Grid(latt_size,simd_layout,mpi_layout);

      std::vector<std::vector<HalfSpinColourVectorD> > xbuf(8,std::vector<HalfSpinColourVectorD>(lat*lat*lat*Ls));
      std::vector<std::vector<HalfSpinColourVectorD> > rbuf(8,std::vector<HalfSpinColourVectorD>(lat*lat*lat*Ls));


      int ncomm;
      int bytes=lat*lat*lat*Ls*sizeof(HalfSpinColourVectorD);

      double start=usecond();
      for(int i=0;i<Nloop;i++){
    
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
      }

      double stop=usecond();
      
      double dbytes    = bytes;
      double xbytes    = Nloop*dbytes*2.0*ncomm;
      double rbytes    = xbytes;
      double bidibytes = xbytes+rbytes;

      double time = stop-start;

      std::cout<<GridLogMessage << lat<<"\t\t"<<Ls<<"\t\t"<<bytes<<"\t\t"<<xbytes/time<<"\t\t"<<bidibytes/time<<std::endl;
    }
  }  



  Grid_finalize();
}
