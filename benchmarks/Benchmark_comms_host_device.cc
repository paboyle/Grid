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


  /////////////////////////////////////////////////////////////////////////////
  // Halo exchange through the PRODUCTION stencil path.
  //
  //   Prepare -> PollDtoH -> CopySynchronise -> Begin -> PollIRecv -> Complete
  //
  // verbatim as CartesianStencil::CommunicateBegin/Complete drive it, so this
  // measures the code that actually runs inside the Dslash rather than a
  // simplified stand-in built on SendToRecvFrom.
  //
  // THREE byte counts are reported, because they differ and conflating them
  // has already produced an impossible number (272 GB/s/node against a
  // 200 GB/s wire):
  //
  //   offered   : the whole halo this rank presents, all 8 directions.
  //   off-node  : the subset whose partner is on another node, classified
  //               GEOMETRICALLY via IsOffNode() -- independent of --shm-mpi.
  //   to-MPI    : what StencilSendToRecvFrom{Prepare,Begin} actually handed to
  //               MPI, i.e. their return value, exactly as the Stencil counts
  //               it.  Under --shm-mpi 1 Grid delegates intranode transfers to
  //               MPI as well, so this approaches "offered"; under --shm-mpi 0
  //               Grid moves them itself and this approaches "off-node".
  //
  // Only the off-node column may be compared against wire speed.
  /////////////////////////////////////////////////////////////////////////////
  std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
  std::cout<<GridLogMessage << "= Benchmarking STENCIL-PATH halo exchange from GPU memory "<<std::endl;
  std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
  std::cout<<GridLogMessage << "  L   Ls   offered   off-node    to-MPI      time   off-node   to-MPI"<<std::endl;
  std::cout<<GridLogMessage << "            MB/nd      MB/nd      MB/nd    us/call   GB/s/nd  GB/s/nd"<<std::endl;

  for(int lat=8;lat<=maxlat;lat+=4){
    for(int Ls=8;Ls<=8;Ls*=2){

      Coordinate latt_size  ({lat*mpi_layout[0],
	                      lat*mpi_layout[1],
	                      lat*mpi_layout[2],
	                      lat*mpi_layout[3]});

      GridCartesian Grid(latt_size,simd_layout,mpi_layout);
      RealD Nrank = Grid._Nprocessors;
      RealD Nnode = Grid.NodeCount();

      uint64_t bytes = (uint64_t)lat*lat*lat*Ls*sizeof(HalfSpinColourVectorD);

      // Buffers MUST come from the shared-memory pool, not acceleratorAllocDevice.
      // The Stencil allocates u_send_buf_p/u_recv_buf_p via ShmBufferMalloc
      // (Stencil.h:1090) and reaches a peer's buffer with ShmBufferTranslate
      // (Stencil.h:423) -- that IPC mapping is what makes an intranode
      // transfer intranode.  Unshared device memory cannot be translated, so
      // the peer path never engages and every packet falls back to MPI, which
      // would silently make this benchmark measure the wrong thing.
      Grid.ShmBufferFreeAll();
      std::vector<HalfSpinColourVectorD *> xbuf(8), rbuf(8);
      for(int d=0;d<8;d++){
	xbuf[d] = (HalfSpinColourVectorD *)Grid.ShmBufferMalloc(bytes);
	rbuf[d] = (HalfSpinColourVectorD *)Grid.ShmBufferMalloc(bytes);
	GRID_ASSERT(xbuf[d]!=nullptr);
	GRID_ASSERT(rbuf[d]!=nullptr);
      }

      // Packet table in the same shape the Stencil builds: 8 directions,
      // to/from ranks from ShiftedRanks, participation gated on mpi_layout.
      // COMPACT packet list, exactly as the Stencil builds it: only the
      // directions that actually communicate get a packet.  Passing a
      // self-send with do_send=0 is not equivalent -- the assert
      // "dest != _processor" in StencilSendToRecvFromBegin is unconditional.
      std::vector<int> to_rank, from_rank, buf_id;
      double offered = 0.0, offnode = 0.0;
      for(int mu=0;mu<4;mu++){
	if ( mpi_layout[mu] > 1 ) {
	  int s,d;
	  Grid.ShiftedRanks(mu,+1,s,d); to_rank.push_back(d); from_rank.push_back(s); buf_id.push_back(mu);
	  Grid.ShiftedRanks(mu,-1,s,d); to_rank.push_back(d); from_rank.push_back(s); buf_id.push_back(mu+4);
	}
      }
      int npkt = to_rank.size();
      for(int i=0;i<npkt;i++){
	offered += 2.0*bytes;                                  // send + receive
	if ( Grid.IsOffNode(to_rank[i]) ) offnode += 2.0*bytes;
      }

      double mpibytes = 0.0, t = 0.0;
      std::vector<CommsRequest_t> reqs;

      for(int n=0;n<Nloop;n++){
	reqs.resize(0);
	double mb = 0.0;
	Grid.StencilBarrier();
	double t0=usecond();
	for(int i=0;i<npkt;i++){
	  mb += Grid.StencilSendToRecvFromPrepare(reqs,
			 (void *)xbuf[buf_id[i]], to_rank[i],   1,
			 (void *)rbuf[buf_id[i]], from_rank[i], 1,
			 bytes,bytes,i);
	}
	Grid.StencilSendToRecvFromPollDtoH(reqs); /* Starts MPI */
	acceleratorCopySynchronise();
	for(int i=0;i<npkt;i++){
	  // xmit and xmit_comp identical: no compression in this benchmark
	  mb += Grid.StencilSendToRecvFromBegin(reqs,
			 (void *)xbuf[buf_id[i]],(void *)xbuf[buf_id[i]], to_rank[i],   1,
			 (void *)rbuf[buf_id[i]],(void *)rbuf[buf_id[i]], from_rank[i], 1,
			 bytes,bytes,i);
	}
	Grid.StencilSendToRecvFromPollIRecv(reqs);
	Grid.StencilSendToRecvFromComplete(reqs,0);
	t += usecond()-t0;
	mpibytes += mb;
      }

      // Sum bytes over the machine, average the time over ranks, report per
      // node.  Reduced once, outside the timed region.  bytes/us == MB/s.
      double sumMpi = mpibytes;        Grid.GlobalSum(sumMpi);
      double sumOff = offnode*Nloop;   Grid.GlobalSum(sumOff);
      double sumAll = offered*Nloop;   Grid.GlobalSum(sumAll);
      double sumT   = t;               Grid.GlobalSum(sumT);
      double avgT   = sumT/Nrank;

      std::cout<<GridLogMessage
	       << std::setw(4) << lat <<" "<< std::setw(4) << Ls
	       <<"  "<< std::setw(9) << sumAll/Nloop/Nnode/1.0e6
	       <<"  "<< std::setw(9) << sumOff/Nloop/Nnode/1.0e6
	       <<"  "<< std::setw(9) << sumMpi/Nloop/Nnode/1.0e6
	       <<"  "<< std::setw(9) << avgT/Nloop
	       <<"  "<< std::setw(8) << sumOff/avgT/Nnode/1000.0
	       <<"  "<< std::setw(8) << sumMpi/avgT/Nnode/1000.0
	       <<std::endl;

      Grid.ShmBufferFreeAll();   // bump allocator: released as a whole
    }
  }

  std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
  std::cout<<GridLogMessage << "= All done; Bye Bye"<<std::endl;
  std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;

  Grid_finalize();
}
