#include <Grid.h>

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  std::vector<int> simd_layout = GridDefaultSimd();
  std::vector<int> mpi_layout  = GridDefaultMpi();

  int Nloop=10;
  int nmu=0;
  for(int mu=0;mu<4;mu++) if (mpi_layout[mu]>1) nmu++;
  std::cout << "===================================================================================================="<<std::endl;
  std::cout << "= Benchmarking concurrent halo exchange in "<<nmu<<" dimensions"<<std::endl;
  std::cout << "===================================================================================================="<<std::endl;
  std::cout << "  L  "<<"\t\t"<<" Ls  "<<"\t\t"<<"bytes"<<"\t\t"<<"MB/s uni"<<"\t\t"<<"MB/s bidi"<<std::endl;



  for(int lat=4;lat<=16;lat+=4){
    for(int Ls=1;Ls<=16;Ls*=2){

      std::vector<int> latt_size  ({lat,lat,lat,lat});

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

      double xbytes    = Nloop*bytes*2*ncomm;
      double rbytes    = xbytes;
      double bidibytes = xbytes+rbytes;

      double time = stop-start;

      std::cout << lat<<"\t\t"<<Ls<<"\t\t"<<bytes<<"\t\t"<<xbytes/time<<"\t\t"<<bidibytes/time<<std::endl;
    }
  }    


  std::cout << "===================================================================================================="<<std::endl;
  std::cout << "= Benchmarking sequential halo exchange in "<<nmu<<" dimensions"<<std::endl;
  std::cout << "===================================================================================================="<<std::endl;
  std::cout << "  L  "<<"\t\t"<<" Ls  "<<"\t\t"<<"bytes"<<"\t\t"<<"MB/s uni"<<"\t\t"<<"MB/s bidi"<<std::endl;



  for(int lat=4;lat<=16;lat+=4){
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

      double xbytes    = Nloop*bytes*2*ncomm;
      double rbytes    = xbytes;
      double bidibytes = xbytes+rbytes;

      double time = stop-start;

      std::cout << lat<<"\t\t"<<Ls<<"\t\t"<<bytes<<"\t\t"<<xbytes/time<<"\t\t"<<bidibytes/time<<std::endl;
    }
  }  



  Grid_finalize();
}
