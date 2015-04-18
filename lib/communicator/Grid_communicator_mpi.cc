#include "Grid.h"
#include <mpi.h>

namespace Grid {

  // Should error check all MPI calls.

CartesianCommunicator::CartesianCommunicator(std::vector<int> &processors)
{
  _ndimension = processors.size();
  std::vector<int> periodic(_ndimension,1);

  _Nprocessors=1;
  _processors = processors;
  _processor_coor.resize(_ndimension);
  
  MPI_Cart_create(MPI_COMM_WORLD, _ndimension,&_processors[0],&periodic[0],1,&communicator);
  MPI_Comm_rank(communicator,&_processor);
  MPI_Cart_coords(communicator,_processor,_ndimension,&_processor_coor[0]);

  for(int i=0;i<_ndimension;i++){
    _Nprocessors*=_processors[i];
  }
  
  int Size; 
  MPI_Comm_size(communicator,&Size);
  
  assert(Size==_Nprocessors);
}

void CartesianCommunicator::GlobalSum(float &f){
  MPI_Allreduce(MPI_IN_PLACE,&f,1,MPI_FLOAT,MPI_SUM,communicator);
}
void CartesianCommunicator::GlobalSumVector(float *f,int N)
{
  MPI_Allreduce(MPI_IN_PLACE,f,N,MPI_FLOAT,MPI_SUM,communicator);
}
void CartesianCommunicator::GlobalSum(double &d)
{
  MPI_Allreduce(MPI_IN_PLACE,&d,1,MPI_DOUBLE,MPI_SUM,communicator);
}
void CartesianCommunicator::GlobalSumVector(double *d,int N)
{
  MPI_Allreduce(MPI_IN_PLACE,d,N,MPI_DOUBLE,MPI_SUM,communicator);
}
void CartesianCommunicator::ShiftedRanks(int dim,int shift,int &source,int &dest)
{
  MPI_Cart_shift(communicator,dim,shift,&source,&dest);
}
int CartesianCommunicator::RankFromProcessorCoor(std::vector<int> &coor)
{
  int rank;
  MPI_Cart_rank  (communicator, &coor[0], &rank);
  return rank;
}
void  CartesianCommunicator::ProcessorCoorFromRank(int rank, std::vector<int> &coor)
{
  coor.resize(_ndimension);
  MPI_Cart_coords  (communicator, rank, _ndimension,&coor[0]);
}

// Basic Halo comms primitive
void CartesianCommunicator::SendToRecvFrom(void *xmit,
					   int dest,
					   void *recv,
					   int from,
					   int bytes)
{
  MPI_Request reqs[2];
  MPI_Status OkeyDokey[2];
  int rank = _processor;
  MPI_Isend(xmit, bytes, MPI_CHAR,dest,_processor,communicator,&reqs[0]);
  MPI_Irecv(recv, bytes, MPI_CHAR,from,from,communicator,&reqs[1]);
  MPI_Waitall(2,reqs,OkeyDokey);

}

void CartesianCommunicator::Barrier(void)
{
  MPI_Barrier(communicator);
}

void CartesianCommunicator::Broadcast(int root,void* data, int bytes)
{
  MPI_Bcast(data,
	    bytes,
	    MPI_BYTE,
	    root,
	    communicator);
}

}

