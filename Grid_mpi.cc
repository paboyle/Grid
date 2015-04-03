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
  printf("Hello world from processor [");
  for(int i=0;i<_ndimension;i++){
    printf("%d ",_processor_coor[i]);
    _Nprocessors*=_processors[i];
  }
  printf("]\n");
  fflush(stdout);
}

void CartesianCommunicator::GlobalSumF(float &f){
  MPI_Allreduce(&f,&f,1,MPI_FLOAT,MPI_SUM,communicator);
}
void CartesianCommunicator::GlobalSumFVector(float *f,int N)
{
  MPI_Allreduce(f,f,N,MPI_FLOAT,MPI_SUM,communicator);
}
void CartesianCommunicator::GlobalSumF(double &d)
{
  MPI_Allreduce(&d,&d,1,MPI_DOUBLE,MPI_SUM,communicator);
}
void CartesianCommunicator::GlobalSumFVector(double *d,int N)
{
  MPI_Allreduce(d,d,N,MPI_DOUBLE,MPI_SUM,communicator);
}

void CartesianCommunicator::ShiftedRanks(int dim,int shift,int &source,int &dest)
{
  MPI_Cart_shift(communicator,dim,shift,&source,&dest);
}
int CartesianCommunicator::Rank(std::vector<int> coor)
{
  int rank;
  MPI_Cart_rank  (communicator, &coor[0], &rank);
  return rank;
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

}

#if 0

// Could possibly do a direct block strided send?
  int MPI_Type_vector(
		      int count,
		      int blocklength,
		      int stride,
		      MPI_Datatype old_type,
  MPI_Datatype *newtype_p
		      );

#endif
