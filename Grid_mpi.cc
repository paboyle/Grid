#include "Grid.h"
#include <mpi.h>

namespace dpo {

  // Should error check all MPI calls.

CartesianCommunicator::CartesianCommunicator(std::vector<int> &processors)
{
  _ndimension = _processors.size();
  std::vector<int> periodic(_ndimension,1);

  _processors = processors;
  _processor_coords.resize(_ndimension);

  MPI_Cart_create(MPI_COMM_WORLD _ndimension,&_processors[0],&periodic[0],1,&communicator);
  MPI_Comm_rank(communicator,&_processor);
  MPI_Cart_coords(communicator,_processor,_ndimension,&_processor_coords[0]);

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

int CartesianCommunicator::ShiftedRank(int dim,int shift)
{
  int rank;
  MPI_Cart_shift(communicator,dim,shift,&_processor,&rank);
  return rank;
}
int CartesianCommunicator::Rank(std::vector<int> coor)
{
  int rank;
  MPI_Cart_rank  (communicator, &coor[0], &rank);
  return rank;
}

// Basic Halo comms primitive
void CartesianCommunicator::SendToRecvFrom(void *xmit,
					   std::vector<int> &to_coordinate,
					   void *recv,
					   std::vector<int> &from_coordinate,
					   int bytes)
{
  int dest = Rank(to_coordinate);
  int from = Rank(from_coordinate);

  MPI_Request x_req;
  MPI_Request r_req;
  MPI_Status OkeyDokey;

  MPI_Isend(xmit, bytes, MPI_CHAR,dest,dest,communicator,&x_req);
  MPI_Irecv(recv, bytes, MPI_CHAR,from,from,communicator,&r_req);
  MPI_Wait(&x_req,&OkeyDokey);
  MPI_Wait(&r_req,&OkeyDokey);
  MPI_Barrier();
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
