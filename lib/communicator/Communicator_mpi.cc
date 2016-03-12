    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/communicator/Communicator_mpi.cc

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
#include "Grid.h"
#include <mpi.h>

namespace Grid {

  // Should error check all MPI calls.
void CartesianCommunicator::Init(int *argc, char ***argv) {
  int flag;
  MPI_Initialized(&flag); // needed to coexist with other libs apparently
  if ( !flag ) {
    MPI_Init(argc,argv);
  }
}

  int Rank(void) {
    int pe;
    MPI_Comm_rank(MPI_COMM_WORLD,&pe);
    return pe;
  }

CartesianCommunicator::CartesianCommunicator(const std::vector<int> &processors)
{
  _ndimension = processors.size();
  std::vector<int> periodic(_ndimension,1);

  _Nprocessors=1;
  _processors = processors;
  _processor_coor.resize(_ndimension);
  std::cout << processors << std::endl;
  
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

void CartesianCommunicator::GlobalSum(uint32_t &u){
  int ierr=MPI_Allreduce(MPI_IN_PLACE,&u,1,MPI_UINT32_T,MPI_SUM,communicator);
  assert(ierr==0);
}
void CartesianCommunicator::GlobalSum(uint64_t &u){
  int ierr=MPI_Allreduce(MPI_IN_PLACE,&u,1,MPI_UINT64_T,MPI_SUM,communicator);
  assert(ierr==0);
}
void CartesianCommunicator::GlobalSum(float &f){
  int ierr=MPI_Allreduce(MPI_IN_PLACE,&f,1,MPI_FLOAT,MPI_SUM,communicator);
  assert(ierr==0);
}
void CartesianCommunicator::GlobalSumVector(float *f,int N)
{
  int ierr=MPI_Allreduce(MPI_IN_PLACE,f,N,MPI_FLOAT,MPI_SUM,communicator);
  assert(ierr==0);
}
void CartesianCommunicator::GlobalSum(double &d)
{
  int ierr = MPI_Allreduce(MPI_IN_PLACE,&d,1,MPI_DOUBLE,MPI_SUM,communicator);
  assert(ierr==0);
}
void CartesianCommunicator::GlobalSumVector(double *d,int N)
{
  int ierr = MPI_Allreduce(MPI_IN_PLACE,d,N,MPI_DOUBLE,MPI_SUM,communicator);
  assert(ierr==0);
}
void CartesianCommunicator::ShiftedRanks(int dim,int shift,int &source,int &dest)
{
  int ierr=MPI_Cart_shift(communicator,dim,shift,&source,&dest);
  assert(ierr==0);
}
int CartesianCommunicator::RankFromProcessorCoor(std::vector<int> &coor)
{
  int rank;
  int ierr=MPI_Cart_rank  (communicator, &coor[0], &rank);
  assert(ierr==0);
  return rank;
}
void  CartesianCommunicator::ProcessorCoorFromRank(int rank, std::vector<int> &coor)
{
  coor.resize(_ndimension);
  int ierr=MPI_Cart_coords  (communicator, rank, _ndimension,&coor[0]);
  assert(ierr==0);
}

// Basic Halo comms primitive
void CartesianCommunicator::SendToRecvFrom(void *xmit,
					   int dest,
					   void *recv,
					   int from,
					   int bytes)
{
  std::vector<CommsRequest_t> reqs(0);
  SendToRecvFromBegin(reqs,xmit,dest,recv,from,bytes);
  SendToRecvFromComplete(reqs);
}

void CartesianCommunicator::SendRecvPacket(void *xmit,
					   void *recv,
					   int sender,
					   int receiver,
					   int bytes)
{
  MPI_Status stat;
  assert(sender != receiver);
  int tag = sender;
  if ( _processor == sender ) {
    MPI_Send(xmit, bytes, MPI_CHAR,receiver,tag,communicator);
  }
  if ( _processor == receiver ) { 
    MPI_Recv(recv, bytes, MPI_CHAR,sender,tag,communicator,&stat);
  }
}

// Basic Halo comms primitive
void CartesianCommunicator::SendToRecvFromBegin(std::vector<CommsRequest_t> &list,
						void *xmit,
						int dest,
						void *recv,
						int from,
						int bytes)
{
  MPI_Request xrq;
  MPI_Request rrq;
  int rank = _processor;
  int ierr;
  ierr =MPI_Isend(xmit, bytes, MPI_CHAR,dest,_processor,communicator,&xrq);
  ierr|=MPI_Irecv(recv, bytes, MPI_CHAR,from,from,communicator,&rrq);
  
  assert(ierr==0);

  list.push_back(xrq);
  list.push_back(rrq);
}
void CartesianCommunicator::SendToRecvFromComplete(std::vector<CommsRequest_t> &list)
{
  int nreq=list.size();
  std::vector<MPI_Status> status(nreq);
  int ierr = MPI_Waitall(nreq,&list[0],&status[0]);

  assert(ierr==0);
}

void CartesianCommunicator::Barrier(void)
{
  int ierr = MPI_Barrier(communicator);
  assert(ierr==0);
}

void CartesianCommunicator::Broadcast(int root,void* data, int bytes)
{
  int ierr=MPI_Bcast(data,
		     bytes,
		     MPI_BYTE,
		     root,
		     communicator);
  assert(ierr==0);
}

void CartesianCommunicator::BroadcastWorld(int root,void* data, int bytes)
{
  int ierr= MPI_Bcast(data,
		      bytes,
		      MPI_BYTE,
		      root,
		      MPI_COMM_WORLD);
  assert(ierr==0);
}

}

