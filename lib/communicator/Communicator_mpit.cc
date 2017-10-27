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
#include <Grid/GridCore.h>
#include <Grid/GridQCDcore.h>
#include <Grid/qcd/action/ActionCore.h>
#include <mpi.h>

namespace Grid {


///////////////////////////////////////////////////////////////////////////////////////////////////
// Info that is setup once and indept of cartesian layout
///////////////////////////////////////////////////////////////////////////////////////////////////
MPI_Comm CartesianCommunicator::communicator_world;

// Should error check all MPI calls.
void CartesianCommunicator::Init(int *argc, char ***argv) {
  int flag;
  int provided;
  MPI_Initialized(&flag); // needed to coexist with other libs apparently
  if ( !flag ) {
    MPI_Init_thread(argc,argv,MPI_THREAD_MULTIPLE,&provided);
    if ( provided != MPI_THREAD_MULTIPLE ) {
      QCD::WilsonKernelsStatic::Comms = QCD::WilsonKernelsStatic::CommsThenCompute;
    }
  }
  MPI_Comm_dup (MPI_COMM_WORLD,&communicator_world);
  ShmInitGeneric();
}

CartesianCommunicator::~CartesianCommunicator()
{
  if (communicator && !MPI::Is_finalized())
    MPI_Comm_free(&communicator);
}


void CartesianCommunicator::GlobalSum(uint32_t &u){
  int ierr=MPI_Allreduce(MPI_IN_PLACE,&u,1,MPI_UINT32_T,MPI_SUM,communicator);
  assert(ierr==0);
}
void CartesianCommunicator::GlobalSum(uint64_t &u){
  int ierr=MPI_Allreduce(MPI_IN_PLACE,&u,1,MPI_UINT64_T,MPI_SUM,communicator);
  assert(ierr==0);
}
void CartesianCommunicator::GlobalXOR(uint32_t &u){
  int ierr=MPI_Allreduce(MPI_IN_PLACE,&u,1,MPI_UINT32_T,MPI_BXOR,communicator);
  assert(ierr==0);
}
void CartesianCommunicator::GlobalXOR(uint64_t &u){
  int ierr=MPI_Allreduce(MPI_IN_PLACE,&u,1,MPI_UINT64_T,MPI_BXOR,communicator);
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
  int myrank = _processor;
  int ierr;
  if ( CommunicatorPolicy == CommunicatorPolicyConcurrent ) { 
    MPI_Request xrq;
    MPI_Request rrq;

    ierr =MPI_Irecv(recv, bytes, MPI_CHAR,from,from,communicator,&rrq);
    ierr|=MPI_Isend(xmit, bytes, MPI_CHAR,dest,_processor,communicator,&xrq);
    
    assert(ierr==0);
    list.push_back(xrq);
    list.push_back(rrq);
  } else { 
    // Give the CPU to MPI immediately; can use threads to overlap optionally
    ierr=MPI_Sendrecv(xmit,bytes,MPI_CHAR,dest,myrank,
		      recv,bytes,MPI_CHAR,from, from,
		      communicator,MPI_STATUS_IGNORE);
    assert(ierr==0);
  }
}
void CartesianCommunicator::SendToRecvFromComplete(std::vector<CommsRequest_t> &list)
{
  if ( CommunicatorPolicy == CommunicatorPolicyConcurrent ) { 
    int nreq=list.size();
    std::vector<MPI_Status> status(nreq);
    int ierr = MPI_Waitall(nreq,&list[0],&status[0]);
    assert(ierr==0);
  }
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
  ///////////////////////////////////////////////////////
  // Should only be used prior to Grid Init finished.
  // Check for this?
  ///////////////////////////////////////////////////////
int CartesianCommunicator::RankWorld(void){ 
  int r; 
  MPI_Comm_rank(communicator_world,&r);
  return r;
}
void CartesianCommunicator::BroadcastWorld(int root,void* data, int bytes)
{
  int ierr= MPI_Bcast(data,
		      bytes,
		      MPI_BYTE,
		      root,
		      communicator_world);
  assert(ierr==0);
}

double CartesianCommunicator::StencilSendToRecvFromBegin(std::vector<CommsRequest_t> &list,
							 void *xmit,
							 int xmit_to_rank,
							 void *recv,
							 int recv_from_rank,
							 int bytes,int dir)
{
  int myrank = _processor;
  int ierr;
  int ncomm  =communicator_halo.size(); 
  int commdir=dir%ncomm;
  
  //  std::cout << " sending on communicator "<<dir<<" " <<communicator_halo[dir]<<std::endl;
  // Give the CPU to MPI immediately; can use threads to overlap optionally
  MPI_Request req[2];
  MPI_Irecv(recv,bytes,MPI_CHAR,recv_from_rank,recv_from_rank, communicator_halo[commdir],&req[1]);
  MPI_Isend(xmit,bytes,MPI_CHAR,xmit_to_rank  ,myrank        , communicator_halo[commdir],&req[0]);

  list.push_back(req[0]);
  list.push_back(req[1]);
  return 2.0*bytes;
}
void CartesianCommunicator::StencilSendToRecvFromComplete(std::vector<CommsRequest_t> &waitall,int dir)
{ 
  int nreq=waitall.size();
  MPI_Waitall(nreq, &waitall[0], MPI_STATUSES_IGNORE);
};
double CartesianCommunicator::StencilSendToRecvFrom(void *xmit,
						    int xmit_to_rank,
						    void *recv,
						    int recv_from_rank,
						    int bytes,int dir)
{
  int myrank = _processor;
  int ierr;
  //  std::cout << " sending on communicator "<<dir<<" " <<communicator_halo.size()<< <std::endl;

  int ncomm  =communicator_halo.size(); 
  int commdir=dir%ncomm;
  // Give the CPU to MPI immediately; can use threads to overlap optionally
  MPI_Request req[2];
  MPI_Irecv(recv,bytes,MPI_CHAR,recv_from_rank,recv_from_rank, communicator_halo[commdir],&req[1]);
  MPI_Isend(xmit,bytes,MPI_CHAR,xmit_to_rank  ,myrank        , communicator_halo[commdir],&req[0]);
  MPI_Waitall(2, req, MPI_STATUSES_IGNORE);
  return 2.0*bytes;
}



}

