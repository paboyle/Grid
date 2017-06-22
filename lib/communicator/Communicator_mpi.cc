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

CartesianCommunicator::CartesianCommunicator(const std::vector<int> &processors) 
{
  InitFromMPICommunicator(processors,communicator_world);
  std::cout << "Passed communicator world to a new communicator" <<std::endl;
}
CartesianCommunicator::CartesianCommunicator(const std::vector<int> &processors,const CartesianCommunicator &parent) 
{
  _ndimension = processors.size();
  assert(_ndimension = parent._ndimension);

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  // split the communicator
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  std::vector<int> ratio(_ndimension);
  std::vector<int> rcoor(_ndimension);
  std::vector<int> scoor(_ndimension);

  int Nsubcomm=1;
  int Nsubrank=1;
  for(int d=0;d<_ndimension;d++) {
    ratio[d] = parent._processors[d] / processors[d];
    rcoor[d] = parent._processor_coor[d] / processors[d];
    scoor[d] = parent._processor_coor[d] % processors[d];
    assert(ratio[d] * processors[d] == parent._processors[d]); // must exactly subdivide
    Nsubcomm *= ratio[d];
    Nsubrank *= processors[d];
  }

  int rlex, slex;
  Lexicographic::IndexFromCoor(rcoor,rlex,ratio);
  Lexicographic::IndexFromCoor(scoor,slex,processors);

  MPI_Comm comm_split;
  if ( Nsubcomm > 1 ) { 
    int ierr= MPI_Comm_split(communicator_world, rlex, slex,&comm_split);
    assert(ierr==0);
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    // Declare victory
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    std::cout << "Divided communicator "<< parent._Nprocessors<<" into "
	      <<Nsubcomm <<" communicators with " << Nsubrank << " ranks"<<std::endl;
  } else {
    comm_split = communicator_world;
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  // Set up from the new split communicator
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  InitFromMPICommunicator(processors,comm_split);

}

//////////////////////////////////////////////////////////////////////////////////////////////////////
// Take an MPI_Comm and self assemble
//////////////////////////////////////////////////////////////////////////////////////////////////////
void CartesianCommunicator::InitFromMPICommunicator(const std::vector<int> &processors, MPI_Comm communicator_base)
{
  if ( communicator_base != communicator_world ) {
    std::cout << "Cartesian communicator created with a non-world communicator"<<std::endl;
  }
  _ndimension = processors.size();
  _processor_coor.resize(_ndimension);

  /////////////////////////////////
  // Count the requested nodes
  /////////////////////////////////
  _Nprocessors=1;
  _processors = processors;
  for(int i=0;i<_ndimension;i++){
    _Nprocessors*=_processors[i];
  }

  std::vector<int> periodic(_ndimension,1);
  MPI_Cart_create(communicator_base, _ndimension,&_processors[0],&periodic[0],1,&communicator);
  MPI_Comm_rank(communicator,&_processor);
  MPI_Cart_coords(communicator,_processor,_ndimension,&_processor_coor[0]);

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

}

