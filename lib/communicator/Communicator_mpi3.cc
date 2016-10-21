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

// Global used by Init and nowhere else. How to hide?
int Rank(void) {
  int pe;
  MPI_Comm_rank(MPI_COMM_WORLD,&pe);
  return pe;
}
  // Should error check all MPI calls.
void CartesianCommunicator::Init(int *argc, char ***argv) {
  int flag;
  MPI_Initialized(&flag); // needed to coexist with other libs apparently
  if ( !flag ) {
    MPI_Init(argc,argv);
  }
}
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Want to implement some magic ... Group sub-cubes into those on same node
  //
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////

void CartesianCommunicator::ShiftedRanks(int dim,int shift,int &source,int &dest)
{
  std::vector<int> coor = _processor_coor;

  assert(std::abs(shift) <_processors[dim]);

  coor[dim] = (_processor_coor[dim] + shift + _processors[dim])%_processors[dim];
  Lexicographic::IndexFromCoor(coor,source,_processors);
  source = LexicographicToWorldRank[source];

  coor[dim] = (_processor_coor[dim] - shift + _processors[dim])%_processors[dim];
  Lexicographic::IndexFromCoor(coor,dest,_processors);
  dest = LexicographicToWorldRank[dest];
}
int CartesianCommunicator::RankFromProcessorCoor(std::vector<int> &coor)
{
  int rank;
  Lexicographic::IndexFromCoor(coor,rank,_processors);
  rank = LexicographicToWorldRank[rank];
  return rank;
}
void  CartesianCommunicator::ProcessorCoorFromRank(int rank, std::vector<int> &coor)
{
  Lexicographic::CoorFromIndex(coor,rank,_processors);
  rank = LexicographicToWorldRank[rank];
}

CartesianCommunicator::CartesianCommunicator(const std::vector<int> &processors)
{ 
  _ndimension = processors.size();
  std::cout << "Creating "<< _ndimension << " dim communicator "<<std::endl;
  for(int d =0;d<_ndimension;d++){
    std::cout << processors[d]<<" ";
  };
  std::cout << std::endl;

  WorldDims = processors;

  communicator = MPI_COMM_WORLD;
  MPI_Comm shmcomm;
  MPI_Comm_split_type(communicator, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL,&shmcomm);
  MPI_Comm_rank(communicator,&WorldRank);
  MPI_Comm_size(communicator,&WorldSize);
  MPI_Comm_rank(shmcomm     ,&ShmRank);
  MPI_Comm_size(shmcomm     ,&ShmSize);
  GroupSize = WorldSize/ShmSize;

  std::cout<< "Ranks per node "<< ShmSize << std::endl;
  std::cout<< "Nodes          "<< GroupSize << std::endl;
  std::cout<< "Ranks          "<< WorldSize << std::endl;
  
  ////////////////////////////////////////////////////////////////
  // Assert power of two shm_size.
  ////////////////////////////////////////////////////////////////
  int log2size = -1;
  for(int i=0;i<=MAXLOG2RANKSPERNODE;i++){  
    if ( (0x1<<i) == ShmSize ) {
      log2size = i;
      break;
    }
  }
  assert(log2size != -1);
  
  ////////////////////////////////////////////////////////////////
  // Identify subblock of ranks on node spreading across dims
  // in a maximally symmetrical way
  ////////////////////////////////////////////////////////////////
  int dim = 0;
  
  ShmDims.resize(_ndimension,1);
  GroupDims.resize(_ndimension);
  
  ShmCoor.resize(_ndimension);
  GroupCoor.resize(_ndimension);
  WorldCoor.resize(_ndimension);
  for(int l2=0;l2<log2size;l2++){
    while ( WorldDims[dim] / ShmDims[dim] <= 1 ) dim=(dim+1)%_ndimension;
    ShmDims[dim]*=2;
    dim=(dim+1)%_ndimension;
  }
  
  std::cout << "Shm group dims "<<std::endl;
  for(int d =0;d<_ndimension;d++){
    std::cout << ShmDims[d]<<" ";
  };
  std::cout << std::endl;

  ////////////////////////////////////////////////////////////////
  // Establish torus of processes and nodes with sub-blockings
  ////////////////////////////////////////////////////////////////
  for(int d=0;d<_ndimension;d++){
    GroupDims[d] = WorldDims[d]/ShmDims[d];
  }
  std::cout << "Group dims "<<std::endl;
  for(int d =0;d<_ndimension;d++){
    std::cout << GroupDims[d]<<" ";
  };
  std::cout << std::endl;
  
  MPI_Group WorldGroup, ShmGroup;
  MPI_Comm_group (communicator, &WorldGroup); 
  MPI_Comm_group (shmcomm, &ShmGroup);
  
  std::vector<int> world_ranks(WorldSize); 
  std::vector<int> group_ranks(WorldSize); 
  std::vector<int> mygroup(GroupSize);
  for(int r=0;r<WorldSize;r++) world_ranks[r]=r;
  
  MPI_Group_translate_ranks (WorldGroup,WorldSize,&world_ranks[0],ShmGroup, &group_ranks[0]); 

  ////////////////////////////////////////////////////////////////
  // Check processor counts match
  ////////////////////////////////////////////////////////////////
  _Nprocessors=1;
  _processors = processors;
  _processor_coor.resize(_ndimension);
  for(int i=0;i<_ndimension;i++){
    std::cout << " p " << _processors[i]<<std::endl;
    _Nprocessors*=_processors[i];
  }
  std::cout << " World " <<WorldSize <<" Nproc "<<_Nprocessors<<std::endl;
  assert(WorldSize==_Nprocessors);
  
  ///////////////////////////////////////////////////////////////////
  // Identify who is in my group and noninate the leader
  ///////////////////////////////////////////////////////////////////
  int g=0;
  for(int rank=0;rank<WorldSize;rank++){
    if(group_ranks[rank]!=MPI_UNDEFINED){
	  mygroup[g] = rank;
    }
  }
  
  std::sort(mygroup.begin(),mygroup.end(),std::greater<int>());
  int myleader = mygroup[0];
  
  std::vector<int> leaders_1hot(WorldSize,0);
  std::vector<int> leaders_group(GroupSize,0);
  leaders_1hot [ myleader ] = 1;
      
  ///////////////////////////////////////////////////////////////////
  // global sum leaders over comm world
  ///////////////////////////////////////////////////////////////////
  int ierr=MPI_Allreduce(MPI_IN_PLACE,&leaders_1hot[0],WorldSize,MPI_INT,MPI_SUM,communicator);
  assert(ierr==0);
  
  ///////////////////////////////////////////////////////////////////
  // find the group leaders world rank
  ///////////////////////////////////////////////////////////////////
  int group=0;
  for(int l=0;l<WorldSize;l++){
    if(leaders_1hot[l]){
      leaders_group[group++] = l;
    }
  }
  
  ///////////////////////////////////////////////////////////////////
  // Identify the rank of the group in which I (and my leader) live
  ///////////////////////////////////////////////////////////////////
  GroupRank=-1;
  for(int g=0;g<GroupSize;g++){
    if (myleader == leaders_group[g]){
      GroupRank=g;
    }
  }
  assert(GroupRank!=-1);
      
  ////////////////////////////////////////////////////////////////
  // Establish mapping between lexico physics coord and WorldRank
  // 
  ////////////////////////////////////////////////////////////////
  LexicographicToWorldRank.resize(WorldSize,0);
  Lexicographic::CoorFromIndex(GroupCoor,GroupRank,GroupDims);
  Lexicographic::CoorFromIndex(ShmCoor,ShmRank,ShmDims);
  for(int d=0;d<_ndimension;d++){
    WorldCoor[d] = GroupCoor[d]*ShmDims[d]+ShmCoor[d];
  }
  _processor_coor = WorldCoor;

  int lexico;
  Lexicographic::IndexFromCoor(WorldCoor,lexico,WorldDims);
  LexicographicToWorldRank[lexico]=WorldRank;
  _processor = lexico;

  ///////////////////////////////////////////////////////////////////
  // global sum Lexico to World mapping
  ///////////////////////////////////////////////////////////////////
  ierr=MPI_Allreduce(MPI_IN_PLACE,&LexicographicToWorldRank[0],WorldSize,MPI_INT,MPI_SUM,communicator);
  assert(ierr==0);
  
};

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

