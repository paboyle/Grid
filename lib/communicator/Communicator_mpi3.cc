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
#include <Grid/Grid.h>
#include <mpi.h>

namespace Grid {

///////////////////////////////////////////////////////////////////////////////////////////////////
// Info that is setup once and indept of cartesian layout
///////////////////////////////////////////////////////////////////////////////////////////////////
int CartesianCommunicator::ShmSetup = 0;

int CartesianCommunicator::ShmRank;
int CartesianCommunicator::ShmSize;
int CartesianCommunicator::GroupRank;
int CartesianCommunicator::GroupSize;
int CartesianCommunicator::WorldRank;
int CartesianCommunicator::WorldSize;

MPI_Comm CartesianCommunicator::communicator_world;
MPI_Comm CartesianCommunicator::ShmComm;
MPI_Win  CartesianCommunicator::ShmWindow;

std::vector<int> CartesianCommunicator::GroupRanks;  
std::vector<int> CartesianCommunicator::MyGroup;
std::vector<void *> CartesianCommunicator::ShmCommBufs;

void *CartesianCommunicator::ShmBufferSelf(void)
{
  return ShmCommBufs[ShmRank];
}
void *CartesianCommunicator::ShmBuffer(int rank)
{
  int gpeer = GroupRanks[rank];
  if (gpeer == MPI_UNDEFINED){
    return NULL;
  } else { 
    return ShmCommBufs[gpeer];
  }
}
void *CartesianCommunicator::ShmBufferTranslate(int rank,void * local_p)
{
  int gpeer = GroupRanks[rank];
  if (gpeer == MPI_UNDEFINED){
    return NULL;
  } else { 
    uint64_t offset = (uint64_t)local_p - (uint64_t)ShmCommBufs[ShmRank];
    uint64_t remote = (uint64_t)ShmCommBufs[gpeer]+offset;
    return (void *) remote;
  }
}

void CartesianCommunicator::Init(int *argc, char ***argv) {
  int flag;
  MPI_Initialized(&flag); // needed to coexist with other libs apparently
  if ( !flag ) {
    MPI_Init(argc,argv);
  }

  MPI_Comm_dup (MPI_COMM_WORLD,&communicator_world);
  MPI_Comm_rank(communicator_world,&WorldRank);
  MPI_Comm_size(communicator_world,&WorldSize);

  /////////////////////////////////////////////////////////////////////
  // Split into groups that can share memory
  /////////////////////////////////////////////////////////////////////
  MPI_Comm_split_type(communicator_world, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL,&ShmComm);
  MPI_Comm_rank(ShmComm     ,&ShmRank);
  MPI_Comm_size(ShmComm     ,&ShmSize);
  GroupSize = WorldSize/ShmSize;

  /////////////////////////////////////////////////////////////////////
  // find world ranks in our SHM group (i.e. which ranks are on our node)
  /////////////////////////////////////////////////////////////////////
  MPI_Group WorldGroup, ShmGroup;
  MPI_Comm_group (communicator_world, &WorldGroup); 
  MPI_Comm_group (ShmComm, &ShmGroup);
  
  std::vector<int> world_ranks(WorldSize); 
  GroupRanks.resize(WorldSize); 
  for(int r=0;r<WorldSize;r++) world_ranks[r]=r;
  
  MPI_Group_translate_ranks (WorldGroup,WorldSize,&world_ranks[0],ShmGroup, &GroupRanks[0]); 

  ///////////////////////////////////////////////////////////////////
  // Identify who is in my group and noninate the leader
  ///////////////////////////////////////////////////////////////////
  int g=0;
  MyGroup.resize(ShmSize);
  for(int rank=0;rank<WorldSize;rank++){
    if(GroupRanks[rank]!=MPI_UNDEFINED){
      assert(g<ShmSize);
      MyGroup[g++] = rank;
    }
  }
  
  std::sort(MyGroup.begin(),MyGroup.end(),std::less<int>());
  int myleader = MyGroup[0];
  
  std::vector<int> leaders_1hot(WorldSize,0);
  std::vector<int> leaders_group(GroupSize,0);
  leaders_1hot [ myleader ] = 1;
    
  ///////////////////////////////////////////////////////////////////
  // global sum leaders over comm world
  ///////////////////////////////////////////////////////////////////
  int ierr=MPI_Allreduce(MPI_IN_PLACE,&leaders_1hot[0],WorldSize,MPI_INT,MPI_SUM,communicator_world);
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
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  // allocate the shared window for our group
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  ShmCommBuf = 0;
  ierr = MPI_Win_allocate_shared(MAX_MPI_SHM_BYTES,1,MPI_INFO_NULL,ShmComm,&ShmCommBuf,&ShmWindow);
  assert(ierr==0);
  // KNL hack -- force to numa-domain 1 in flat
#if 0
  //#include <numaif.h>
  for(uint64_t page=0;page<MAX_MPI_SHM_BYTES;page+=4096){
    void *pages = (void *) ( page + ShmCommBuf );
    int status;
    int flags=MPOL_MF_MOVE_ALL;
    int nodes=1; // numa domain == MCDRAM
    unsigned long count=1;
    ierr= move_pages(0,count, &pages,&nodes,&status,flags);
    if (ierr && (page==0)) perror("numa relocate command failed");
  }
#endif
  MPI_Win_lock_all (MPI_MODE_NOCHECK, ShmWindow);
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Plan: allocate a fixed SHM region. Scratch that is just used via some scheme during stencil comms, with no allocate free.
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ShmCommBufs.resize(ShmSize);
  for(int r=0;r<ShmSize;r++){
    MPI_Aint sz;
    int dsp_unit;
    MPI_Win_shared_query (ShmWindow, r, &sz, &dsp_unit, &ShmCommBufs[r]);
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Verbose for now
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  if (WorldRank == 0){
    std::cout<<GridLogMessage<< "Grid MPI-3 configuration: detected ";
    std::cout<< WorldSize << " Ranks " ;
    std::cout<< GroupSize << " Nodes " ;
    std::cout<<  ShmSize  << " with ranks-per-node "<<std::endl;
    
    std::cout<<GridLogMessage     <<"Grid MPI-3 configuration: allocated shared memory region of size ";
    std::cout<<std::hex << MAX_MPI_SHM_BYTES <<" ShmCommBuf address = "<<ShmCommBuf << std::dec<<std::endl;

    for(int g=0;g<GroupSize;g++){
      std::cout<<GridLogMessage<<" Node "<<g<<" led by MPI rank "<<leaders_group[g]<<std::endl;
    }

    std::cout<<GridLogMessage<<" Boss Node Shm Pointers are {";
    for(int g=0;g<ShmSize;g++){
      std::cout<<std::hex<<ShmCommBufs[g]<<std::dec;
      if(g!=ShmSize-1) std::cout<<",";
      else std::cout<<"}"<<std::endl;
    }

  }
  
  for(int g=0;g<GroupSize;g++){
    if ( (ShmRank == 0) && (GroupRank==g) )  std::cout<<GridLogMessage<<"["<<g<<"] Node Group "<<g<<" is ranks {";
    for(int r=0;r<ShmSize;r++){
      if ( (ShmRank == 0) && (GroupRank==g) ) {
	std::cout<<MyGroup[r];
	if(r<ShmSize-1) std::cout<<",";
	else std::cout<<"}"<<std::endl;
      }
      MPI_Barrier(communicator_world);
    }
  }
  
  assert(ShmSetup==0);  ShmSetup=1;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Want to implement some magic ... Group sub-cubes into those on same node
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
  int ierr;

  communicator=communicator_world;

  _ndimension = processors.size();
  
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
  
  std::vector<int> WorldDims = processors;

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

  ////////////////////////////////////////////////////////////////
  // Establish torus of processes and nodes with sub-blockings
  ////////////////////////////////////////////////////////////////
  for(int d=0;d<_ndimension;d++){
    GroupDims[d] = WorldDims[d]/ShmDims[d];
  }

  ////////////////////////////////////////////////////////////////
  // Check processor counts match
  ////////////////////////////////////////////////////////////////
  _Nprocessors=1;
  _processors = processors;
  _processor_coor.resize(_ndimension);
  for(int i=0;i<_ndimension;i++){
    _Nprocessors*=_processors[i];
  }
  assert(WorldSize==_Nprocessors);
      
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
#if 0
  this->StencilBarrier();

  MPI_Request xrq;
  MPI_Request rrq;
  
  static int sequence;

  int ierr;
  int tag;
  int check;

  assert(dest != _processor);
  assert(from != _processor);
  
  int gdest = GroupRanks[dest];
  int gfrom = GroupRanks[from];
  int gme   = GroupRanks[_processor];

  sequence++;
  
  char *from_ptr = (char *)ShmCommBufs[ShmRank];

  int small = (bytes<MAX_MPI_SHM_BYTES);

  typedef uint64_t T;
  int words = bytes/sizeof(T);

  assert(((size_t)bytes &(sizeof(T)-1))==0);
  assert(gme == ShmRank);

  if ( small && (gdest !=MPI_UNDEFINED) ) {

    char *to_ptr   = (char *)ShmCommBufs[gdest];

    assert(gme != gdest);

    T *ip = (T *)xmit;
    T *op = (T *)to_ptr;
PARALLEL_FOR_LOOP 
    for(int w=0;w<words;w++) {
      op[w]=ip[w];
    }

    bcopy(&_processor,&to_ptr[bytes],sizeof(_processor));
    bcopy(&  sequence,&to_ptr[bytes+4],sizeof(sequence));
  } else { 
    ierr =MPI_Isend(xmit, bytes, MPI_CHAR,dest,_processor,communicator,&xrq);
    assert(ierr==0);
    list.push_back(xrq);
  }

  this->StencilBarrier();
  
  if (small && (gfrom !=MPI_UNDEFINED) ) {
    T *ip = (T *)from_ptr;
    T *op = (T *)recv;
PARALLEL_FOR_LOOP 
    for(int w=0;w<words;w++) {
      op[w]=ip[w];
    }
    bcopy(&from_ptr[bytes]  ,&tag  ,sizeof(tag));
    bcopy(&from_ptr[bytes+4],&check,sizeof(check));
    assert(check==sequence);
    assert(tag==from);
  } else { 
    ierr=MPI_Irecv(recv, bytes, MPI_CHAR,from,from,communicator,&rrq);
    assert(ierr==0);
    list.push_back(rrq);
  }

  this->StencilBarrier();

#else
  MPI_Request xrq;
  MPI_Request rrq;
  int rank = _processor;
  int ierr;
  ierr =MPI_Isend(xmit, bytes, MPI_CHAR,dest,_processor,communicator,&xrq);
  ierr|=MPI_Irecv(recv, bytes, MPI_CHAR,from,from,communicator,&rrq);
  
  assert(ierr==0);

  list.push_back(xrq);
  list.push_back(rrq);
#endif
}

void CartesianCommunicator::StencilSendToRecvFromBegin(std::vector<CommsRequest_t> &list,
						       void *xmit,
						       int dest,
						       void *recv,
						       int from,
						       int bytes)
{
  MPI_Request xrq;
  MPI_Request rrq;

  int ierr;

  assert(dest != _processor);
  assert(from != _processor);
  
  int gdest = GroupRanks[dest];
  int gfrom = GroupRanks[from];
  int gme   = GroupRanks[_processor];

  assert(gme == ShmRank);

  if ( gdest == MPI_UNDEFINED ) {
    ierr =MPI_Isend(xmit, bytes, MPI_CHAR,dest,_processor,communicator,&xrq);
    assert(ierr==0);
    list.push_back(xrq);
  }
  
  if ( gfrom ==MPI_UNDEFINED) {
    ierr=MPI_Irecv(recv, bytes, MPI_CHAR,from,from,communicator,&rrq);
    assert(ierr==0);
    list.push_back(rrq);
  }

}


void CartesianCommunicator::StencilSendToRecvFromComplete(std::vector<CommsRequest_t> &list)
{
  SendToRecvFromComplete(list);
}

void CartesianCommunicator::StencilBarrier(void)
{
  MPI_Win_sync (ShmWindow);   
  MPI_Barrier  (ShmComm);
  MPI_Win_sync (ShmWindow);   
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
		      communicator_world);
  assert(ierr==0);
}

}

