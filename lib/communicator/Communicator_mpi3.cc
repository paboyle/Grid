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

#include <mpi.h>

#include <semaphore.h>
#include <fcntl.h>
#include <unistd.h>
#include <limits.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <sys/mman.h>
#include <zlib.h>
#ifdef HAVE_NUMAIF_H
#include <numaif.h>
#endif


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

int CartesianCommunicator::NodeCount(void)    { return GroupSize;};
int CartesianCommunicator::RankCount(void)    { return WorldSize;};


#undef FORCE_COMMS
void *CartesianCommunicator::ShmBufferSelf(void)
{
  return ShmCommBufs[ShmRank];
}
void *CartesianCommunicator::ShmBuffer(int rank)
{
  int gpeer = GroupRanks[rank];
#ifdef FORCE_COMMS
  return NULL;
#endif
  if (gpeer == MPI_UNDEFINED){
    return NULL;
  } else { 
    return ShmCommBufs[gpeer];
  }
}
void *CartesianCommunicator::ShmBufferTranslate(int rank,void * local_p)
{
  static int count =0;
  int gpeer = GroupRanks[rank];
  assert(gpeer!=ShmRank); // never send to self
  assert(rank!=WorldRank);// never send to self
#ifdef FORCE_COMMS
  return NULL;
#endif
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
  int provided;
  //  mtrace();

  MPI_Initialized(&flag); // needed to coexist with other libs apparently
  if ( !flag ) {
    MPI_Init_thread(argc,argv,MPI_THREAD_MULTIPLE,&provided);
    assert (provided == MPI_THREAD_MULTIPLE);
  }

  Grid_quiesce_nodes();

  MPI_Comm_dup (MPI_COMM_WORLD,&communicator_world);
  MPI_Comm_rank(communicator_world,&WorldRank);
  MPI_Comm_size(communicator_world,&WorldSize);

  if ( WorldRank == 0 ) {
    std::cout << GridLogMessage<< "Initialising MPI "<< WorldRank <<"/"<<WorldSize <<std::endl;
  }

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
  MPI_Barrier(ShmComm);

  ShmCommBuf = 0;
  ShmCommBufs.resize(ShmSize);

  ////////////////////////////////////////////////////////////////////////////////////////////
  // Hugetlbf and others map filesystems as mappable huge pages
  ////////////////////////////////////////////////////////////////////////////////////////////
#ifdef GRID_MPI3_SHMMMAP
  char shm_name [NAME_MAX];
  for(int r=0;r<ShmSize;r++){
    
    size_t size = CartesianCommunicator::MAX_MPI_SHM_BYTES;
    sprintf(shm_name,GRID_SHM_PATH "/Grid_mpi3_shm_%d_%d",GroupRank,r);
    //sprintf(shm_name,"/var/lib/hugetlbfs/group/wheel/pagesize-2MB/" "Grid_mpi3_shm_%d_%d",GroupRank,r);
    //    printf("Opening file %s \n",shm_name);
    int fd=open(shm_name,O_RDWR|O_CREAT,0666);
    if ( fd == -1) { 
      printf("open %s failed\n",shm_name);
      perror("open hugetlbfs");
      exit(0);
    }
    int mmap_flag = MAP_SHARED ;
#ifdef MAP_POPULATE    
    mmap_flag|=MAP_POPULATE;
#endif
#ifdef MAP_HUGETLB
    if ( Hugepages ) mmap_flag |= MAP_HUGETLB;
#endif
    void *ptr = (void *) mmap(NULL, MAX_MPI_SHM_BYTES, PROT_READ | PROT_WRITE, mmap_flag,fd, 0); 
    if ( ptr == (void *)MAP_FAILED ) {    
      printf("mmap %s failed\n",shm_name);
      perror("failed mmap");      assert(0);    
    }
    assert(((uint64_t)ptr&0x3F)==0);
    ShmCommBufs[r] =ptr;
    
  }
#endif
  ////////////////////////////////////////////////////////////////////////////////////////////
  // POSIX SHMOPEN ; as far as I know Linux does not allow EXPLICIT HugePages with this case
  // tmpfs (Larry Meadows says) does not support explicit huge page, and this is used for 
  // the posix shm virtual file system
  ////////////////////////////////////////////////////////////////////////////////////////////
#ifdef GRID_MPI3_SHMOPEN
  char shm_name [NAME_MAX];
  if ( ShmRank == 0 ) {
    for(int r=0;r<ShmSize;r++){

      size_t size = CartesianCommunicator::MAX_MPI_SHM_BYTES;

      sprintf(shm_name,"/Grid_mpi3_shm_%d_%d",GroupRank,r);

      shm_unlink(shm_name);
      int fd=shm_open(shm_name,O_RDWR|O_CREAT,0666);
      if ( fd < 0 ) {	perror("failed shm_open");	assert(0);      }
      ftruncate(fd, size);
      
      int mmap_flag = MAP_SHARED;
#ifdef MAP_POPULATE 
      mmap_flag |= MAP_POPULATE;
#endif
#ifdef MAP_HUGETLB
      if (Hugepages) mmap_flag |= MAP_HUGETLB;
#endif
      void * ptr =  mmap(NULL,size, PROT_READ | PROT_WRITE, mmap_flag, fd, 0);

      if ( ptr == (void * )MAP_FAILED ) {       perror("failed mmap");      assert(0);    }
      assert(((uint64_t)ptr&0x3F)==0);

// Experiments; Experiments; Try to force numa domain on the shm segment if we have numaif.h
#if 0
//#ifdef HAVE_NUMAIF_H
	int status;
	int flags=MPOL_MF_MOVE;
#ifdef KNL
	int nodes=1; // numa domain == MCDRAM
	// Find out if in SNC2,SNC4 mode ?
#else
	int nodes=r; // numa domain == MPI ID
#endif
	unsigned long count=1;
	for(uint64_t page=0;page<size;page+=4096){
	  void *pages = (void *) ( page + (uint64_t)ptr );
	  uint64_t *cow_it = (uint64_t *)pages;	*cow_it = 1;
	  ierr= move_pages(0,count, &pages,&nodes,&status,flags);
	  if (ierr && (page==0)) perror("numa relocate command failed");
	}
#endif
	ShmCommBufs[r] =ptr;
      
    }
  }

  MPI_Barrier(ShmComm);

  if ( ShmRank != 0 ) { 
    for(int r=0;r<ShmSize;r++){
      size_t size = CartesianCommunicator::MAX_MPI_SHM_BYTES ;
    
      sprintf(shm_name,"/Grid_mpi3_shm_%d_%d",GroupRank,r);

      int fd=shm_open(shm_name,O_RDWR,0666);
      if ( fd<0 ) {	perror("failed shm_open");	assert(0);      }

      void * ptr =  mmap(NULL,size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
      if ( ptr == MAP_FAILED ) {       perror("failed mmap");      assert(0);    }
      assert(((uint64_t)ptr&0x3F)==0);
      ShmCommBufs[r] =ptr;
    }
  }
#endif
  ////////////////////////////////////////////////////////////////////////////////////////////
  // SHMGET SHMAT and SHM_HUGETLB flag
  ////////////////////////////////////////////////////////////////////////////////////////////
#ifdef GRID_MPI3_SHMGET
  std::vector<int> shmids(ShmSize);

  if ( ShmRank == 0 ) {
    for(int r=0;r<ShmSize;r++){
      size_t size = CartesianCommunicator::MAX_MPI_SHM_BYTES;
      key_t key   = IPC_PRIVATE;
      int flags = IPC_CREAT | SHM_R | SHM_W;
#ifdef SHM_HUGETLB
      if (Hugepages) flags|=SHM_HUGETLB;
#endif
      if ((shmids[r]= shmget(key,size, flags)) ==-1) {
	int errsv = errno;
	printf("Errno %d\n",errsv);
	printf("key   %d\n",key);
	printf("size  %lld\n",size);
	printf("flags %d\n",flags);
	perror("shmget");
	exit(1);
      } else { 
	printf("shmid: 0x%x\n", shmids[r]);
      }
    }
  }
  MPI_Barrier(ShmComm);
  MPI_Bcast(&shmids[0],ShmSize*sizeof(int),MPI_BYTE,0,ShmComm);
  MPI_Barrier(ShmComm);

  for(int r=0;r<ShmSize;r++){
    ShmCommBufs[r] = (uint64_t *)shmat(shmids[r], NULL,0);
    if (ShmCommBufs[r] == (uint64_t *)-1) {
      perror("Shared memory attach failure");
      shmctl(shmids[r], IPC_RMID, NULL);
      exit(2);
    }
    printf("shmaddr: %p\n", ShmCommBufs[r]);
  }
  MPI_Barrier(ShmComm);
  // Mark for clean up
  for(int r=0;r<ShmSize;r++){
    shmctl(shmids[r], IPC_RMID,(struct shmid_ds *)NULL);
  }
  MPI_Barrier(ShmComm);

#endif
  ShmCommBuf         = ShmCommBufs[ShmRank];

  MPI_Barrier(ShmComm);
  if ( ShmRank == 0 ) {
    for(int r=0;r<ShmSize;r++){
      uint64_t * check = (uint64_t *) ShmCommBufs[r];
      check[0] = GroupRank;
      check[1] = r;
      check[2] = 0x5A5A5A;
    }
  }

  MPI_Barrier(ShmComm);
  for(int r=0;r<ShmSize;r++){
    uint64_t * check = (uint64_t *) ShmCommBufs[r];
    
    assert(check[0]==GroupRank);
    assert(check[1]==r);
    assert(check[2]==0x5A5A5A);

  }
  MPI_Barrier(ShmComm);

  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Verbose for now
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  if (WorldRank == 0){
    std::cout<<GridLogMessage<< "Grid MPI-3 configuration: detected ";
    std::cout<< WorldSize << " Ranks " ;
    std::cout<< GroupSize << " Nodes " ;
    std::cout<< " with "<< ShmSize  << " ranks-per-node "<<std::endl;
    
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
	else std::cout<<"}"<<std::endl<<std::flush;
      }
      MPI_Barrier(communicator_world);
    }
  }

  assert(ShmSetup==0);  ShmSetup=1;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Want to implement some magic ... Group sub-cubes into those on same node
////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CartesianCommunicator::ShiftedRanks(int dim,int shift,int &dest,int &source)
{
  std::vector<int> coor = _processor_coor; // my coord
  assert(std::abs(shift) <_processors[dim]);

  coor[dim] = (_processor_coor[dim] + shift + _processors[dim])%_processors[dim];
  Lexicographic::IndexFromCoor(coor,source,_processors);
  source = LexicographicToWorldRank[source];

  coor[dim] = (_processor_coor[dim] - shift + _processors[dim])%_processors[dim];
  Lexicographic::IndexFromCoor(coor,dest,_processors);
  dest = LexicographicToWorldRank[dest];

}// rank is world rank.

int CartesianCommunicator::RankFromProcessorCoor(std::vector<int> &coor)
{
  int rank;
  Lexicographic::IndexFromCoor(coor,rank,_processors);
  rank = LexicographicToWorldRank[rank];
  return rank;
}// rank is world rank

void  CartesianCommunicator::ProcessorCoorFromRank(int rank, std::vector<int> &coor)
{
  int lr=-1;
  for(int r=0;r<WorldSize;r++){// map world Rank to lexico and then to coor
    if( LexicographicToWorldRank[r]==rank) lr = r;
  }
  assert(lr!=-1);
  Lexicographic::CoorFromIndex(coor,lr,_processors);
}

//////////////////////////////////
// Try to subdivide communicator
//////////////////////////////////
CartesianCommunicator::CartesianCommunicator(const std::vector<int> &processors,const CartesianCommunicator &parent) 
  : CartesianCommunicator(processors) 
{
  std::cout << "Attempts to split MPI3 communicators will fail until implemented" <<std::endl;
}
CartesianCommunicator::CartesianCommunicator(const std::vector<int> &processors)
{ 
  int ierr;
  communicator=communicator_world;

  _ndimension = processors.size();

  communicator_halo.resize (2*_ndimension);
  for(int i=0;i<_ndimension*2;i++){
    MPI_Comm_dup(communicator,&communicator_halo[i]);
  }

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
  std::vector<int> WorldDims = processors;

  ShmDims.resize  (_ndimension,1);
  GroupDims.resize(_ndimension);
  ShmCoor.resize  (_ndimension);
  GroupCoor.resize(_ndimension);
  WorldCoor.resize(_ndimension);

  int dim = 0;
  for(int l2=0;l2<log2size;l2++){
    while ( (WorldDims[dim] / ShmDims[dim]) <= 1 ) dim=(dim+1)%_ndimension;
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
  // Verbose
  ////////////////////////////////////////////////////////////////
#if 0
  std::cout<< GridLogMessage << "MPI-3 usage "<<std::endl;
  std::cout<< GridLogMessage << "SHM   ";
  for(int d=0;d<_ndimension;d++){
    std::cout<< ShmDims[d] <<" ";
  }
  std::cout<< std::endl;

  std::cout<< GridLogMessage << "Group ";
  for(int d=0;d<_ndimension;d++){
    std::cout<< GroupDims[d] <<" ";
  }
  std::cout<< std::endl;

  std::cout<< GridLogMessage<<"World ";
  for(int d=0;d<_ndimension;d++){
    std::cout<< WorldDims[d] <<" ";
  }
  std::cout<< std::endl;
#endif
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
  ////////////////////////////////////////////////////////////////
  Lexicographic::CoorFromIndex(GroupCoor,GroupRank,GroupDims);
  Lexicographic::CoorFromIndex(ShmCoor,ShmRank,ShmDims);
  for(int d=0;d<_ndimension;d++){
    WorldCoor[d] = GroupCoor[d]*ShmDims[d]+ShmCoor[d];
  }
  _processor_coor = WorldCoor;
  _processor      = WorldRank;

  ///////////////////////////////////////////////////////////////////
  // global sum Lexico to World mapping
  ///////////////////////////////////////////////////////////////////
  int lexico;
  LexicographicToWorldRank.resize(WorldSize,0);
  Lexicographic::IndexFromCoor(WorldCoor,lexico,WorldDims);
  LexicographicToWorldRank[lexico] = WorldRank;
  ierr=MPI_Allreduce(MPI_IN_PLACE,&LexicographicToWorldRank[0],WorldSize,MPI_INT,MPI_SUM,communicator);
  assert(ierr==0);

  for(int i=0;i<WorldSize;i++){

    int wr = LexicographicToWorldRank[i];
    //    int wr = i;

    std::vector<int> coor(_ndimension);
    ProcessorCoorFromRank(wr,coor); // from world rank
    int ck = RankFromProcessorCoor(coor);
    assert(ck==wr);

    if ( wr == WorldRank ) { 
      for(int j=0;j<coor.size();j++) {
	assert(coor[j] == _processor_coor[j]);
      }
    }
    /*
    std::cout << GridLogMessage<< " Lexicographic "<<i;
    std::cout << " MPI rank      "<<wr;
    std::cout << " Coor          ";
    for(int j=0;j<coor.size();j++) std::cout << coor[j];
    std::cout<< std::endl;
    */
    /////////////////////////////////////////////////////
    // Check everyone agrees on everyone elses coords
    /////////////////////////////////////////////////////
    std::vector<int> mcoor = coor;
    this->Broadcast(0,(void *)&mcoor[0],mcoor.size()*sizeof(int));
    for(int d = 0 ; d< _ndimension; d++) {
      assert(coor[d] == mcoor[d]);
    }
  }
};
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
// Basic Halo comms primitive
void CartesianCommunicator::SendToRecvFrom(void *xmit,
					   int dest,
					   void *recv,
					   int from,
					   int bytes)
{
  std::vector<CommsRequest_t> reqs(0);
  //    unsigned long  xcrc = crc32(0L, Z_NULL, 0);
  //    unsigned long  rcrc = crc32(0L, Z_NULL, 0);
  //    xcrc = crc32(xcrc,(unsigned char *)xmit,bytes);
  SendToRecvFromBegin(reqs,xmit,dest,recv,from,bytes);
  SendToRecvFromComplete(reqs);
  //    rcrc = crc32(rcrc,(unsigned char *)recv,bytes);
  //    printf("proc %d SendToRecvFrom %d bytes %lx %lx\n",_processor,bytes,xcrc,rcrc);
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

double CartesianCommunicator::StencilSendToRecvFrom( void *xmit,
						     int dest,
						     void *recv,
						     int from,
						     int bytes,int dir)
{
  std::vector<CommsRequest_t> list;
  double offbytes = StencilSendToRecvFromBegin(list,xmit,dest,recv,from,bytes,dir);
  StencilSendToRecvFromComplete(list,dir);
  return offbytes;
}

double CartesianCommunicator::StencilSendToRecvFromBegin(std::vector<CommsRequest_t> &list,
							 void *xmit,
							 int dest,
							 void *recv,
							 int from,
							 int bytes,int dir)
{
  int ncomm  =communicator_halo.size(); 
  int commdir=dir%ncomm;

  MPI_Request xrq;
  MPI_Request rrq;

  int ierr;
  int gdest = GroupRanks[dest];
  int gfrom = GroupRanks[from];
  int gme   = GroupRanks[_processor];

  assert(dest != _processor);
  assert(from != _processor);
  assert(gme  == ShmRank);
  double off_node_bytes=0.0;

#ifdef FORCE_COMMS
  gdest = MPI_UNDEFINED;
  gfrom = MPI_UNDEFINED;
#endif
  if ( gfrom ==MPI_UNDEFINED) {
    ierr=MPI_Irecv(recv, bytes, MPI_CHAR,from,from,communicator_halo[commdir],&rrq);
    assert(ierr==0);
    list.push_back(rrq);
    off_node_bytes+=bytes;
  }

  if ( gdest == MPI_UNDEFINED ) {
    ierr =MPI_Isend(xmit, bytes, MPI_CHAR,dest,_processor,communicator_halo[commdir],&xrq);
    assert(ierr==0);
    list.push_back(xrq);
    off_node_bytes+=bytes;
  }

  if ( CommunicatorPolicy == CommunicatorPolicySequential ) { 
    this->StencilSendToRecvFromComplete(list,dir);
  }

  return off_node_bytes;
}
void CartesianCommunicator::StencilSendToRecvFromComplete(std::vector<CommsRequest_t> &waitall,int dir)
{
  SendToRecvFromComplete(waitall);
}
void CartesianCommunicator::StencilBarrier(void)
{
  MPI_Barrier  (ShmComm);
}
void CartesianCommunicator::SendToRecvFromComplete(std::vector<CommsRequest_t> &list)
{
  int nreq=list.size();

  if (nreq==0) return;

  std::vector<MPI_Status> status(nreq);
  int ierr = MPI_Waitall(nreq,&list[0],&status[0]);
  assert(ierr==0);
  list.resize(0);
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

