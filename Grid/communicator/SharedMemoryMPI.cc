/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/communicator/SharedMemory.cc

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: Christoph Lehner <christoph@lhnr.de>

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
#include <pwd.h>

#ifdef GRID_CUDA
#include <cuda_runtime_api.h>
#endif
#ifdef GRID_HIP
#include <hip/hip_runtime_api.h>
#endif
#ifdef GRID_SYCl

#endif

NAMESPACE_BEGIN(Grid); 
#define header "SharedMemoryMpi: "
/*Construct from an MPI communicator*/
void GlobalSharedMemory::Init(Grid_MPI_Comm comm)
{
  assert(_ShmSetup==0);
  WorldComm = comm;
  MPI_Comm_rank(WorldComm,&WorldRank);
  MPI_Comm_size(WorldComm,&WorldSize);
  // WorldComm, WorldSize, WorldRank

  /////////////////////////////////////////////////////////////////////
  // Split into groups that can share memory
  /////////////////////////////////////////////////////////////////////
#ifndef GRID_MPI3_SHM_NONE
  MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL,&WorldShmComm);
#else
  MPI_Comm_split(comm, WorldRank, 0, &WorldShmComm);
#endif

  MPI_Comm_rank(WorldShmComm     ,&WorldShmRank);
  MPI_Comm_size(WorldShmComm     ,&WorldShmSize);

  if ( WorldRank == 0) {
    std::cout << header " World communicator of size " <<WorldSize << std::endl;  
    std::cout << header " Node  communicator of size " <<WorldShmSize << std::endl;
  }
  // WorldShmComm, WorldShmSize, WorldShmRank

  // WorldNodes
  WorldNodes = WorldSize/WorldShmSize;
  assert( (WorldNodes * WorldShmSize) == WorldSize );

  // FIXME: Check all WorldShmSize are the same ?

  /////////////////////////////////////////////////////////////////////
  // find world ranks in our SHM group (i.e. which ranks are on our node)
  /////////////////////////////////////////////////////////////////////
  MPI_Group WorldGroup, ShmGroup;
  MPI_Comm_group (WorldComm, &WorldGroup); 
  MPI_Comm_group (WorldShmComm, &ShmGroup);

  std::vector<int> world_ranks(WorldSize);   for(int r=0;r<WorldSize;r++) world_ranks[r]=r;

  WorldShmRanks.resize(WorldSize); 
  MPI_Group_translate_ranks (WorldGroup,WorldSize,&world_ranks[0],ShmGroup, &WorldShmRanks[0]); 

  ///////////////////////////////////////////////////////////////////
  // Identify who is in my group and nominate the leader
  ///////////////////////////////////////////////////////////////////
  int g=0;
  std::vector<int> MyGroup;
  MyGroup.resize(WorldShmSize);
  for(int rank=0;rank<WorldSize;rank++){
    if(WorldShmRanks[rank]!=MPI_UNDEFINED){
      assert(g<WorldShmSize);
      MyGroup[g++] = rank;
    }
  }
  
  std::sort(MyGroup.begin(),MyGroup.end(),std::less<int>());
  int myleader = MyGroup[0];
  
  std::vector<int> leaders_1hot(WorldSize,0);
  std::vector<int> leaders_group(WorldNodes,0);
  leaders_1hot [ myleader ] = 1;
    
  ///////////////////////////////////////////////////////////////////
  // global sum leaders over comm world
  ///////////////////////////////////////////////////////////////////
  int ierr=MPI_Allreduce(MPI_IN_PLACE,&leaders_1hot[0],WorldSize,MPI_INT,MPI_SUM,WorldComm);
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
  // Identify the node of the group in which I (and my leader) live
  ///////////////////////////////////////////////////////////////////
  WorldNode=-1;
  for(int g=0;g<WorldNodes;g++){
    if (myleader == leaders_group[g]){
      WorldNode=g;
    }
  }
  assert(WorldNode!=-1);
  _ShmSetup=1;
}
// Gray encode support 
int BinaryToGray (int  binary) {
  int gray = (binary>>1)^binary;
  return gray;
}
int Log2Size(int TwoToPower,int MAXLOG2)
{
  int log2size = -1;
  for(int i=0;i<=MAXLOG2;i++){
    if ( (0x1<<i) == TwoToPower ) {
      log2size = i;
      break;
    }
  }
  return log2size;
}
void GlobalSharedMemory::OptimalCommunicator(const Coordinate &processors,Grid_MPI_Comm & optimal_comm)
{
  //////////////////////////////////////////////////////////////////////////////
  // Look and see if it looks like an HPE 8600 based on hostname conventions
  //////////////////////////////////////////////////////////////////////////////
  const int namelen = _POSIX_HOST_NAME_MAX;
  char name[namelen];
  int R;
  int I;
  int N;
  gethostname(name,namelen);
  int nscan = sscanf(name,"r%di%dn%d",&R,&I,&N) ;

  if(nscan==3 && HPEhypercube ) OptimalCommunicatorHypercube(processors,optimal_comm);
  else                          OptimalCommunicatorSharedMemory(processors,optimal_comm);
}
static inline int divides(int a,int b)
{
  return ( b == ( (b/a)*a ) );
}
void GlobalSharedMemory::GetShmDims(const Coordinate &WorldDims,Coordinate &ShmDims)
{
  ////////////////////////////////////////////////////////////////
  // Allow user to configure through environment variable
  ////////////////////////////////////////////////////////////////
  char* str = getenv(("GRID_SHM_DIMS_" + std::to_string(ShmDims.size())).c_str());
  if ( str ) {
    std::vector<int> IntShmDims;
    GridCmdOptionIntVector(std::string(str),IntShmDims);
    assert(IntShmDims.size() == WorldDims.size());
    long ShmSize = 1;
    for (int dim=0;dim<WorldDims.size();dim++) {
      ShmSize *= (ShmDims[dim] = IntShmDims[dim]);
      assert(divides(ShmDims[dim],WorldDims[dim]));
    }
    assert(ShmSize == WorldShmSize);
    return;
  }
  
  ////////////////////////////////////////////////////////////////
  // Powers of 2,3,5 only in prime decomposition for now
  ////////////////////////////////////////////////////////////////
  int ndimension = WorldDims.size();
  ShmDims=Coordinate(ndimension,1);

  std::vector<int> primes({2,3,5});

  int dim = 0;
  int last_dim = ndimension - 1;
  int AutoShmSize = 1;
  while(AutoShmSize != WorldShmSize) {
    int p;
    for(p=0;p<primes.size();p++) {
      int prime=primes[p];
      if ( divides(prime,WorldDims[dim]/ShmDims[dim])
        && divides(prime,WorldShmSize/AutoShmSize)  ) {
	AutoShmSize*=prime;
	ShmDims[dim]*=prime;
	last_dim = dim;
	break;
      }
    }
    if (p == primes.size() && last_dim == dim) {
      std::cerr << "GlobalSharedMemory::GetShmDims failed" << std::endl;
      exit(EXIT_FAILURE);
    }
    dim=(dim+1) %ndimension;
  }
}
void GlobalSharedMemory::OptimalCommunicatorHypercube(const Coordinate &processors,Grid_MPI_Comm & optimal_comm)
{
  ////////////////////////////////////////////////////////////////
  // Assert power of two shm_size.
  ////////////////////////////////////////////////////////////////
  int log2size = Log2Size(WorldShmSize,MAXLOG2RANKSPERNODE);
  assert(log2size != -1);

  ////////////////////////////////////////////////////////////////
  // Identify the hypercube coordinate of this node using hostname
  ////////////////////////////////////////////////////////////////
  // n runs 0...7 9...16 18...25 27...34     (8*4)  5 bits
  // i runs 0..7                                    3 bits
  // r runs 0..3                                    2 bits
  // 2^10 = 1024 nodes
  const int maxhdim = 10; 
  std::vector<int> HyperCubeCoords(maxhdim,0);
  std::vector<int> RootHyperCubeCoords(maxhdim,0);
  int R;
  int I;
  int N;
  const int namelen = _POSIX_HOST_NAME_MAX;
  char name[namelen];

  // Parse ICE-XA hostname to get hypercube location
  gethostname(name,namelen);
  int nscan = sscanf(name,"r%di%dn%d",&R,&I,&N) ;
  assert(nscan==3);

  int nlo = N%9;
  int nhi = N/9;
  uint32_t hypercoor = (R<<8)|(I<<5)|(nhi<<3)|nlo ;
  uint32_t rootcoor  = hypercoor;

  //////////////////////////////////////////////////////////////////
  // Print debug info
  //////////////////////////////////////////////////////////////////
  for(int d=0;d<maxhdim;d++){
    HyperCubeCoords[d] = (hypercoor>>d)&0x1;
  }

  std::string hname(name);
  //  std::cout << "hostname "<<hname<<std::endl;
  //  std::cout << "R " << R << " I " << I << " N "<< N
  //            << " hypercoor 0x"<<std::hex<<hypercoor<<std::dec<<std::endl;

  //////////////////////////////////////////////////////////////////
  // broadcast node 0's base coordinate for this partition.
  //////////////////////////////////////////////////////////////////
  MPI_Bcast(&rootcoor, sizeof(rootcoor), MPI_BYTE, 0, WorldComm); 
  hypercoor=hypercoor-rootcoor;
  assert(hypercoor<WorldSize);
  assert(hypercoor>=0);

  //////////////////////////////////////
  // Printing
  //////////////////////////////////////
  for(int d=0;d<maxhdim;d++){
    HyperCubeCoords[d] = (hypercoor>>d)&0x1;
  }

  ////////////////////////////////////////////////////////////////
  // Identify subblock of ranks on node spreading across dims
  // in a maximally symmetrical way
  ////////////////////////////////////////////////////////////////
  int ndimension              = processors.size();
  Coordinate processor_coor(ndimension);
  Coordinate WorldDims = processors;
  Coordinate ShmDims  (ndimension);  Coordinate NodeDims (ndimension);
  Coordinate ShmCoor  (ndimension);    Coordinate NodeCoor (ndimension);    Coordinate WorldCoor(ndimension);
  Coordinate HyperCoor(ndimension);

  GetShmDims(WorldDims,ShmDims);

  ////////////////////////////////////////////////////////////////
  // Establish torus of processes and nodes with sub-blockings
  ////////////////////////////////////////////////////////////////
  for(int d=0;d<ndimension;d++){
    NodeDims[d] = WorldDims[d]/ShmDims[d];
  }
  ////////////////////////////////////////////////////////////////
  // Map Hcube according to physical lattice 
  // must partition. Loop over dims and find out who would join.
  ////////////////////////////////////////////////////////////////
  int hcoor = hypercoor;
  for(int d=0;d<ndimension;d++){
     int bits = Log2Size(NodeDims[d],MAXLOG2RANKSPERNODE);
     int msk  = (0x1<<bits)-1;
     HyperCoor[d]=hcoor & msk;  
     HyperCoor[d]=BinaryToGray(HyperCoor[d]); // Space filling curve magic
     hcoor = hcoor >> bits;
  }
  ////////////////////////////////////////////////////////////////
  // Check processor counts match
  ////////////////////////////////////////////////////////////////
  int Nprocessors=1;
  for(int i=0;i<ndimension;i++){
    Nprocessors*=processors[i];
  }
  assert(WorldSize==Nprocessors);

  ////////////////////////////////////////////////////////////////
  // Establish mapping between lexico physics coord and WorldRank
  ////////////////////////////////////////////////////////////////
  int rank;

  Lexicographic::CoorFromIndexReversed(NodeCoor,WorldNode   ,NodeDims);

  for(int d=0;d<ndimension;d++) NodeCoor[d]=HyperCoor[d];

  Lexicographic::CoorFromIndexReversed(ShmCoor ,WorldShmRank,ShmDims);
  for(int d=0;d<ndimension;d++) WorldCoor[d] = NodeCoor[d]*ShmDims[d]+ShmCoor[d];
  Lexicographic::IndexFromCoorReversed(WorldCoor,rank,WorldDims);

  /////////////////////////////////////////////////////////////////
  // Build the new communicator
  /////////////////////////////////////////////////////////////////
  int ierr= MPI_Comm_split(WorldComm,0,rank,&optimal_comm);
  assert(ierr==0);
}
void GlobalSharedMemory::OptimalCommunicatorSharedMemory(const Coordinate &processors,Grid_MPI_Comm & optimal_comm)
{
  ////////////////////////////////////////////////////////////////
  // Identify subblock of ranks on node spreading across dims
  // in a maximally symmetrical way
  ////////////////////////////////////////////////////////////////
  int ndimension              = processors.size();
  Coordinate processor_coor(ndimension);
  Coordinate WorldDims = processors; Coordinate ShmDims(ndimension);  Coordinate NodeDims (ndimension);
  Coordinate ShmCoor(ndimension);    Coordinate NodeCoor(ndimension);   Coordinate WorldCoor(ndimension);

  GetShmDims(WorldDims,ShmDims);
  ////////////////////////////////////////////////////////////////
  // Establish torus of processes and nodes with sub-blockings
  ////////////////////////////////////////////////////////////////
  for(int d=0;d<ndimension;d++){
    NodeDims[d] = WorldDims[d]/ShmDims[d];
  }

  ////////////////////////////////////////////////////////////////
  // Check processor counts match
  ////////////////////////////////////////////////////////////////
  int Nprocessors=1;
  for(int i=0;i<ndimension;i++){
    Nprocessors*=processors[i];
  }
  assert(WorldSize==Nprocessors);

  ////////////////////////////////////////////////////////////////
  // Establish mapping between lexico physics coord and WorldRank
  ////////////////////////////////////////////////////////////////
  int rank;

  Lexicographic::CoorFromIndexReversed(NodeCoor,WorldNode   ,NodeDims);
  Lexicographic::CoorFromIndexReversed(ShmCoor ,WorldShmRank,ShmDims);
  for(int d=0;d<ndimension;d++) WorldCoor[d] = NodeCoor[d]*ShmDims[d]+ShmCoor[d];
  Lexicographic::IndexFromCoorReversed(WorldCoor,rank,WorldDims);

  /////////////////////////////////////////////////////////////////
  // Build the new communicator
  /////////////////////////////////////////////////////////////////
  int ierr= MPI_Comm_split(WorldComm,0,rank,&optimal_comm);
  assert(ierr==0);
}
////////////////////////////////////////////////////////////////////////////////////////////
// SHMGET
////////////////////////////////////////////////////////////////////////////////////////////
#ifdef GRID_MPI3_SHMGET
void GlobalSharedMemory::SharedMemoryAllocate(uint64_t bytes, int flags)
{
  std::cout << header "SharedMemoryAllocate "<< bytes<< " shmget implementation "<<std::endl;
  assert(_ShmSetup==1);
  assert(_ShmAlloc==0);

  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  // allocate the shared windows for our group
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  MPI_Barrier(WorldShmComm);
  WorldShmCommBufs.resize(WorldShmSize);
  std::vector<int> shmids(WorldShmSize);

  if ( WorldShmRank == 0 ) {
    for(int r=0;r<WorldShmSize;r++){
      size_t size = bytes;
      key_t key   = IPC_PRIVATE;
      int flags = IPC_CREAT | SHM_R | SHM_W;
#ifdef SHM_HUGETLB
      if (Hugepages) flags|=SHM_HUGETLB;
#endif
      if ((shmids[r]= shmget(key,size, flags)) ==-1) {
        int errsv = errno;
        printf("Errno %d\n",errsv);
        printf("key   %d\n",key);
        printf("size  %ld\n",size);
        printf("flags %d\n",flags);
        perror("shmget");
        exit(1);
      }
    }
  }
  MPI_Barrier(WorldShmComm);
  MPI_Bcast(&shmids[0],WorldShmSize*sizeof(int),MPI_BYTE,0,WorldShmComm);
  MPI_Barrier(WorldShmComm);

  for(int r=0;r<WorldShmSize;r++){
    WorldShmCommBufs[r] = (uint64_t *)shmat(shmids[r], NULL,0);
    if (WorldShmCommBufs[r] == (uint64_t *)-1) {
      perror("Shared memory attach failure");
      shmctl(shmids[r], IPC_RMID, NULL);
      exit(2);
    }
  }
  MPI_Barrier(WorldShmComm);
  ///////////////////////////////////
  // Mark for clean up
  ///////////////////////////////////
  for(int r=0;r<WorldShmSize;r++){
    shmctl(shmids[r], IPC_RMID,(struct shmid_ds *)NULL);
  }
  MPI_Barrier(WorldShmComm);

  _ShmAlloc=1;
  _ShmAllocBytes  = bytes;
}
#endif

////////////////////////////////////////////////////////////////////////////////////////////
// Hugetlbfs mapping intended
////////////////////////////////////////////////////////////////////////////////////////////
#if defined(GRID_CUDA) ||defined(GRID_HIP)  || defined(GRID_SYCL)

#if defined(GRID_SYCL) 
void GlobalSharedMemory::SharedMemoryAllocate(uint64_t bytes, int flags)
{
  void * ShmCommBuf ; 
  assert(_ShmSetup==1);
  assert(_ShmAlloc==0);

  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  // allocate the pointer array for shared windows for our group
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  MPI_Barrier(WorldShmComm);
  WorldShmCommBufs.resize(WorldShmSize);

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Each MPI rank should allocate our own buffer
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////
  ShmCommBuf = acceleratorAllocDevice(bytes);

  if (ShmCommBuf == (void *)NULL ) {
    std::cerr << " SharedMemoryMPI.cc acceleratorAllocDevice failed NULL pointer for " << bytes<<" bytes " << std::endl;
    exit(EXIT_FAILURE);  
  }
  //  if ( WorldRank == 0 ){
  if ( 1 ){
    std::cout << WorldRank << header " SharedMemoryMPI.cc acceleratorAllocDevice "<< bytes 
	      << "bytes at "<< std::hex<< ShmCommBuf <<std::dec<<" for comms buffers " <<std::endl;
  }
  SharedMemoryZero(ShmCommBuf,bytes);

  for(int r=0;r<WorldShmSize;r++){
    WorldShmCommBufs[r] = ShmCommBuf;
  }
  _ShmAllocBytes=bytes;
  _ShmAlloc=1;
}
#endif

#if defined(GRID_CUDA) ||defined(GRID_HIP) 
void GlobalSharedMemory::SharedMemoryAllocate(uint64_t bytes, int flags)
{
  void * ShmCommBuf ; 
  assert(_ShmSetup==1);
  assert(_ShmAlloc==0);

  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  // allocate the pointer array for shared windows for our group
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  MPI_Barrier(WorldShmComm);
  WorldShmCommBufs.resize(WorldShmSize);

  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  // TODO/FIXME : NOT ALL NVLINK BOARDS have full Peer to peer connectivity.
  // The annoyance is that they have partial peer 2 peer. This occurs on the 8 GPU blades.
  // e.g. DGX1, supermicro board, 
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  //  cudaDeviceGetP2PAttribute(&perfRank, cudaDevP2PAttrPerformanceRank, device1, device2);

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Each MPI rank should allocate our own buffer
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////
  ShmCommBuf = acceleratorAllocDevice(bytes);

  if (ShmCommBuf == (void *)NULL ) {
    std::cerr << " SharedMemoryMPI.cc acceleratorAllocDevice failed NULL pointer for " << bytes<<" bytes " << std::endl;
    exit(EXIT_FAILURE);  
  }
  //  if ( WorldRank == 0 ){
  if ( 1 ){
    std::cout << WorldRank << header " SharedMemoryMPI.cc acceleratorAllocDevice "<< bytes 
	      << "bytes at "<< std::hex<< ShmCommBuf <<std::dec<<" for comms buffers " <<std::endl;
  }
  SharedMemoryZero(ShmCommBuf,bytes);

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Loop over ranks/gpu's on our node
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////
  for(int r=0;r<WorldShmSize;r++){

#ifndef GRID_MPI3_SHM_NONE
    //////////////////////////////////////////////////
    // If it is me, pass around the IPC access key
    //////////////////////////////////////////////////
#ifdef GRID_CUDA
    cudaIpcMemHandle_t handle;
    if ( r==WorldShmRank ) { 
      auto err = cudaIpcGetMemHandle(&handle,ShmCommBuf);
      if ( err !=  cudaSuccess) {
	std::cerr << " SharedMemoryMPI.cc cudaIpcGetMemHandle failed for rank" << r <<" "<<cudaGetErrorString(err)<< std::endl;
	exit(EXIT_FAILURE);
      }
    }
#endif
#ifdef GRID_HIP
    hipIpcMemHandle_t handle;    
    if ( r==WorldShmRank ) { 
      auto err = hipIpcGetMemHandle(&handle,ShmCommBuf);
      if ( err !=  hipSuccess) {
	std::cerr << " SharedMemoryMPI.cc hipIpcGetMemHandle failed for rank" << r <<" "<<hipGetErrorString(err)<< std::endl;
	exit(EXIT_FAILURE);
      }
    }
#endif
    //////////////////////////////////////////////////
    // Share this IPC handle across the Shm Comm
    //////////////////////////////////////////////////
    { 
      int ierr=MPI_Bcast(&handle,
			 sizeof(handle),
			 MPI_BYTE,
			 r,
			 WorldShmComm);
      assert(ierr==0);
    }
    
    ///////////////////////////////////////////////////////////////
    // If I am not the source, overwrite thisBuf with remote buffer
    ///////////////////////////////////////////////////////////////
    void * thisBuf = ShmCommBuf;
#ifdef GRID_CUDA
    if ( r!=WorldShmRank ) { 
      auto err = cudaIpcOpenMemHandle(&thisBuf,handle,cudaIpcMemLazyEnablePeerAccess);
      if ( err !=  cudaSuccess) {
	std::cerr << " SharedMemoryMPI.cc cudaIpcOpenMemHandle failed for rank" << r <<" "<<cudaGetErrorString(err)<< std::endl;
	exit(EXIT_FAILURE);
      }
    }
#endif
#ifdef GRID_HIP
    if ( r!=WorldShmRank ) { 
      auto err = hipIpcOpenMemHandle(&thisBuf,handle,hipIpcMemLazyEnablePeerAccess);
      if ( err !=  hipSuccess) {
	std::cerr << " SharedMemoryMPI.cc hipIpcOpenMemHandle failed for rank" << r <<" "<<hipGetErrorString(err)<< std::endl;
	exit(EXIT_FAILURE);
      }
    }
#endif
    ///////////////////////////////////////////////////////////////
    // Save a copy of the device buffers
    ///////////////////////////////////////////////////////////////
    WorldShmCommBufs[r] = thisBuf;
#else
    WorldShmCommBufs[r] = ShmCommBuf;
#endif
  }

  _ShmAllocBytes=bytes;
  _ShmAlloc=1;
}
#endif

#else 
#ifdef GRID_MPI3_SHMMMAP
void GlobalSharedMemory::SharedMemoryAllocate(uint64_t bytes, int flags)
{
  std::cout << header "SharedMemoryAllocate "<< bytes<< " MMAP implementation "<< GRID_SHM_PATH <<std::endl;
  assert(_ShmSetup==1);
  assert(_ShmAlloc==0);
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  // allocate the shared windows for our group
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  MPI_Barrier(WorldShmComm);
  WorldShmCommBufs.resize(WorldShmSize);
  
  ////////////////////////////////////////////////////////////////////////////////////////////
  // Hugetlbfs and others map filesystems as mappable huge pages
  ////////////////////////////////////////////////////////////////////////////////////////////
  char shm_name [NAME_MAX];
  for(int r=0;r<WorldShmSize;r++){
    
    sprintf(shm_name,GRID_SHM_PATH "/Grid_mpi3_shm_%d_%d",WorldNode,r);
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
    if ( flags ) mmap_flag |= MAP_HUGETLB;
#endif
    void *ptr = (void *) mmap(NULL, bytes, PROT_READ | PROT_WRITE, mmap_flag,fd, 0); 
    if ( ptr == (void *)MAP_FAILED ) {    
      printf("mmap %s failed\n",shm_name);
      perror("failed mmap");      assert(0);    
    }
    assert(((uint64_t)ptr&0x3F)==0);
    close(fd);
    WorldShmCommBufs[r] =ptr;
    //    std::cout << header "Set WorldShmCommBufs["<<r<<"]="<<ptr<< "("<< bytes<< "bytes)"<<std::endl;
  }
  _ShmAlloc=1;
  _ShmAllocBytes  = bytes;
};
#endif // MMAP

#ifdef GRID_MPI3_SHM_NONE
void GlobalSharedMemory::SharedMemoryAllocate(uint64_t bytes, int flags)
{
  std::cout << header "SharedMemoryAllocate "<< bytes<< " MMAP anonymous implementation "<<std::endl;
  assert(_ShmSetup==1);
  assert(_ShmAlloc==0);
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  // allocate the shared windows for our group
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  MPI_Barrier(WorldShmComm);
  WorldShmCommBufs.resize(WorldShmSize);
  
  ////////////////////////////////////////////////////////////////////////////////////////////
  // Hugetlbf and others map filesystems as mappable huge pages
  ////////////////////////////////////////////////////////////////////////////////////////////
  char shm_name [NAME_MAX];
  assert(WorldShmSize == 1);
  for(int r=0;r<WorldShmSize;r++){
    
    int fd=-1;
    int mmap_flag = MAP_SHARED |MAP_ANONYMOUS ;
#ifdef MAP_POPULATE    
    mmap_flag|=MAP_POPULATE;
#endif
#ifdef MAP_HUGETLB
    if ( flags ) mmap_flag |= MAP_HUGETLB;
#endif
    void *ptr = (void *) mmap(NULL, bytes, PROT_READ | PROT_WRITE, mmap_flag,fd, 0); 
    if ( ptr == (void *)MAP_FAILED ) {    
      printf("mmap %s failed\n",shm_name);
      perror("failed mmap");      assert(0);    
    }
    assert(((uint64_t)ptr&0x3F)==0);
    close(fd);
    WorldShmCommBufs[r] =ptr;
    //    std::cout << "Set WorldShmCommBufs["<<r<<"]="<<ptr<< "("<< bytes<< "bytes)"<<std::endl;
  }
  _ShmAlloc=1;
  _ShmAllocBytes  = bytes;
};
#endif // MMAP

#ifdef GRID_MPI3_SHMOPEN
////////////////////////////////////////////////////////////////////////////////////////////
// POSIX SHMOPEN ; as far as I know Linux does not allow EXPLICIT HugePages with this case
// tmpfs (Larry Meadows says) does not support explicit huge page, and this is used for 
// the posix shm virtual file system
////////////////////////////////////////////////////////////////////////////////////////////
void GlobalSharedMemory::SharedMemoryAllocate(uint64_t bytes, int flags)
{ 
  std::cout << header "SharedMemoryAllocate "<< bytes<< " SHMOPEN implementation "<<std::endl;
  assert(_ShmSetup==1);
  assert(_ShmAlloc==0); 
  MPI_Barrier(WorldShmComm);
  WorldShmCommBufs.resize(WorldShmSize);

  char shm_name [NAME_MAX];
  if ( WorldShmRank == 0 ) {
    for(int r=0;r<WorldShmSize;r++){
	
      size_t size = bytes;
      
      struct passwd *pw = getpwuid (getuid());
      sprintf(shm_name,"/Grid_%s_mpi3_shm_%d_%d",pw->pw_name,WorldNode,r);
      
      shm_unlink(shm_name);
      int fd=shm_open(shm_name,O_RDWR|O_CREAT,0666);
      if ( fd < 0 ) {	perror("failed shm_open");	assert(0);      }
      ftruncate(fd, size);
	
      int mmap_flag = MAP_SHARED;
#ifdef MAP_POPULATE 
      mmap_flag |= MAP_POPULATE;
#endif
#ifdef MAP_HUGETLB
      if (flags) mmap_flag |= MAP_HUGETLB;
#endif
      void * ptr =  mmap(NULL,size, PROT_READ | PROT_WRITE, mmap_flag, fd, 0);
      
      if ( ptr == (void * )MAP_FAILED ) {       
	perror("failed mmap");     
	assert(0);    
      }
      assert(((uint64_t)ptr&0x3F)==0);
      
      WorldShmCommBufs[r] =ptr;
      close(fd);
    }
  }

  MPI_Barrier(WorldShmComm);
  
  if ( WorldShmRank != 0 ) { 
    for(int r=0;r<WorldShmSize;r++){

      size_t size = bytes ;
      
      struct passwd *pw = getpwuid (getuid());
      sprintf(shm_name,"/Grid_%s_mpi3_shm_%d_%d",pw->pw_name,WorldNode,r);
      
      int fd=shm_open(shm_name,O_RDWR,0666);
      if ( fd<0 ) {	perror("failed shm_open");	assert(0);      }
      
      void * ptr =  mmap(NULL,size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
      if ( ptr == MAP_FAILED ) {       perror("failed mmap");      assert(0);    }
      assert(((uint64_t)ptr&0x3F)==0);
      WorldShmCommBufs[r] =ptr;

      close(fd);
    }
  }
  _ShmAlloc=1;
  _ShmAllocBytes = bytes;
}
#endif
#endif // End NVCC case for GPU device buffers

/////////////////////////////////////////////////////////////////////////
// Routines accessing shared memory should route through for GPU safety
/////////////////////////////////////////////////////////////////////////
void GlobalSharedMemory::SharedMemoryZero(void *dest,size_t bytes)
{
#if defined(GRID_CUDA) || defined(GRID_HIP) || defined(GRID_SYCL)
  acceleratorMemSet(dest,0,bytes);
#else
  bzero(dest,bytes);
#endif
}
void GlobalSharedMemory::SharedMemoryCopy(void *dest,void *src,size_t bytes)
{
#if defined(GRID_CUDA) || defined(GRID_HIP) || defined(GRID_SYCL)
  acceleratorCopyToDevice(src,dest,bytes);
#else   
  bcopy(src,dest,bytes);
#endif
}
////////////////////////////////////////////////////////
// Global shared functionality finished
// Now move to per communicator functionality
////////////////////////////////////////////////////////
void SharedMemory::SetCommunicator(Grid_MPI_Comm comm)
{
  int rank, size;
  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&size);
  ShmRanks.resize(size);

  /////////////////////////////////////////////////////////////////////
  // Split into groups that can share memory
  /////////////////////////////////////////////////////////////////////
#ifndef GRID_MPI3_SHM_NONE
  MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL,&ShmComm);
#else
  MPI_Comm_split(comm, rank, 0, &ShmComm);
#endif
  MPI_Comm_rank(ShmComm     ,&ShmRank);
  MPI_Comm_size(ShmComm     ,&ShmSize);
  ShmCommBufs.resize(ShmSize);

  //////////////////////////////////////////////////////////////////////
  // Map ShmRank to WorldShmRank and use the right buffer
  //////////////////////////////////////////////////////////////////////
  assert (GlobalSharedMemory::ShmAlloc()==1);
  heap_size = GlobalSharedMemory::ShmAllocBytes();
  for(int r=0;r<ShmSize;r++){

    uint32_t wsr = (r==ShmRank) ? GlobalSharedMemory::WorldShmRank : 0 ;

    MPI_Allreduce(MPI_IN_PLACE,&wsr,1,MPI_UINT32_T,MPI_SUM,ShmComm);

    ShmCommBufs[r] = GlobalSharedMemory::WorldShmCommBufs[wsr];
  }
  ShmBufferFreeAll();

  /////////////////////////////////////////////////////////////////////
  // find comm ranks in our SHM group (i.e. which ranks are on our node)
  /////////////////////////////////////////////////////////////////////
  MPI_Group FullGroup, ShmGroup;
  MPI_Comm_group (comm   , &FullGroup); 
  MPI_Comm_group (ShmComm, &ShmGroup);

  std::vector<int> ranks(size);   for(int r=0;r<size;r++) ranks[r]=r;
  MPI_Group_translate_ranks (FullGroup,size,&ranks[0],ShmGroup, &ShmRanks[0]); 

#ifdef GRID_SHM_FORCE_MPI
  // Hide the shared memory path between ranks
  {
    for(int r=0;r<size;r++){
      if ( r!=rank ) {
	ShmRanks[r] = MPI_UNDEFINED;
      }
    }
  }
#endif

  SharedMemoryTest();
}
//////////////////////////////////////////////////////////////////
// On node barrier
//////////////////////////////////////////////////////////////////
void SharedMemory::ShmBarrier(void)
{
  MPI_Barrier  (ShmComm);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////
// Test the shared memory is working
//////////////////////////////////////////////////////////////////////////////////////////////////////////
void SharedMemory::SharedMemoryTest(void)
{
  ShmBarrier();
  uint64_t check[3];
  uint64_t magic = 0x5A5A5A;
  if ( ShmRank == 0 ) {
    for(uint64_t r=0;r<ShmSize;r++){
       check[0]=GlobalSharedMemory::WorldNode;
       check[1]=r;
       check[2]=magic;
       GlobalSharedMemory::SharedMemoryCopy( ShmCommBufs[r], check, 3*sizeof(uint64_t));
    }
  }
  ShmBarrier();
  for(uint64_t r=0;r<ShmSize;r++){
    ShmBarrier();
    GlobalSharedMemory::SharedMemoryCopy(check,ShmCommBufs[r], 3*sizeof(uint64_t));
    ShmBarrier();
    assert(check[0]==GlobalSharedMemory::WorldNode);
    assert(check[1]==r);
    assert(check[2]==magic);
    ShmBarrier();
  }
}

void *SharedMemory::ShmBuffer(int rank)
{
  int gpeer = ShmRanks[rank];
  if (gpeer == MPI_UNDEFINED){
    return NULL;
  } else { 
    return ShmCommBufs[gpeer];
  }
}
void *SharedMemory::ShmBufferTranslate(int rank,void * local_p)
{
  int gpeer = ShmRanks[rank];
  assert(gpeer!=ShmRank); // never send to self
  if (gpeer == MPI_UNDEFINED){
    return NULL;
  } else { 
    uint64_t offset = (uint64_t)local_p - (uint64_t)ShmCommBufs[ShmRank];
    uint64_t remote = (uint64_t)ShmCommBufs[gpeer]+offset;
    return (void *) remote;
  }
}
SharedMemory::~SharedMemory()
{
  int MPI_is_finalised;  MPI_Finalized(&MPI_is_finalised);
  if ( !MPI_is_finalised ) { 
    MPI_Comm_free(&ShmComm);
  }
};

NAMESPACE_END(Grid); 

