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

#define Mheader "SharedMemoryMpi: "

#include <Grid/GridCore.h>
#include <pwd.h>

#ifdef GRID_FJMPI
#include <mpi-ext.h>
#endif

#ifdef GRID_CUDA
#include <cuda_runtime_api.h>
#endif
#ifdef GRID_HIP
#include <hip/hip_runtime_api.h>
#endif
#ifdef GRID_SYCL
#ifdef ACCELERATOR_AWARE_MPI
#define GRID_SYCL_LEVEL_ZERO_IPC
#define SHM_SOCKETS
#else
#endif 
#include <syscall.h>
#endif

#include <sys/socket.h>
#include <sys/un.h>

NAMESPACE_BEGIN(Grid); 

#ifdef SHM_SOCKETS

/*
 * Barbaric extra intranode communication route in case we need sockets to pass FDs
 * Forced by level_zero not being nicely designed
 */
static int sock;
static const char *sock_path_fmt = "/tmp/GridUnixSocket.%d";
static char sock_path[256];
class UnixSockets {
public:
  static void Open(int rank)
  {
    int errnum;

    sock = socket(AF_UNIX, SOCK_DGRAM, 0);  GRID_ASSERT(sock>0);

    struct sockaddr_un sa_un = { 0 };
    sa_un.sun_family = AF_UNIX;
    snprintf(sa_un.sun_path, sizeof(sa_un.sun_path),sock_path_fmt,rank);
    unlink(sa_un.sun_path);
    if (bind(sock, (struct sockaddr *)&sa_un, sizeof(sa_un))) {
      perror("bind failure");
      exit(EXIT_FAILURE);
    }
  }

  static int RecvFileDescriptor(void)
  {
    int n;
    int fd;
    char buf[1];
    struct iovec iov;
    struct msghdr msg;
    struct cmsghdr *cmsg;
    char cms[CMSG_SPACE(sizeof(int))];

    iov.iov_base = buf;
    iov.iov_len = 1;

    memset(&msg, 0, sizeof msg);
    msg.msg_name = 0;
    msg.msg_namelen = 0;
    msg.msg_iov = &iov;
    msg.msg_iovlen = 1;

    msg.msg_control = (caddr_t)cms;
    msg.msg_controllen = sizeof cms;

    if((n=recvmsg(sock, &msg, 0)) < 0) {
      perror("recvmsg failed");
      return -1;
    }
    if(n == 0){
      perror("recvmsg returned 0");
      return -1;
    }
    cmsg = CMSG_FIRSTHDR(&msg);

    memmove(&fd, CMSG_DATA(cmsg), sizeof(int));

    return fd;
  }

  static void SendFileDescriptor(int fildes,int xmit_to_rank)
  {
    struct msghdr msg;
    struct iovec iov;
    struct cmsghdr *cmsg = NULL;
    char ctrl[CMSG_SPACE(sizeof(int))];
    char data = ' ';

    memset(&msg, 0, sizeof(struct msghdr));
    memset(ctrl, 0, CMSG_SPACE(sizeof(int)));
    iov.iov_base = &data;
    iov.iov_len = sizeof(data);
    
    sprintf(sock_path,sock_path_fmt,xmit_to_rank);
    
    struct sockaddr_un sa_un = { 0 };
    sa_un.sun_family = AF_UNIX;
    snprintf(sa_un.sun_path, sizeof(sa_un.sun_path),sock_path_fmt,xmit_to_rank);

    msg.msg_name = (void *)&sa_un;
    msg.msg_namelen = sizeof(sa_un);
    msg.msg_iov = &iov;
    msg.msg_iovlen = 1;
    msg.msg_controllen =  CMSG_SPACE(sizeof(int));
    msg.msg_control = ctrl;

    cmsg = CMSG_FIRSTHDR(&msg);
    cmsg->cmsg_level = SOL_SOCKET;
    cmsg->cmsg_type = SCM_RIGHTS;
    cmsg->cmsg_len = CMSG_LEN(sizeof(int));

    *((int *) CMSG_DATA(cmsg)) = fildes;

    sendmsg(sock, &msg, 0);
  };
};
#endif


/*Construct from an MPI communicator*/
void GlobalSharedMemory::Init(Grid_MPI_Comm comm)
{
  GRID_ASSERT(_ShmSetup==0);
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
    std::cout << Mheader " World communicator of size " <<WorldSize << std::endl;  
    std::cout << Mheader " Node  communicator of size " <<WorldShmSize << std::endl;
  }
  // WorldShmComm, WorldShmSize, WorldShmRank

  // WorldNodes
  WorldNodes = WorldSize/WorldShmSize;
  GRID_ASSERT( (WorldNodes * WorldShmSize) == WorldSize );


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
      GRID_ASSERT(g<WorldShmSize);
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
  GRID_ASSERT(ierr==0);

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
  GRID_ASSERT(WorldNode!=-1);
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
void GlobalSharedMemory::OptimalCommunicator(const Coordinate &processors,Grid_MPI_Comm & optimal_comm,Coordinate &SHM)
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

  if ( TofuPlacement ) {
#ifdef GRID_FJMPI
    OptimalCommunicatorTofu(processors,optimal_comm,SHM);
#else
    if ( WorldRank==0 ) {
      std::cerr << "Grid --tofu: not built with Fujitsu MPI" << std::endl;
    }
    MPI_Barrier(WorldComm);
    MPI_Abort(WorldComm,EXIT_FAILURE);
#endif
    return;
  }
  if(nscan==3 && HPEhypercube ) OptimalCommunicatorHypercube(processors,optimal_comm,SHM);
  else                          OptimalCommunicatorSharedMemory(processors,optimal_comm,SHM);
}

void GlobalSharedMemory::OptimalCommunicatorHypercube(const Coordinate &processors,Grid_MPI_Comm & optimal_comm,Coordinate &SHM)
{
  ////////////////////////////////////////////////////////////////
  // Assert power of two shm_size.
  ////////////////////////////////////////////////////////////////
  int log2size = Log2Size(WorldShmSize,MAXLOG2RANKSPERNODE);
  GRID_ASSERT(log2size != -1);

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
  GRID_ASSERT(nscan==3);

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
  GRID_ASSERT(hypercoor<WorldSize);
  GRID_ASSERT(hypercoor>=0);

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
  SHM = ShmDims;
  
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
  GRID_ASSERT(WorldSize==Nprocessors);

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
  GRID_ASSERT(ierr==0);
}
void GlobalSharedMemory::OptimalCommunicatorSharedMemory(const Coordinate &processors,Grid_MPI_Comm & optimal_comm,Coordinate &SHM)
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
  SHM=ShmDims;

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
  //  std::cerr << " WorldSize "<<WorldSize << " Nprocessors "<<Nprocessors<<" "<<processors<<std::endl; 
  GRID_ASSERT(WorldSize==Nprocessors);

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
  GRID_ASSERT(ierr==0);
}
#ifdef GRID_FJMPI
////////////////////////////////////////////////////////////////////////////////////////////
// Tofu placement (--tofu): lay the node grid onto the job's logical torus (pjsub -L node=AxBxC:torus),
// each direction along a job axis of the same extent. Errors are detected on all ranks alike.
////////////////////////////////////////////////////////////////////////////////////////////
static void TofuPlacementError(const std::string &msg)
{
  if ( GlobalSharedMemory::WorldRank==0 ) std::cerr << "Grid --tofu: " << msg << std::endl;
  MPI_Barrier(GlobalSharedMemory::WorldComm);
  MPI_Abort(GlobalSharedMemory::WorldComm,EXIT_FAILURE);
}
static std::string TofuJoin(const std::vector<int> &v,const char *sep)
{
  std::ostringstream os;
  for(int i=0;i<v.size();i++) os << (i ? sep : "") << v[i];
  return os.str();
}
void GlobalSharedMemory::OptimalCommunicatorTofu(const Coordinate &processors,Grid_MPI_Comm & optimal_comm,Coordinate &SHM)
{
  int ndimension = processors.size();
  Coordinate WorldDims = processors; Coordinate ShmDims(ndimension);  Coordinate NodeDims (ndimension);
  Coordinate ShmCoor(ndimension);    Coordinate NodeCoor(ndimension);   Coordinate WorldCoor(ndimension);

  GetShmDims(WorldDims,ShmDims);
  SHM=ShmDims;
  int nodes=1;
  for(int d=0;d<ndimension;d++){
    NodeDims[d] = WorldDims[d]/ShmDims[d];
    nodes *= NodeDims[d];
  }
  int Nprocessors=1;
  for(int i=0;i<ndimension;i++){
    Nprocessors*=processors[i];
  }
  GRID_ASSERT(WorldSize==Nprocessors);

  std::vector<int> pvec(processors.begin(),processors.end());
  std::vector<int> nvec(NodeDims.begin(),NodeDims.end());
  std::vector<int> svec(ShmDims.begin(),ShmDims.end());
  std::string layout = "--mpi " + TofuJoin(pvec,".") + " gives a node grid " + TofuJoin(nvec,".")
    + " (" + std::to_string(nodes) + " nodes, ranks per node " + TofuJoin(svec,".") + ")";

  ////////////////////////////////////////////////////////////////
  // Shape of the job's logical node torus and this node's coordinate in it
  ////////////////////////////////////////////////////////////////
  int jdim=0;
  if ( FJMPI_Topology_get_dimension(&jdim)!=FJMPI_SUCCESS || jdim<1 || jdim>3 ) {
    TofuPlacementError("cannot read the job's node shape (FJMPI_Topology_get_dimension); "
                       "is this a Fujitsu MPI job on whole nodes?");
  }
  std::vector<int> jshape(3,1), jcoor(3,0);
  if ( FJMPI_Topology_get_shape(&jshape[0],&jshape[1],&jshape[2])!=FJMPI_SUCCESS ) {
    TofuPlacementError("cannot read the job's node shape (FJMPI_Topology_get_shape)");
  }
  for(int a=jdim;a<3;a++) jshape[a]=1;
  if ( FJMPI_Topology_get_coords(WorldComm,WorldRank,FJMPI_LOGICAL,jdim,&jcoor[0])!=FJMPI_SUCCESS ) {
    TofuPlacementError("cannot read the node coordinates (FJMPI_Topology_get_coords)");
  }
  std::vector<int> jvec(jshape.begin(),jshape.begin()+jdim);
  int jnodes = jshape[0]*jshape[1]*jshape[2];

  ////////////////////////////////////////////////////////////////
  // Lay each node grid direction of extent > 1 along an unused job axis of the same extent
  ////////////////////////////////////////////////////////////////
  std::vector<int> axis(ndimension,-1);
  std::vector<int> used(3,0);
  bool solvable = (nodes==WorldNodes) && (jnodes==WorldNodes);
  for(int d=0;d<ndimension && solvable;d++){
    if ( NodeDims[d]==1 ) continue;
    for(int a=0;a<3;a++){
      if ( !used[a] && jshape[a]==NodeDims[d] ) { axis[d]=a; used[a]=1; break; }
    }
    if ( axis[d]==-1 ) solvable=false;
  }
  for(int a=0;a<3;a++) if ( jshape[a]>1 && !used[a] ) solvable=false;
  if ( !solvable ) {
    std::vector<int> want;
    for(int d=0;d<ndimension;d++) if ( NodeDims[d]>1 ) want.push_back(NodeDims[d]);
    std::ostringstream os;
    os << layout << ", which cannot be laid onto the job's node shape " << TofuJoin(jvec,"x")
       << " (" << jnodes << " nodes, " << WorldNodes << " in use): every node grid extent larger than 1"
       << " must equal a different extent of the job shape. ";
    if ( want.size()>=1 && want.size()<=3 ) {
      while ( want.size()<3 ) want.push_back(1);
      os << "Request a job shape with the extents " << TofuJoin(want,", ") << " in any order"
         << " (pjsub -L \"node=" << TofuJoin(want,"x") << ":torus\"), or change --mpi.";
    } else {
      os << "At most 3 node grid directions may have more than one node; change --mpi.";
    }
    TofuPlacementError(os.str());
  }

  ////////////////////////////////////////////////////////////////
  // Processor coordinate: node block from the torus, on-node block as in the default layout
  ////////////////////////////////////////////////////////////////
  for(int d=0;d<ndimension;d++) NodeCoor[d] = (axis[d]>=0) ? jcoor[axis[d]] : 0;
  Lexicographic::CoorFromIndexReversed(ShmCoor ,WorldShmRank,ShmDims);
  for(int d=0;d<ndimension;d++) WorldCoor[d] = NodeCoor[d]*ShmDims[d]+ShmCoor[d];
  int rank;
  Lexicographic::IndexFromCoorReversed(WorldCoor,rank,WorldDims);

  ////////////////////////////////////////////////////////////////
  // Every processor coordinate must be taken exactly once
  ////////////////////////////////////////////////////////////////
  std::vector<int> ranks(WorldSize);
  MPI_Allgather(&rank,1,MPI_INT,&ranks[0],1,MPI_INT,WorldComm);
  std::vector<int> owner(WorldSize,-1);
  for(int r=0;r<WorldSize;r++){
    int k = ranks[r];
    if ( owner[k]!=-1 ) {
      std::ostringstream os;
      os << layout << ": world ranks " << owner[k] << " and " << r << " get the same position; the job's"
         << " logical node coordinates are not one per node, or the ranks per node differ from " << WorldShmSize;
      TofuPlacementError(os.str());
    }
    owner[k]=r;
  }

  ////////////////////////////////////////////////////////////////
  // Log the placement once per processor grid: summary, then (Debug) one line per rank
  // and the Tofu hop counts between neighbouring ranks on different nodes
  ////////////////////////////////////////////////////////////////
  static std::vector<Coordinate> logged;
  if ( std::find(logged.begin(),logged.end(),processors)==logged.end() ) {
    logged.push_back(processors);
    struct TofuRankInfo {
      int tofu[6], logical[3];
      char host[32];
    } me;
    for(int i=0;i<6;i++) me.tofu[i]=-1;
    FJMPI_Topology_get_coords(WorldComm,WorldRank,FJMPI_TOFU_SYS,6,me.tofu);
    for(int i=0;i<3;i++) me.logical[i]=jcoor[i];
    memset(me.host,0,sizeof(me.host));
    gethostname(me.host,sizeof(me.host)-1);
    std::vector<TofuRankInfo> all(WorldSize);
    MPI_Allgather(&me,sizeof(me),MPI_BYTE,&all[0],sizeof(me),MPI_BYTE,WorldComm);

    std::vector<std::string> along;
    for(int d=0;d<ndimension;d++){
      if ( axis[d]>=0 ) along.push_back("dimension " + std::to_string(d) + " along job axis " + std::to_string(axis[d]));
    }
    std::ostringstream os;
    for(int i=0;i<along.size();i++) os << (i ? ", " : "") << along[i];
    std::cout << GridLogMessage << "Tofu placement: " << layout << ", job shape " << TofuJoin(jvec,"x")
              << ": " << os.str() << std::endl;

    // Tofu hops: sum over X,Y,Z,A,B,C, B a ring of 3; approximate, X,Y,Z are taken without wrap
    auto hops = [](const TofuRankInfo &p,const TofuRankInfo &q) {
      int h=0;
      for(int i=0;i<6;i++){
        int dd = std::abs(p.tofu[i]-q.tofu[i]);
        if ( i==4 ) dd = std::min(dd,3-dd);
        h += dd;
      }
      return h;
    };
    for(int k=0;k<WorldSize;k++){
      const TofuRankInfo &r = all[owner[k]];
      Coordinate c;
      Lexicographic::CoorFromIndexReversed(c,k,WorldDims);
      std::vector<int> cv(c.begin(),c.end());
      char line[256];
      snprintf(line,sizeof(line),"rank %5d coor %-14s world %5d host %-16s tofu %2d %2d %2d %d %d %d logical %2d %2d %2d",
               k,TofuJoin(cv,".").c_str(),owner[k],r.host,r.tofu[0],r.tofu[1],r.tofu[2],r.tofu[3],r.tofu[4],r.tofu[5],
               r.logical[0],r.logical[1],r.logical[2]);
      std::cout << GridLogDebug << "Tofu placement: " << line << std::endl;
    }
    for(int d=0;d<ndimension;d++){
      if ( WorldDims[d]==1 ) continue;
      std::map<int,int> hist;
      int onnode=0;
      for(int k=0;k<WorldSize;k++){
        Coordinate c;
        Lexicographic::CoorFromIndexReversed(c,k,WorldDims);
        c[d] = (c[d]+1)%WorldDims[d];
        int n;
        Lexicographic::IndexFromCoorReversed(c,n,WorldDims);
        const TofuRankInfo &p = all[owner[k]], &q = all[owner[n]];
        if ( strncmp(p.host,q.host,sizeof(p.host))==0 ) onnode++;
        else hist[hops(p,q)]++;
      }
      std::ostringstream hs;
      for(auto &h : hist) hs << " " << h.second << " at " << h.first << (h.first==1 ? " hop;" : " hops;");
      std::cout << GridLogDebug << "Tofu placement: dimension " << d << " (" << WorldDims[d] << " ranks): "
                << onnode << " of " << WorldSize << " +1 neighbours on the same node;" << hs.str()
                << " (hops approximate, X/Y/Z not wrapped)" << std::endl;
    }
  }

  /////////////////////////////////////////////////////////////////
  // Build the new communicator
  /////////////////////////////////////////////////////////////////
  int ierr= MPI_Comm_split(WorldComm,0,rank,&optimal_comm);
  GRID_ASSERT(ierr==0);
}
#endif
////////////////////////////////////////////////////////////////////////////////////////////
// SHMGET
////////////////////////////////////////////////////////////////////////////////////////////
#ifdef GRID_MPI3_SHMGET
void GlobalSharedMemory::SharedMemoryAllocate(uint64_t bytes, int flags)
{
  std::cout << Mheader "SharedMemoryAllocate "<< bytes<< " shmget implementation "<<std::endl;
  GRID_ASSERT(_ShmSetup==1);
  GRID_ASSERT(_ShmAlloc==0);

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
void GlobalSharedMemory::SharedMemoryAllocate(uint64_t bytes, int flags)
{
  void * ShmCommBuf ; 
  GRID_ASSERT(_ShmSetup==1);
  GRID_ASSERT(_ShmAlloc==0);

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
#ifndef ACCELERATOR_AWARE_MPI
  // printf("Host buffer allocate for GPU non-aware MPI\n");
  HostCommBuf= malloc(bytes); /// CHANGE THIS TO malloc_host
#endif  
  ShmCommBuf = acceleratorAllocDevice(bytes);
  if (ShmCommBuf == (void *)NULL ) {
    std::cerr << "SharedMemoryMPI.cc acceleratorAllocDevice failed NULL pointer for " << bytes<<" bytes " << std::endl;
    exit(EXIT_FAILURE);  
  }
  if ( WorldRank == 0 ){
    std::cout << Mheader " acceleratorAllocDevice "<< bytes 
	      << "bytes at "<< std::hex<< ShmCommBuf << " - "<<(bytes-1+(uint64_t)ShmCommBuf) <<std::dec<<" for comms buffers " <<std::endl;
  }
  SharedMemoryZero(ShmCommBuf,bytes);
  if ( WorldRank == 0 ){
    std::cout<< Mheader "Setting up IPC"<<std::endl;
  }
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Loop over ranks/gpu's on our node
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef SHM_SOCKETS
  UnixSockets::Open(WorldShmRank);
#endif
  for(int r=0;r<WorldShmSize;r++){

    MPI_Barrier(WorldShmComm);

#ifndef GRID_MPI3_SHM_NONE
    //////////////////////////////////////////////////
    // If it is me, pass around the IPC access key
    //////////////////////////////////////////////////
    void * thisBuf = ShmCommBuf;
    if(!Stencil_force_mpi) {
#ifdef GRID_SYCL_LEVEL_ZERO_IPC
    typedef struct { int fd; pid_t pid ; ze_ipc_mem_handle_t ze; } clone_mem_t;

    auto zeDevice    = sycl::get_native<sycl::backend::ext_oneapi_level_zero>(theGridAccelerator->get_device());
    auto zeContext   = sycl::get_native<sycl::backend::ext_oneapi_level_zero>(theGridAccelerator->get_context());
      
    ze_ipc_mem_handle_t ihandle;
    clone_mem_t handle;
    
    if ( r==WorldShmRank ) { 
      auto err = zeMemGetIpcHandle(zeContext,ShmCommBuf,&ihandle);
      if ( err != ZE_RESULT_SUCCESS ) {
	std::cerr << "SharedMemoryMPI.cc zeMemGetIpcHandle failed for rank "<<r<<" "<<std::hex<<err<<std::dec<<std::endl;
	exit(EXIT_FAILURE);
      }
      memcpy((void *)&handle.fd,(void *)&ihandle,sizeof(int));
      handle.pid = getpid();
      memcpy((void *)&handle.ze,(void *)&ihandle,sizeof(ihandle));
#ifdef SHM_SOCKETS
      for(int rr=0;rr<WorldShmSize;rr++){
	if(rr!=r){
	  UnixSockets::SendFileDescriptor(handle.fd,rr);
	}
      }
#endif
    }
#endif
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
      MPI_Barrier(WorldShmComm);
      int ierr=MPI_Bcast(&handle,
			 sizeof(handle),
			 MPI_BYTE,
			 r,
			 WorldShmComm);
      GRID_ASSERT(ierr==0);
    }
    
    ///////////////////////////////////////////////////////////////
    // If I am not the source, overwrite thisBuf with remote buffer
    ///////////////////////////////////////////////////////////////

#ifdef GRID_SYCL_LEVEL_ZERO_IPC
    if ( r!=WorldShmRank ) {
      thisBuf = nullptr;
      int myfd;
#ifdef SHM_SOCKETS
      myfd=UnixSockets::RecvFileDescriptor();
#else
      //      std::cout<<"mapping seeking remote pid/fd "
      //	       <<handle.pid<<"/"
      //	       <<handle.fd<<std::endl;

      int pidfd = syscall(SYS_pidfd_open,handle.pid,0);
      //      std::cout<<"Using IpcHandle pidfd "<<pidfd<<"\n";
      //      int myfd  = syscall(SYS_pidfd_getfd,pidfd,handle.fd,0);
      myfd  = syscall(438,pidfd,handle.fd,0);
      int err_t = errno;
      if (myfd < 0) {
        fprintf(stderr,"pidfd_getfd returned %d errno was %d\n", myfd,err_t); fflush(stderr);
	perror("pidfd_getfd failed ");
	assert(0);
      }
#endif
      //      std::cout<<"Using IpcHandle mapped remote pid "<<handle.pid <<" FD "<<handle.fd <<" to myfd "<<myfd<<"\n";
      memcpy((void *)&ihandle,(void *)&handle.ze,sizeof(ihandle));
      memcpy((void *)&ihandle,(void *)&myfd,sizeof(int));

      auto err = zeMemOpenIpcHandle(zeContext,zeDevice,ihandle,0,&thisBuf);
      if ( err != ZE_RESULT_SUCCESS ) {
	std::cerr << "SharedMemoryMPI.cc "<<zeContext<<" "<<zeDevice<<std::endl;
	std::cerr << "SharedMemoryMPI.cc zeMemOpenIpcHandle failed for rank "<<r<<" "<<std::hex<<err<<std::dec<<std::endl; 
	exit(EXIT_FAILURE);
      }
      GRID_ASSERT(thisBuf!=nullptr);
    }
#endif
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
    }
    WorldShmCommBufs[r] = thisBuf;
#else
    WorldShmCommBufs[r] = ShmCommBuf;
#endif
    MPI_Barrier(WorldShmComm);
  }

  _ShmAllocBytes=bytes;
  _ShmAlloc=1;
}

#else 
#ifdef GRID_MPI3_SHMMMAP
void GlobalSharedMemory::SharedMemoryAllocate(uint64_t bytes, int flags)
{
  std::cout << Mheader "SharedMemoryAllocate "<< bytes<< " MMAP implementation "<< GRID_SHM_PATH <<std::endl;
  GRID_ASSERT(_ShmSetup==1);
  GRID_ASSERT(_ShmAlloc==0);
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
      perror("failed mmap");      GRID_ASSERT(0);    
    }
    GRID_ASSERT(((uint64_t)ptr&0x3F)==0);
    close(fd);
    WorldShmCommBufs[r] =ptr;
    //    std::cout << Mheader "Set WorldShmCommBufs["<<r<<"]="<<ptr<< "("<< bytes<< "bytes)"<<std::endl;
  }
  std::cout<< Mheader " Intra-node IPC setup is complete "<<std::endl;
  _ShmAlloc=1;
  _ShmAllocBytes  = bytes;
};
#endif // MMAP

#ifdef GRID_MPI3_SHM_NONE
void GlobalSharedMemory::SharedMemoryAllocate(uint64_t bytes, int flags)
{
  std::cout << Mheader "SharedMemoryAllocate "<< bytes<< " MMAP anonymous implementation "<<std::endl;
  GRID_ASSERT(_ShmSetup==1);
  GRID_ASSERT(_ShmAlloc==0);
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  // allocate the shared windows for our group
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  MPI_Barrier(WorldShmComm);
  WorldShmCommBufs.resize(WorldShmSize);
  
  ////////////////////////////////////////////////////////////////////////////////////////////
  // Hugetlbf and others map filesystems as mappable huge pages
  ////////////////////////////////////////////////////////////////////////////////////////////
  char shm_name [NAME_MAX];
  GRID_ASSERT(WorldShmSize == 1);
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
      perror("failed mmap");      GRID_ASSERT(0);    
    }
    GRID_ASSERT(((uint64_t)ptr&0x3F)==0);
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
  std::cout << Mheader "SharedMemoryAllocate "<< bytes<< " SHMOPEN implementation "<<std::endl;
  GRID_ASSERT(_ShmSetup==1);
  GRID_ASSERT(_ShmAlloc==0); 
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
      GRID_ASSERT(((uint64_t)ptr&0x3F)==0);
      
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
      if ( ptr == MAP_FAILED ) {       perror("failed mmap");      GRID_ASSERT(0);    }
      GRID_ASSERT(((uint64_t)ptr&0x3F)==0);
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
//void GlobalSharedMemory::SharedMemoryCopy(void *dest,void *src,size_t bytes)
//{
//#if defined(GRID_CUDA) || defined(GRID_HIP) || defined(GRID_SYCL)
//  acceleratorCopyToDevice(src,dest,bytes);
//#else   
//  bcopy(src,dest,bytes);
//#endif
//}
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
  GRID_ASSERT (GlobalSharedMemory::ShmAlloc()==1);
  heap_size = GlobalSharedMemory::ShmAllocBytes();
  for(int r=0;r<ShmSize;r++){

    uint32_t wsr = (r==ShmRank) ? GlobalSharedMemory::WorldShmRank : 0 ;

    MPI_Allreduce(MPI_IN_PLACE,&wsr,1,MPI_UINT32_T,MPI_SUM,ShmComm);

    ShmCommBufs[r] = GlobalSharedMemory::WorldShmCommBufs[wsr];
    //    std::cerr << " SetCommunicator rank "<<r<<" comm "<<ShmCommBufs[r] <<std::endl;
  }
  ShmBufferFreeAll();

#ifndef ACCELERATOR_AWARE_MPI
  host_heap_size = heap_size;
  HostCommBuf= GlobalSharedMemory::HostCommBuf;
  HostBufferFreeAll();
#endif  

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

  //  SharedMemoryTest();
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
       acceleratorCopyToDevice(check,ShmCommBufs[r],3*sizeof(uint64_t));
    }
  }
  ShmBarrier();
  for(uint64_t r=0;r<ShmSize;r++){
    acceleratorCopyFromDevice(ShmCommBufs[r],check,3*sizeof(uint64_t));
    GRID_ASSERT(check[0]==GlobalSharedMemory::WorldNode);
    GRID_ASSERT(check[1]==r);
    GRID_ASSERT(check[2]==magic);
  }
  ShmBarrier();
  std::cout << GridLogDebug << " SharedMemoryTest has passed "<<std::endl;
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
  GRID_ASSERT(gpeer!=ShmRank); // never send to self
  //  std::cout << "ShmBufferTranslate for rank " << rank<<" peer "<<gpeer<<std::endl;
  if (gpeer == MPI_UNDEFINED){
    return NULL;
  } else { 
    uint64_t offset = (uint64_t)local_p - (uint64_t)ShmCommBufs[ShmRank];
    uint64_t remote = (uint64_t)ShmCommBufs[gpeer]+offset;
    //    std::cout << "ShmBufferTranslate : local,offset,remote "<<std::hex<<local_p<<" "<<offset<<" "<<remote<<std::dec<<std::endl;
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

