#include <cassert>
#include <complex>
#include <memory>
#include <vector>
#include <algorithm>
#include <array>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <ctime>
#include <sys/time.h>

#include <mpi.h>

/**************************************************************
 * GPU - GPU memory cartesian halo exchange benchmark
 * Config: what is the target
 **************************************************************
 */
#undef ACC_CUDA
#define  ACC_HIP
#undef  ACC_SYCL
#undef  ACC_NONE

/**************************************************************
 * Some MPI globals
 **************************************************************
 */
MPI_Comm WorldComm;
MPI_Comm WorldShmComm;

int WorldSize;
int WorldRank;

int WorldShmSize;
int WorldShmRank;

/**************************************************************
 * Allocate buffers on the GPU, SYCL needs an init call and context
 **************************************************************
 */
#ifdef ACC_CUDA
#include <cuda.h>
void acceleratorInit(void){}
void *acceleratorAllocDevice(size_t bytes)
{
  void *ptr=NULL;
  auto err = cudaMalloc((void **)&ptr,bytes);
  GRID_ASSERT(err==cudaSuccess);
  return ptr;
}
void acceleratorFreeDevice(void *ptr){  cudaFree(ptr);}
#endif
#ifdef ACC_HIP
#include <hip/hip_runtime.h>
void acceleratorInit(void){}
inline void *acceleratorAllocDevice(size_t bytes)
{
  void *ptr=NULL;
  auto err = hipMalloc((void **)&ptr,bytes);
  if( err != hipSuccess ) {
    ptr = (void *) NULL;
    printf(" hipMalloc failed for %ld %s \n",bytes,hipGetErrorString(err));
  }
  return ptr;
};
inline void acceleratorFreeDevice(void *ptr){ auto r=hipFree(ptr);};
#endif
#ifdef ACC_SYCL
#include <sycl/CL/sycl.hpp>
#include <sycl/usm.hpp>
cl::sycl::queue *theAccelerator;
void acceleratorInit(void)
{
  int nDevices = 1;
#if 1
  cl::sycl::gpu_selector selector;
  cl::sycl::device selectedDevice { selector };
  theAccelerator = new sycl::queue (selectedDevice);
#else
  cl::sycl::device selectedDevice {cl::sycl::gpu_selector_v  };
  theAccelerator = new sycl::queue (selectedDevice);
#endif
  auto name = theAccelerator->get_device().get_info<sycl::info::device::name>();
  printf("AcceleratorSyclInit: Selected device is %s\n",name.c_str()); fflush(stdout);
}
inline void *acceleratorAllocDevice(size_t bytes){ return malloc_device(bytes,*theAccelerator);};
inline void acceleratorFreeDevice(void *ptr){free(ptr,*theAccelerator);};
#endif
#ifdef ACC_NONE
void acceleratorInit(void){}
inline void *acceleratorAllocDevice(size_t bytes){ return malloc(bytes);};
inline void acceleratorFreeDevice(void *ptr){free(ptr);};
#endif


/**************************************************************
 * Microsecond timer
 **************************************************************
 */
inline double usecond(void) {
  struct timeval tv;
  gettimeofday(&tv,NULL);
  return 1.0e6*tv.tv_sec + 1.0*tv.tv_usec;
}
/**************************************************************
 * Point to point cost, per Cartesian direction.
 *
 * PaddedCell::Face_exchange issues one MPI_Sendrecv per dimension and waits,
 * so the quantity that bounds a halo round is the sendrecv time for that
 * direction, not a one way ping-pong.  Sweeping the size gives the latency
 * intercept and the bandwidth slope separately.
 *
 * Each direction is labelled on- or off-node, which is what decides whether a
 * halo round rides the intra-node fabric or the network, and therefore which
 * dimension ordering costs least.  Times are reduced across ranks: the max is
 * the one that bounds a collective halo round.
 **************************************************************
 */
void PingPong(std::vector<int> cart_geom,bool use_device,int ncall)
{
  int Nd=cart_geom.size();
  std::vector<int> periodic(Nd,1);
  std::vector<int> coor(Nd);
  int rank;

  MPI_Comm communicator;
  MPI_Cart_create(WorldComm,Nd,&cart_geom[0],&periodic[0],0,&communicator);
  MPI_Comm_rank(communicator,&rank);
  MPI_Cart_coords(communicator,rank,Nd,&coor[0]);

  // An identifier every rank on a node agrees on, so a neighbour can be
  // classified by comparing it.
  int node_id = WorldRank;
  MPI_Bcast(&node_id,1,MPI_INT,0,WorldShmComm);

  size_t max_bytes = 2*1024*1024;
  void *xmit, *recv;
  if ( use_device ) {
    xmit = acceleratorAllocDevice(max_bytes);
    recv = acceleratorAllocDevice(max_bytes);
  } else {
    xmit = malloc(max_bytes);
    recv = malloc(max_bytes);
  }

  if ( !WorldRank ) {
    printf("= dim  dir       bytes      us(max)      us(min)        MB/s   off-node ranks\n");
    fflush(stdout);
  }

  for(int d=0;d<Nd;d++){

    // An undecomposed dimension is a local wrap in PaddedCell and carries no
    // message at all; a shift here would only time a rank talking to itself.
    if ( cart_geom[d] == 1 ) continue;

    for(int sign=-1;sign<=1;sign+=2){

      int from,to;
      MPI_Cart_shift(communicator,d,sign,&from,&to);

      // Fetch the identifier of the neighbour we SEND to, which means sending
      // ours the other way round the shift.
      int peer_node=node_id;
      MPI_Sendrecv(&node_id, 1,MPI_INT,from,rank,
		   &peer_node,1,MPI_INT,to,  to,
		   communicator,MPI_STATUS_IGNORE);
      int offnode = (peer_node != node_id) ? 1 : 0;
      int offnode_count=0;
      MPI_Reduce(&offnode,&offnode_count,1,MPI_INT,MPI_SUM,0,communicator);

      for(size_t bytes=8; bytes<=max_bytes; bytes*=8){

	// Latency needs the repeats; bandwidth does not, and large messages
	// would otherwise dominate the run time.
	int nc = (bytes<=32768) ? ncall : (ncall/20 > 10 ? ncall/20 : 10);

	// Untimed warm-up: the first exchange on a fresh buffer pays connection
	// setup and registration, and would otherwise land in the first row.
	MPI_Sendrecv(xmit,bytes,MPI_CHAR,to,rank,
		     recv,bytes,MPI_CHAR,from,from,
		     communicator,MPI_STATUS_IGNORE);

	MPI_Barrier(communicator);
	double t0=usecond();
	for(int i=0;i<nc;i++){
	  MPI_Sendrecv(xmit,bytes,MPI_CHAR,to,rank,
		       recv,bytes,MPI_CHAR,from,from,
		       communicator,MPI_STATUS_IGNORE);
	}
	double t1=usecond();

	double us = (t1-t0)/(double)nc;
	double us_max, us_min;
	MPI_Reduce(&us,&us_max,1,MPI_DOUBLE,MPI_MAX,0,communicator);
	MPI_Reduce(&us,&us_min,1,MPI_DOUBLE,MPI_MIN,0,communicator);

	if ( !WorldRank ) {
	  printf("= %3d  %s  %10zu  %11.2f  %11.2f  %10.1f   %d\n",
		 d, sign>0?"+":"-", bytes, us_max, us_min,
		 (double)bytes/us_max, offnode_count);
	  fflush(stdout);
	}
      }
    }
  }

  if ( use_device ) {
    acceleratorFreeDevice(xmit);
    acceleratorFreeDevice(recv);
  } else {
    free(xmit);
    free(recv);
  }
  MPI_Comm_free(&communicator);
}

/**************************************************************
 * Main benchmark routine
 **************************************************************
 */
void Benchmark(int64_t L,std::vector<int> cart_geom,bool use_device,int ncall)
{
  int64_t words = 3*4*2;
  int64_t face,vol;
  int Nd=cart_geom.size();
  
  /**************************************************************
   * L^Nd volume, L^(Nd-1) faces, 12 complex per site
   * Allocate memory for these
   **************************************************************
   */
  face=1; for( int d=0;d<Nd-1;d++) face = face*L;
  vol=1;  for( int d=0;d<Nd;d++) vol = vol*L;

  
  std::vector<void *> send_bufs;
  std::vector<void *> recv_bufs;
  size_t vw = face*words;
  size_t bytes = face*words*sizeof(double);

  if ( use_device ) {
    for(int d=0;d<2*Nd;d++){
      send_bufs.push_back(acceleratorAllocDevice(bytes));
      recv_bufs.push_back(acceleratorAllocDevice(bytes));
    }
  } else {
    for(int d=0;d<2*Nd;d++){
      send_bufs.push_back(malloc(bytes));
      recv_bufs.push_back(malloc(bytes));
    }
  }
  /*********************************************************
   * Build cartesian communicator
   *********************************************************
   */
  int ierr;
  int rank;
  std::vector<int> coor(Nd);
  MPI_Comm communicator;
  std::vector<int> periodic(Nd,1);
  MPI_Cart_create(WorldComm,Nd,&cart_geom[0],&periodic[0],0,&communicator);
  MPI_Comm_rank(communicator,&rank);
  MPI_Cart_coords(communicator,rank,Nd,&coor[0]);

  static int reported;
  if ( ! reported ) { 
    printf("World Rank %d Shm Rank %d CartCoor %d %d %d %d\n",WorldRank,WorldShmRank,
	 coor[0],coor[1],coor[2],coor[3]); fflush(stdout);
    reported =1 ;
  }
  /*********************************************************
   * Perform halo exchanges
   *********************************************************
   */
  for(int d=0;d<Nd;d++){
    if ( cart_geom[d]>1 ) {
      double t0=usecond();

      int from,to;
      
      MPI_Barrier(communicator);
      for(int n=0;n<ncall;n++){
	
	void *xmit = (void *)send_bufs[d];
	void *recv = (void *)recv_bufs[d];
	
	ierr=MPI_Cart_shift(communicator,d,1,&from,&to);
	assert(ierr==0);
	
	ierr=MPI_Sendrecv(xmit,bytes,MPI_CHAR,to,rank,
			  recv,bytes,MPI_CHAR,from, from,
			  communicator,MPI_STATUS_IGNORE);
	assert(ierr==0);
	
	xmit = (void *)send_bufs[Nd+d];
	recv = (void *)recv_bufs[Nd+d];
	
	ierr=MPI_Cart_shift(communicator,d,-1,&from,&to);
	assert(ierr==0);
	
	ierr=MPI_Sendrecv(xmit,bytes,MPI_CHAR,to,rank,
			  recv,bytes,MPI_CHAR,from, from,
			  communicator,MPI_STATUS_IGNORE);
	assert(ierr==0);
      }
      MPI_Barrier(communicator);

      double t1=usecond();
      
      double dbytes    = bytes*WorldShmSize;
      double xbytes    = dbytes*2.0*ncall;
      double rbytes    = xbytes;
      double bidibytes = xbytes+rbytes;

      if ( ! WorldRank ) {
	printf("\t%12ld\t %12ld %16.0lf\n",L,bytes,bidibytes/(t1-t0)); fflush(stdout);
      }
    }
  }
  /*********************************************************
   * Free memory
   *********************************************************
   */
  if ( use_device ) {
    for(int d=0;d<2*Nd;d++){
      acceleratorFreeDevice(send_bufs[d]);
      acceleratorFreeDevice(recv_bufs[d]);
    }
  } else {
    for(int d=0;d<2*Nd;d++){
      free(send_bufs[d]);
      free(recv_bufs[d]);
    }
  }

}

/**************************************
 * Command line junk
 **************************************/

std::string CmdOptionPayload(char ** begin, char ** end, const std::string & option)
{
  char ** itr = std::find(begin, end, option);
  if (itr != end && ++itr != end) {
    std::string payload(*itr);
    return payload;
  }
  return std::string("");
}
bool CmdOptionExists(char** begin, char** end, const std::string& option)
{
  return std::find(begin, end, option) != end;
}
void CmdOptionIntVector(const std::string &str,std::vector<int> & vec)
{
  vec.resize(0);
  std::stringstream ss(str);
  int i;
  while (ss >> i){
    vec.push_back(i);
    if(std::ispunct(ss.peek()))
      ss.ignore();
  }
  return;
}
/**************************************
 * Command line junk
 **************************************/
int main(int argc, char **argv)
{
  std::string arg;

  acceleratorInit();

  MPI_Init(&argc,&argv);

  WorldComm = MPI_COMM_WORLD;
  
  MPI_Comm_split_type(WorldComm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL,&WorldShmComm);

  MPI_Comm_rank(WorldComm     ,&WorldRank);
  MPI_Comm_size(WorldComm     ,&WorldSize);

  MPI_Comm_rank(WorldShmComm     ,&WorldShmRank);
  MPI_Comm_size(WorldShmComm     ,&WorldShmSize);

  if ( WorldSize/WorldShmSize > 2) {
    printf("This benchmark is meant to run on at most two nodes only\n");
  }

  auto mpi =std::vector<int>({1,1,1,1});

  if( CmdOptionExists(argv,argv+argc,"--mpi") ){
    arg = CmdOptionPayload(argv,argv+argc,"--mpi");
    CmdOptionIntVector(arg,mpi);
  } else {
    printf("Must specify --mpi <n1.n2.n3.n4> command line argument\n");
    exit(0);
  }

  if( !WorldRank ) {
    printf("***********************************\n");
    printf("%d ranks\n",WorldSize); 
    printf("%d ranks-per-node\n",WorldShmSize);
    printf("%d nodes\n",WorldSize/WorldShmSize);fflush(stdout);
    printf("Cartesian layout: ");
    for(int d=0;d<mpi.size();d++){
      printf("%d ",mpi[d]);
    }
    printf("\n");fflush(stdout);
    printf("***********************************\n");
  }

  
  if( !WorldRank ) {
    printf("=========================================================\n");
    printf("= Point to point cost per direction, HOST memory         \n");
    printf("=========================================================\n");fflush(stdout);
  }
  PingPong(mpi,false,1000);

  if( !WorldRank ) {
    printf("=========================================================\n");
    printf("= Point to point cost per direction, DEVICE memory       \n");
    printf("=========================================================\n");fflush(stdout);
  }
  PingPong(mpi,true,1000);

  if( !WorldRank ) {
    printf("=========================================================\n");
    printf("= Benchmarking HOST memory MPI performance               \n");
    printf("=========================================================\n");fflush(stdout);
    printf("= L\t pkt bytes\t MB/s           \n");
    printf("=========================================================\n");fflush(stdout);
  }

  for(int L=16;L<=64;L+=4){
    Benchmark(L,mpi,false,100);
  }  

  if( !WorldRank ) {
    printf("=========================================================\n");
    printf("= Benchmarking DEVICE memory MPI performance             \n");
    printf("=========================================================\n");fflush(stdout);
  }
  for(int L=16;L<=64;L+=4){
    Benchmark(L,mpi,true,100);
  }  

  if( !WorldRank ) {
    printf("=========================================================\n");
    printf("= DONE   \n");
    printf("=========================================================\n");
  }
  MPI_Finalize();
}
