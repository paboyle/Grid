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
#undef  ACC_HIP
#define  ACC_SYCL
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
  assert(err==cudaSuccess);
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
