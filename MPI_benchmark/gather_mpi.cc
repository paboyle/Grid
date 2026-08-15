#include <cassert>
#include <complex>
#include <memory>
#include <vector>
#include <algorithm>
#include <array>
#include <string>
#include <sstream>
#include <cstring>
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <ctime>
#include <sys/time.h>

#include <mpi.h>

/**************************************************************
 * Distributed dense-panel gather benchmark / reproducer.
 *
 * Pattern under test: P ranks each OWN a contiguous block of rows
 * of a (krows x ncols) fp64-complex panel; every rank must end
 * holding the WHOLE panel.  This is the setup-time collective of a
 * distributed recursive Schur inversion (row-distributed dense
 * matrices, rank-major ordering => owned rows contiguous).
 *
 * Three implementations of identical semantics:
 *   A  zero-fill + MPI_Allreduce(SUM)   (simple; non-owners send zeros)
 *   B  MPI_Allgatherv                   (owners-only send)
 *   C  ring allgather via MPI_Sendrecv  (halo-exchange primitive)
 *
 * Measured motivation (Frontier, 288 ranks, ~340MB panels): pattern A
 * on host buffers delivered 0.48 GB/s effective payload -- ~8% of its
 * own per-node wire floor, ~13x below the efficiency the Sendrecv
 * halo exchange achieves on the same NICs (185 GB/s bidirectional,
 * see halo_mpi.cc).
 *
 * Config: what is the target
 **************************************************************
 */
#undef  ACC_CUDA
#define ACC_HIP
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
  assert(err==cudaSuccess);
  return ptr;
}
void acceleratorFreeDevice(void *ptr){  cudaFree(ptr);}
void acceleratorMemSet(void *ptr,int val,size_t bytes){ cudaMemset(ptr,val,bytes);}
void acceleratorCopyToDevice(const void *from,void *to,size_t bytes){ cudaMemcpy(to,from,bytes,cudaMemcpyHostToDevice);}
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
inline void acceleratorMemSet(void *ptr,int val,size_t bytes){ auto r=hipMemset(ptr,val,bytes);};
inline void acceleratorCopyToDevice(const void *from,void *to,size_t bytes){ auto r=hipMemcpy(to,from,bytes,hipMemcpyHostToDevice);};
#endif
#ifdef ACC_SYCL
#include <sycl/CL/sycl.hpp>
#include <sycl/usm.hpp>
cl::sycl::queue *theAccelerator;
void acceleratorInit(void)
{
  cl::sycl::gpu_selector selector;
  cl::sycl::device selectedDevice { selector };
  theAccelerator = new sycl::queue (selectedDevice);
  auto name = theAccelerator->get_device().get_info<sycl::info::device::name>();
  printf("AcceleratorSyclInit: Selected device is %s\n",name.c_str()); fflush(stdout);
}
inline void *acceleratorAllocDevice(size_t bytes){ return malloc_device(bytes,*theAccelerator);};
inline void acceleratorFreeDevice(void *ptr){free(ptr,*theAccelerator);};
inline void acceleratorMemSet(void *ptr,int val,size_t bytes){ theAccelerator->memset(ptr,val,bytes).wait();};
inline void acceleratorCopyToDevice(const void *from,void *to,size_t bytes){ theAccelerator->memcpy(to,from,bytes).wait();};
#endif
#ifdef ACC_NONE
void acceleratorInit(void){}
inline void *acceleratorAllocDevice(size_t bytes){ return malloc(bytes);};
inline void acceleratorFreeDevice(void *ptr){free(ptr);};
inline void acceleratorMemSet(void *ptr,int val,size_t bytes){ memset(ptr,val,bytes);};
inline void acceleratorCopyToDevice(const void *from,void *to,size_t bytes){ memcpy(to,from,bytes);};
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
 * Main benchmark routine.
 *
 * panel = krows x ncols complex<double>, column major, krows = P*myrows.
 * Rank r owns rows [r*myrows, (r+1)*myrows): with column-major layout
 * the owned data is strided; for B and C the owners' contribution is
 * packed contiguously (rank-major block layout), which is exactly how
 * the production code would call it.  Effective payload = full panel
 * bytes; every rank must hold it all at the end.
 **************************************************************
 */
void Benchmark(size_t panel_bytes,bool use_device,int ncall)
{
  size_t words = panel_bytes/sizeof(double); // treat as flat doubles
  size_t mywords = words/WorldSize;
  words = mywords*WorldSize;                  // exact division
  size_t bytes  = words*sizeof(double);
  size_t mybytes= mywords*sizeof(double);

  double *panel;
  double *contrib;
  if ( use_device ) {
    panel   = (double *)acceleratorAllocDevice(bytes);
    contrib = (double *)acceleratorAllocDevice(mybytes);
    if ( panel==NULL || contrib==NULL ) { printf("alloc failed\n"); return; }
  } else {
    panel   = (double *)malloc(bytes);
    contrib = (double *)malloc(mybytes);
  }
  // Owners deposit non-trivial data once (content is irrelevant to timing)
  std::vector<double> init(mywords,1.0*WorldRank);
  acceleratorCopyToDevice(&init[0],contrib,mybytes);
  if ( !use_device ) memcpy(contrib,&init[0],mybytes);

  double tA,tB,tC;

  /*********************************************************
   * A: zero-fill + Allreduce(SUM) -- the simple idiom
   *********************************************************/
  {
    MPI_Barrier(WorldComm);
    double t0=usecond();
    for(int n=0;n<ncall;n++){
      if ( use_device ) acceleratorMemSet(panel,0,bytes);
      else              memset(panel,0,bytes);
      // owner deposits its block at its rank-major offset
      if ( use_device ) {
        // device-to-device deposit stands in for the real kernel
#ifdef ACC_HIP
        auto r=hipMemcpy(panel+WorldRank*mywords,contrib,mybytes,hipMemcpyDeviceToDevice);
#else
        acceleratorCopyToDevice(contrib,panel+WorldRank*mywords,mybytes);
#endif
      } else {
        memcpy(panel+WorldRank*mywords,contrib,mybytes);
      }
      int ierr=MPI_Allreduce(MPI_IN_PLACE,panel,words,MPI_DOUBLE,MPI_SUM,WorldComm);
      assert(ierr==0);
    }
    MPI_Barrier(WorldComm);
    tA=(usecond()-t0)/ncall;
  }

  /*********************************************************
   * A2: as A but with a FRESH buffer every call -- the
   * registration-cache killer (unpinned pages each call).
   * Host-memory case only; isolates the production bug.
   *********************************************************/
  double tA2 = 0.0;
  if ( !use_device ) {
    MPI_Barrier(WorldComm);
    double t0=usecond();
    for(int n=0;n<ncall;n++){
      double *fresh = (double *)malloc(bytes);
      memset(fresh,0,bytes);
      memcpy(fresh+WorldRank*mywords,contrib,mybytes);
      int ierr=MPI_Allreduce(MPI_IN_PLACE,fresh,words,MPI_DOUBLE,MPI_SUM,WorldComm);
      assert(ierr==0);
      free(fresh);
    }
    MPI_Barrier(WorldComm);
    tA2=(usecond()-t0)/ncall;
  }

  /*********************************************************
   * D: chunked Allreduce on the persistent buffer -- hand
   * decomposition of the primitive (the split-K analogue):
   * does MPI fail to pipeline internally?
   *********************************************************/
  double tD;
  {
    size_t chunkwords = (32ULL<<20)/sizeof(double);   // 32MB chunks
    MPI_Barrier(WorldComm);
    double t0=usecond();
    for(int n=0;n<ncall;n++){
      if ( use_device ) acceleratorMemSet(panel,0,bytes);
      else              memset(panel,0,bytes);
      if ( use_device ) {
#ifdef ACC_HIP
        auto r=hipMemcpy(panel+WorldRank*mywords,contrib,mybytes,hipMemcpyDeviceToDevice);
#else
        acceleratorCopyToDevice(contrib,panel+WorldRank*mywords,mybytes);
#endif
      } else {
        memcpy(panel+WorldRank*mywords,contrib,mybytes);
      }
      for(size_t off=0; off<words; off+=chunkwords){
        size_t cw = std::min(chunkwords, words-off);
        int ierr=MPI_Allreduce(MPI_IN_PLACE,panel+off,(int)cw,MPI_DOUBLE,MPI_SUM,WorldComm);
        assert(ierr==0);
      }
    }
    MPI_Barrier(WorldComm);
    tD=(usecond()-t0)/ncall;
  }

  /*********************************************************
   * B: Allgatherv -- owners-only send, identical result
   *********************************************************/
  {
    std::vector<int> counts(WorldSize,(int)mywords);
    std::vector<int> displs(WorldSize);
    for(int r=0;r<WorldSize;r++) displs[r]=r*(int)mywords;

    MPI_Barrier(WorldComm);
    double t0=usecond();
    for(int n=0;n<ncall;n++){
      int ierr=MPI_Allgatherv(contrib,(int)mywords,MPI_DOUBLE,
                              panel,&counts[0],&displs[0],MPI_DOUBLE,WorldComm);
      assert(ierr==0);
    }
    MPI_Barrier(WorldComm);
    tB=(usecond()-t0)/ncall;
  }

  /*********************************************************
   * C: ring allgather from Sendrecv -- the halo primitive.
   * P-1 steps; step s passes the block received at step s-1 to
   * the right neighbour.  Per-step message = panel/P.
   *********************************************************/
  {
    int right=(WorldRank+1)%WorldSize;
    int left =(WorldRank-1+WorldSize)%WorldSize;

    MPI_Barrier(WorldComm);
    double t0=usecond();
    for(int n=0;n<ncall;n++){
      // start with my own block resident at my slot
      if ( use_device ) {
#ifdef ACC_HIP
        auto r=hipMemcpy(panel+WorldRank*mywords,contrib,mybytes,hipMemcpyDeviceToDevice);
#else
        acceleratorCopyToDevice(contrib,panel+WorldRank*mywords,mybytes);
#endif
      } else {
        memcpy(panel+WorldRank*mywords,contrib,mybytes);
      }
      for(int s=0;s<WorldSize-1;s++){
        int sendblock=(WorldRank-s+WorldSize)%WorldSize;
        int recvblock=(WorldRank-s-1+WorldSize)%WorldSize;
        int ierr=MPI_Sendrecv(panel+sendblock*mywords,(int)mywords,MPI_DOUBLE,right,s,
                              panel+recvblock*mywords,(int)mywords,MPI_DOUBLE,left, s,
                              WorldComm,MPI_STATUS_IGNORE);
        assert(ierr==0);
      }
    }
    MPI_Barrier(WorldComm);
    tC=(usecond()-t0)/ncall;
  }

  if ( !WorldRank ) {
    double GB=bytes/1024./1024./1024.;
    printf("\t%10.3f GB\t A allreduce %8.2f\t A2 fresh-buf %8.2f\t D chunked %8.2f\t B allgatherv %8.2f\t C ring %8.2f  GB/s\n",
           GB,
           GB/(tA/1.0e6),
           (tA2>0.0) ? GB/(tA2/1.0e6) : 0.0,
           GB/(tD/1.0e6),
           GB/(tB/1.0e6),
           GB/(tC/1.0e6));
    fflush(stdout);
  }

  if ( use_device ) {
    acceleratorFreeDevice(panel);
    acceleratorFreeDevice(contrib);
  } else {
    free(panel);
    free(contrib);
  }
}

/**************************************
 * Command line junk
 **************************************/
int main(int argc, char **argv)
{
  acceleratorInit();

  MPI_Init(&argc,&argv);

  WorldComm = MPI_COMM_WORLD;

  MPI_Comm_split_type(WorldComm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL,&WorldShmComm);

  MPI_Comm_rank(WorldComm     ,&WorldRank);
  MPI_Comm_size(WorldComm     ,&WorldSize);

  MPI_Comm_rank(WorldShmComm     ,&WorldShmRank);
  MPI_Comm_size(WorldShmComm     ,&WorldShmSize);

  if( !WorldRank ) {
    printf("***********************************\n");
    printf("%d ranks\n",WorldSize);
    printf("%d ranks-per-node\n",WorldShmSize);
    printf("%d nodes\n",WorldSize/WorldShmSize);fflush(stdout);
    printf("***********************************\n");
    printf("Panel gather: every rank owns 1/%d of the rows;\n",WorldSize);
    printf("every rank must finish holding the whole panel.\n");
    printf("Effective payload rate = panel bytes / wall.  Same semantics x3:\n");
    printf("  A zero-fill+Allreduce   B Allgatherv   C ring Sendrecv\n");
  }

  std::vector<size_t> sizes({ (size_t)8<<20, (size_t)64<<20, (size_t)256<<20, (size_t)1<<30 });

  if( !WorldRank ) {
    printf("=========================================================\n");
    printf("= HOST memory                                            \n");
    printf("=========================================================\n");fflush(stdout);
  }
  for( auto sz : sizes ) Benchmark(sz,false,5);

  if( !WorldRank ) {
    printf("=========================================================\n");
    printf("= DEVICE memory                                          \n");
    printf("=========================================================\n");fflush(stdout);
  }
  for( auto sz : sizes ) Benchmark(sz,true,5);

  if( !WorldRank ) {
    printf("=========================================================\n");
    printf("= DONE   \n");
    printf("=========================================================\n");
  }
  MPI_Finalize();
}
