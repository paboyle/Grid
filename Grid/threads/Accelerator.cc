#include <Grid/GridCore.h>

NAMESPACE_BEGIN(Grid);
uint32_t accelerator_threads=8;
uint32_t acceleratorThreads(void)       {return accelerator_threads;};
void     acceleratorThreads(uint32_t t) {accelerator_threads = t;};

#ifdef GRID_CUDA
cudaDeviceProp *gpu_props;
void acceleratorInit(void)
{
  int nDevices = 1;
  cudaGetDeviceCount(&nDevices);
  gpu_props = new cudaDeviceProp[nDevices];

  char * localRankStr = NULL;
  int rank = 0, world_rank=0; 
#define ENV_LOCAL_RANK_OMPI    "OMPI_COMM_WORLD_LOCAL_RANK"
#define ENV_LOCAL_RANK_MVAPICH "MV2_COMM_WORLD_LOCAL_RANK"
#define ENV_RANK_OMPI          "OMPI_COMM_WORLD_RANK"
#define ENV_RANK_MVAPICH       "MV2_COMM_WORLD_RANK"
  // We extract the local rank initialization using an environment variable
  if ((localRankStr = getenv(ENV_LOCAL_RANK_OMPI)) != NULL)
  {
    rank = atoi(localRankStr);		
  }
  if ((localRankStr = getenv(ENV_LOCAL_RANK_MVAPICH)) != NULL)
  {
    rank = atoi(localRankStr);		
  }
  if ((localRankStr = getenv(ENV_RANK_OMPI   )) != NULL) { world_rank = atoi(localRankStr);}
  if ((localRankStr = getenv(ENV_RANK_MVAPICH)) != NULL) { world_rank = atoi(localRankStr);}

  for (int i = 0; i < nDevices; i++) {

#define GPU_PROP_FMT(canMapHostMemory,FMT)     printf("GpuInit:   " #canMapHostMemory ": " FMT" \n",prop.canMapHostMemory);
#define GPU_PROP(canMapHostMemory)             GPU_PROP_FMT(canMapHostMemory,"%d");
    
    cudaGetDeviceProperties(&gpu_props[i], i);
    if ( world_rank == 0) {
      cudaDeviceProp prop; 
      prop = gpu_props[i];
      printf("GpuInit: ========================\n");
      printf("GpuInit: Device Number    : %d\n", i);
      printf("GpuInit: ========================\n");
      printf("GpuInit: Device identifier: %s\n", prop.name);

      GPU_PROP(managedMemory);
      GPU_PROP(isMultiGpuBoard);
      GPU_PROP(warpSize);
      //      GPU_PROP(unifiedAddressing);
      //      GPU_PROP(l2CacheSize);
      //      GPU_PROP(singleToDoublePrecisionPerfRatio);
    }
  }
#ifdef GRID_IBM_SUMMIT
  // IBM Jsrun makes cuda Device numbering screwy and not match rank
  if ( world_rank == 0 )  printf("GpuInit: IBM Summit or similar - NOT setting device to node rank\n");
#else
  if ( world_rank == 0 )  printf("GpuInit: setting device to node rank\n");
  cudaSetDevice(rank);
#endif
  if ( world_rank == 0 )  printf("GpuInit: ================================================\n");
}
#endif

#ifdef GRID_HIP
hipDeviceProp_t *gpu_props;
void acceleratorInit(void)
{
  int nDevices = 1;
  hipGetDeviceCount(&nDevices);
  gpu_props = new hipDeviceProp_t[nDevices];

  char * localRankStr = NULL;
  int rank = 0, world_rank=0; 
#define ENV_LOCAL_RANK_OMPI    "OMPI_COMM_WORLD_LOCAL_RANK"
#define ENV_LOCAL_RANK_MVAPICH "MV2_COMM_WORLD_LOCAL_RANK"
#define ENV_RANK_OMPI          "OMPI_COMM_WORLD_RANK"
#define ENV_RANK_MVAPICH       "MV2_COMM_WORLD_RANK"
  // We extract the local rank initialization using an environment variable
  if ((localRankStr = getenv(ENV_LOCAL_RANK_OMPI)) != NULL)
  {
    rank = atoi(localRankStr);		
  }
  if ((localRankStr = getenv(ENV_LOCAL_RANK_MVAPICH)) != NULL)
  {
    rank = atoi(localRankStr);		
  }
  if ((localRankStr = getenv(ENV_RANK_OMPI   )) != NULL) { world_rank = atoi(localRankStr);}
  if ((localRankStr = getenv(ENV_RANK_MVAPICH)) != NULL) { world_rank = atoi(localRankStr);}

  for (int i = 0; i < nDevices; i++) {

#define GPU_PROP_FMT(canMapHostMemory,FMT)     printf("GpuInit:   " #canMapHostMemory ": " FMT" \n",prop.canMapHostMemory);
#define GPU_PROP(canMapHostMemory)             GPU_PROP_FMT(canMapHostMemory,"%d");
    
    hipGetDeviceProperties(&gpu_props[i], i);
    if ( world_rank == 0) {
      hipDeviceProp_t prop; 
      prop = gpu_props[i];
      printf("GpuInit: ========================\n");
      printf("GpuInit: Device Number    : %d\n", i);
      printf("GpuInit: ========================\n");
      printf("GpuInit: Device identifier: %s\n", prop.name);

      //      GPU_PROP(managedMemory);
      GPU_PROP(isMultiGpuBoard);
      GPU_PROP(warpSize);
      //      GPU_PROP(unifiedAddressing);
      //      GPU_PROP(l2CacheSize);
      //      GPU_PROP(singleToDoublePrecisionPerfRatio);
    }
  }
#ifdef GRID_IBM_SUMMIT
  // IBM Jsrun makes cuda Device numbering screwy and not match rank
  if ( world_rank == 0 )  printf("GpuInit: IBM Summit or similar - NOT setting device to node rank\n");
#else
  if ( world_rank == 0 )  printf("GpuInit: setting device to node rank\n");
  cudaSetDevice(rank);
#endif
  if ( world_rank == 0 )  printf("GpuInit: ================================================\n");
}
#endif


#ifdef GRID_SYCL

cl::sycl::queue *theGridAccelerator;

void acceleratorInit(void)
{
  int nDevices = 1;
  cl::sycl::gpu_selector selector;
  cl::sycl::device selectedDevice { selector };
  theGridAccelerator = new sycl::queue (selectedDevice);

  char * localRankStr = NULL;
  int rank = 0, world_rank=0; 
#define ENV_LOCAL_RANK_OMPI    "OMPI_COMM_WORLD_LOCAL_RANK"
#define ENV_LOCAL_RANK_MVAPICH "MV2_COMM_WORLD_LOCAL_RANK"
#define ENV_RANK_OMPI          "OMPI_COMM_WORLD_RANK"
#define ENV_RANK_MVAPICH       "MV2_COMM_WORLD_RANK"
  // We extract the local rank initialization using an environment variable
  if ((localRankStr = getenv(ENV_LOCAL_RANK_OMPI)) != NULL)
  {
    rank = atoi(localRankStr);		
  }
  if ((localRankStr = getenv(ENV_LOCAL_RANK_MVAPICH)) != NULL)
  {
    rank = atoi(localRankStr);		
  }
  if ((localRankStr = getenv(ENV_RANK_OMPI   )) != NULL) { world_rank = atoi(localRankStr);}
  if ((localRankStr = getenv(ENV_RANK_MVAPICH)) != NULL) { world_rank = atoi(localRankStr);}

  if ( world_rank == 0 ) {
    GridBanner();
  }
  /*
  for (int i = 0; i < nDevices; i++) {

#define GPU_PROP_FMT(canMapHostMemory,FMT)     printf("GpuInit:   " #canMapHostMemory ": " FMT" \n",prop.canMapHostMemory);
#define GPU_PROP(canMapHostMemory)             GPU_PROP_FMT(canMapHostMemory,"%d");
    
    cudaGetDeviceProperties(&gpu_props[i], i);
    if ( world_rank == 0) {
      cudaDeviceProp prop; 
      prop = gpu_props[i];
      printf("GpuInit: ========================\n");
      printf("GpuInit: Device Number    : %d\n", i);
      printf("GpuInit: ========================\n");
      printf("GpuInit: Device identifier: %s\n", prop.name);
    }
  }
  */
  if ( world_rank == 0 ) {
    printf("GpuInit: ================================================\n");
  }
}
#endif

#if (!defined(GRID_CUDA)) && (!defined(GRID_SYCL))&& (!defined(GRID_HIP))
void acceleratorInit(void){}
#endif

NAMESPACE_END(Grid);
