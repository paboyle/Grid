#include <Grid/GridCore.h>

NAMESPACE_BEGIN(Grid);
int      acceleratorAbortOnGpuError=1;
uint32_t accelerator_threads=2;
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
#define ENV_RANK_OMPI          "OMPI_COMM_WORLD_RANK"
#define ENV_LOCAL_RANK_SLURM   "SLURM_LOCALID"
#define ENV_RANK_SLURM         "SLURM_PROCID"
#define ENV_LOCAL_RANK_MVAPICH "MV2_COMM_WORLD_LOCAL_RANK"
#define ENV_RANK_MVAPICH       "MV2_COMM_WORLD_RANK"
  if ((localRankStr = getenv(ENV_RANK_OMPI   )) != NULL) { world_rank = atoi(localRankStr);}
  if ((localRankStr = getenv(ENV_RANK_MVAPICH)) != NULL) { world_rank = atoi(localRankStr);}
  if ((localRankStr = getenv(ENV_RANK_SLURM  )) != NULL) { world_rank = atoi(localRankStr);}
  // We extract the local rank initialization using an environment variable
  if ((localRankStr = getenv(ENV_LOCAL_RANK_OMPI)) != NULL) {
    if (!world_rank)
      printf("OPENMPI detected\n");
    rank = atoi(localRankStr);		
  } else if ((localRankStr = getenv(ENV_LOCAL_RANK_MVAPICH)) != NULL) {
    if (!world_rank)
      printf("MVAPICH detected\n");
    rank = atoi(localRankStr);		
  } else if ((localRankStr = getenv(ENV_LOCAL_RANK_SLURM)) != NULL) {
    if (!world_rank)
      printf("SLURM detected\n");
    rank = atoi(localRankStr);		
  } else { 
    if (!world_rank)
      printf("MPI version is unknown - bad things may happen\n");
  }

  size_t totalDeviceMem=0;
  for (int i = 0; i < nDevices; i++) {

#define GPU_PROP_FMT(canMapHostMemory,FMT)     printf("AcceleratorCudaInit[%d]:   " #canMapHostMemory ": " FMT" \n",rank,prop.canMapHostMemory);
#define GPU_PROP(canMapHostMemory)             GPU_PROP_FMT(canMapHostMemory,"%d");
    cudaGetDeviceProperties(&gpu_props[i], i);
    cudaDeviceProp prop; 
    prop = gpu_props[i];
    totalDeviceMem = prop.totalGlobalMem;
    if ( world_rank == 0) {
#ifndef GRID_DEFAULT_GPU
      if ( i==rank ) {
	printf("AcceleratorCudaInit[%d]: ========================\n",rank);
	printf("AcceleratorCudaInit[%d]: Device Number    : %d\n", rank,i);
	printf("AcceleratorCudaInit[%d]: ========================\n",rank);
	printf("AcceleratorCudaInit[%d]: Device identifier: %s\n",rank, prop.name);


	GPU_PROP_FMT(totalGlobalMem,"%lld");
	GPU_PROP(managedMemory);
	GPU_PROP(isMultiGpuBoard);
	GPU_PROP(warpSize);
	GPU_PROP(pciBusID);
	GPU_PROP(pciDeviceID);
      }
#endif
      //      GPU_PROP(unifiedAddressing);
      //      GPU_PROP(l2CacheSize);
      //      GPU_PROP(singleToDoublePrecisionPerfRatio);
    }
  }
  MemoryManager::DeviceMaxBytes = (8*totalDeviceMem)/10; // Assume 80% ours
#undef GPU_PROP_FMT    
#undef GPU_PROP

#ifdef GRID_DEFAULT_GPU
  // IBM Jsrun makes cuda Device numbering screwy and not match rank
  if ( world_rank == 0 ) {
    printf("AcceleratorCudaInit: using default device \n");
    printf("AcceleratorCudaInit: assume user either uses a) IBM jsrun, or \n");
    printf("AcceleratorCudaInit: b) invokes through a wrapping script to set CUDA_VISIBLE_DEVICES, UCX_NET_DEVICES, and numa binding \n");
    printf("AcceleratorCudaInit: Configure options --enable-summit, --enable-select-gpu=no \n");
  }
#else
  printf("AcceleratorCudaInit: rank %d setting device to node rank %d\n",world_rank,rank);
  printf("AcceleratorCudaInit: Configure options --enable-select-gpu=yes \n");
  cudaSetDevice(rank);
#endif
  if ( world_rank == 0 )  printf("AcceleratorCudaInit: ================================================\n");
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

  printf("world_rank %d has %d devices\n",world_rank,nDevices);
  size_t totalDeviceMem=0;
  for (int i = 0; i < nDevices; i++) {

#define GPU_PROP_FMT(canMapHostMemory,FMT)     printf("AcceleratorHipInit:   " #canMapHostMemory ": " FMT" \n",prop.canMapHostMemory);
#define GPU_PROP(canMapHostMemory)             GPU_PROP_FMT(canMapHostMemory,"%d");
    
    hipGetDeviceProperties(&gpu_props[i], i);
    hipDeviceProp_t prop; 
    prop = gpu_props[i];
    totalDeviceMem = prop.totalGlobalMem;
    if ( world_rank == 0) {
      printf("AcceleratorHipInit: ========================\n");
      printf("AcceleratorHipInit: Device Number    : %d\n", i);
      printf("AcceleratorHipInit: ========================\n");
      printf("AcceleratorHipInit: Device identifier: %s\n", prop.name);

      GPU_PROP_FMT(totalGlobalMem,"%lu");
      //      GPU_PROP(managedMemory);
      GPU_PROP(isMultiGpuBoard);
      GPU_PROP(warpSize);
      //      GPU_PROP(unifiedAddressing);
      //      GPU_PROP(l2CacheSize);
      //      GPU_PROP(singleToDoublePrecisionPerfRatio);
    }
  }
  MemoryManager::DeviceMaxBytes = (8*totalDeviceMem)/10; // Assume 80% ours
#undef GPU_PROP_FMT    
#undef GPU_PROP

#ifdef GRID_DEFAULT_GPU
  if ( world_rank == 0 ) {
    printf("AcceleratorHipInit: using default device \n");
    printf("AcceleratorHipInit: assume user either uses a wrapping script to set CUDA_VISIBLE_DEVICES, UCX_NET_DEVICES, and numa binding \n");
    printf("AcceleratorHipInit: Configure options --enable-summit, --enable-select-gpu=no \n");
  }
#else
  if ( world_rank == 0 ) {
    printf("AcceleratorHipInit: rank %d setting device to node rank %d\n",world_rank,rank);
    printf("AcceleratorHipInit: Configure options --enable-select-gpu=yes \n");
  }
  hipSetDevice(rank);
#endif
  if ( world_rank == 0 )  printf("AcceleratorHipInit: ================================================\n");
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

  auto devices = cl::sycl::device::get_devices();
  for(int d = 0;d<devices.size();d++){

#define GPU_PROP_STR(prop) \
    printf("AcceleratorSyclInit:   " #prop ": %s \n",devices[d].get_info<cl::sycl::info::device::prop>().c_str());

#define GPU_PROP_FMT(prop,FMT) \
    printf("AcceleratorSyclInit:   " #prop ": " FMT" \n",devices[d].get_info<cl::sycl::info::device::prop>());

#define GPU_PROP(prop)             GPU_PROP_FMT(prop,"%ld");

    GPU_PROP_STR(vendor);
    GPU_PROP_STR(version);
    //    GPU_PROP_STR(device_type);
    /*
    GPU_PROP(max_compute_units);
    GPU_PROP(native_vector_width_char);
    GPU_PROP(native_vector_width_short);
    GPU_PROP(native_vector_width_int);
    GPU_PROP(native_vector_width_long);
    GPU_PROP(native_vector_width_float);
    GPU_PROP(native_vector_width_double);
    GPU_PROP(native_vector_width_half);
    GPU_PROP(address_bits);
    GPU_PROP(half_fp_config);
    GPU_PROP(single_fp_config);
    */
    //    GPU_PROP(double_fp_config);
    GPU_PROP(global_mem_size);

  }
  if ( world_rank == 0 ) {
    auto name = theGridAccelerator->get_device().get_info<sycl::info::device::name>();
    printf("AcceleratorSyclInit: Selected device is %s\n",name.c_str());
    printf("AcceleratorSyclInit: ================================================\n");
  }
}
#endif

#if (!defined(GRID_CUDA)) && (!defined(GRID_SYCL))&& (!defined(GRID_HIP))
void acceleratorInit(void){}
#endif

NAMESPACE_END(Grid);
