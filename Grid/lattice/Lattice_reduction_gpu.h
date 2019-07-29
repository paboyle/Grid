NAMESPACE_BEGIN(Grid);

#define WARP_SIZE 32
extern cudaDeviceProp *gpu_props;
__device__ unsigned int retirementCount = 0;

template <class Iterator>
unsigned int nextPow2(Iterator x) {
  --x;
  x |= x >> 1;
  x |= x >> 2;
  x |= x >> 4;
  x |= x >> 8;
  x |= x >> 16;
  return ++x;
}

template <class Iterator>
void getNumBlocksAndThreads(const Iterator n, const size_t sizeofsobj, Iterator &threads, Iterator &blocks) {
  
  int device;
  cudaGetDevice(&device);
  
  Iterator warpSize            = gpu_props[device].warpSize;
  Iterator sharedMemPerBlock   = gpu_props[device].sharedMemPerBlock;
  Iterator maxThreadsPerBlock  = gpu_props[device].maxThreadsPerBlock;
  Iterator multiProcessorCount = gpu_props[device].multiProcessorCount;
  
  std::cout << GridLogDebug << "GPU has:" << std::endl;
  std::cout << GridLogDebug << "\twarpSize            = " << warpSize << std::endl;
  std::cout << GridLogDebug << "\tsharedMemPerBlock   = " << sharedMemPerBlock << std::endl;
  std::cout << GridLogDebug << "\tmaxThreadsPerBlock  = " << maxThreadsPerBlock << std::endl;
  std::cout << GridLogDebug << "\tmaxThreadsPerBlock  = " << warpSize << std::endl;
  std::cout << GridLogDebug << "\tmultiProcessorCount = " << multiProcessorCount << std::endl;
  
  if (warpSize != WARP_SIZE) {
    std::cout << GridLogError << "The warp size of the GPU in use does not match the warp size set when compiling Grid." << std::endl;
    exit(EXIT_FAILURE);
  }
  
  // let the number of threads in a block be a multiple of 2, starting from warpSize
  threads = warpSize;
  while( 2*threads*sizeofsobj < sharedMemPerBlock && 2*threads <= maxThreadsPerBlock ) threads *= 2;
  // keep all the streaming multiprocessors busy
  blocks = nextPow2(multiProcessorCount);
  
}

template <class sobj, class Iterator>
__device__ void reduceBlock(volatile sobj *sdata, sobj mySum, const Iterator tid) {
  
  Iterator blockSize = blockDim.x;
  
  // cannot use overloaded operators for sobj as they are not volatile-qualified
  memcpy((void *)&sdata[tid], (void *)&mySum, sizeof(sobj));
  __syncwarp();
  
  const Iterator VEC = WARP_SIZE;
  const Iterator vid = tid & (VEC-1);
  
  sobj beta, temp;
  memcpy((void *)&beta, (void *)&mySum, sizeof(sobj));
  
  for (int i = VEC/2; i > 0; i>>=1) {
    if (vid < i) {
      memcpy((void *)&temp, (void *)&sdata[tid+i], sizeof(sobj));
      beta += temp;
      memcpy((void *)&sdata[tid], (void *)&beta, sizeof(sobj));
    }
    __syncwarp();
  }
  __syncthreads();
  
  if (threadIdx.x == 0) {
    beta  = Zero();
    for (Iterator i = 0; i < blockSize; i += VEC) {
      memcpy((void *)&temp, (void *)&sdata[i], sizeof(sobj));
      beta  += temp;
    }
    memcpy((void *)&sdata[0], (void *)&beta, sizeof(sobj));
  }
  __syncthreads();
}

template <class vobj, class Iterator>
__device__ void reduceBlocks(const LatticeView<vobj> g_idata, typename vobj::scalar_object *g_odata, Iterator n) {
  
  typedef typename vobj::scalar_type scalar_type;
  typedef typename vobj::vector_type vector_type;
  typedef typename vobj::scalar_object sobj;
  constexpr Iterator nsimd = sizeof(vector_type)/sizeof(scalar_type);
  
  Iterator blockSize = blockDim.x;
  
  // force shared memory alignment
  extern __shared__ __align__(COALESCE_GRANULARITY) unsigned char shmem_pointer[];
  // it's not possible to have two extern __shared__ arrays with same name
  // but different types in different scopes -- need to cast each time
  sobj *sdata = (sobj *)shmem_pointer;
  
  // first level of reduction,
  // each thread writes result in mySum
  Iterator tid = threadIdx.x;
  Iterator i = blockIdx.x*(blockSize*2) + threadIdx.x;
  Iterator gridSize = blockSize*2*gridDim.x;
  sobj mySum = Zero();
  
  while (i < n) {
    Iterator lane = i % nsimd;
    Iterator ss   = i / nsimd;
    auto tmp = extractLane(lane,g_idata[ss]);
    mySum += tmp;
    
    if (i + blockSize < n) {
      lane = (i+blockSize) % nsimd;
      ss   = (i+blockSize) / nsimd;
      tmp = extractLane(lane,g_idata[ss]);
      mySum += tmp;
    }
    i += gridSize;
  }
  
  // copy mySum to shared memory and perform
  // reduction for all threads in this block
  reduceBlock(sdata, mySum, tid);
  if (tid == 0) g_odata[blockIdx.x] = sdata[0];
}

template <class vobj, class Iterator>
__global__ void reduceKernel(const LatticeView<vobj> lat, typename vobj::scalar_object *buffer, Iterator n) {
  typedef typename vobj::scalar_object sobj;
  
  Iterator blockSize = blockDim.x;
  
  // perform reduction for this block and
  // write result to global memory buffer
  reduceBlocks(lat, buffer, n);
  
  if (gridDim.x > 1) {
    
    const Iterator tid = threadIdx.x;
    __shared__ bool amLast;
    // force shared memory alignment
    extern __shared__ __align__(COALESCE_GRANULARITY) unsigned char shmem_pointer[];
    // it's not possible to have two extern __shared__ arrays with same name
    // but different types in different scopes -- need to cast each time
    sobj *smem = (sobj *)shmem_pointer;
    
    // wait until all outstanding memory instructions in this thread are finished
    __threadfence();
    
    if (tid==0) {
      unsigned int ticket = atomicInc(&retirementCount, gridDim.x);
      // true if this block is the last block to be done
      amLast = (ticket == gridDim.x-1);
    }
    
    // each thread must read the correct value of amLast
    __syncthreads();
    
    if (amLast) {
      // reduce buffer[0], ..., buffer[gridDim.x-1]
      Iterator i = tid;
      sobj mySum = Zero();
      
      while (i < gridDim.x) {
        mySum += buffer[i];
        i += blockSize;
      }
      
      reduceBlock(smem, mySum, tid);
      
      if (tid==0) {
        buffer[0] = smem[0];
        // reset count variable
        retirementCount = 0;
      }
    }
  }
}

template <class vobj>
inline typename vobj::scalar_object sum(const Lattice<vobj> &lat) 
{
  
  LatticeView<vobj> lat_v = lat.View();
  
  typedef typename vobj::scalar_object sobj;
  typedef decltype(lat_v.begin()) Iterator;
  
  Iterator size = lat.Grid()->lSites();
  
  Iterator numThreads, numBlocks;
  getNumBlocksAndThreads(size, sizeof(sobj), numThreads, numBlocks);
  Iterator smemSize = numThreads * sizeof(sobj);
  /*  
  std::cout << GridLogDebug << "Starting reduction with:" << std::endl;
  std::cout << GridLogDebug << "\tsize       = " << size << std::endl;
  std::cout << GridLogDebug << "\tnumThreads = " << numThreads << std::endl;
  std::cout << GridLogDebug << "\tnumBlocks  = " << numBlocks << std::endl;
  std::cout << GridLogDebug << "\tsmemSize   = " << smemSize << std::endl;
  */

  Vector<sobj> buffer(numBlocks);
  sobj *buffer_v = &buffer[0];
  
  reduceKernel<<< numBlocks, numThreads, smemSize >>>(lat_v, buffer_v, size);
  cudaDeviceSynchronize();
  
  cudaError err = cudaGetLastError();
  if ( cudaSuccess != err ) {
    printf("Cuda error %s\n",cudaGetErrorString( err ));
    exit(0);
  }
  
  auto result = buffer_v[0];
  lat.Grid()->GlobalSum(result);
  return result;
}

NAMESPACE_END(Grid);
