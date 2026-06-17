NAMESPACE_BEGIN(Grid);

#ifdef GRID_HIP
extern hipDeviceProp_t *gpu_props;
#define WARP_SIZE 64
#endif
#ifdef GRID_CUDA
extern cudaDeviceProp *gpu_props;
#define WARP_SIZE 32
#endif

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
int getNumBlocksAndThreads(const Iterator n, const size_t sizeofsobj, Iterator &threads, Iterator &blocks) {
  
  int device;
#ifdef GRID_CUDA
  cudaGetDevice(&device);
#endif
#ifdef GRID_HIP
  auto r=hipGetDevice(&device);
#endif
  
  Iterator warpSize            = gpu_props[device].warpSize;
  Iterator sharedMemPerBlock   = gpu_props[device].sharedMemPerBlock;
  Iterator maxThreadsPerBlock  = gpu_props[device].maxThreadsPerBlock;
  Iterator multiProcessorCount = gpu_props[device].multiProcessorCount;
  /*  
  std::cout << GridLogDebug << "GPU has:" << std::endl;
  std::cout << GridLogDebug << "\twarpSize            = " << warpSize << std::endl;
  std::cout << GridLogDebug << "\tsharedMemPerBlock   = " << sharedMemPerBlock << std::endl;
  std::cout << GridLogDebug << "\tmaxThreadsPerBlock  = " << maxThreadsPerBlock << std::endl;
  std::cout << GridLogDebug << "\tmultiProcessorCount = " << multiProcessorCount << std::endl;
  */  
  if (warpSize != WARP_SIZE) {
    std::cout << GridLogError << "The warp size of the GPU in use does not match the warp size set when compiling Grid." << std::endl;
    exit(EXIT_FAILURE);
  }
  
  // let the number of threads in a block be a multiple of 2, starting from warpSize
  threads = warpSize;
  if ( threads*sizeofsobj > sharedMemPerBlock ) {
    std::cout << GridLogError << "The object is too large for the shared memory." << std::endl;
    return 0;
  }
  while( 2*threads*sizeofsobj < sharedMemPerBlock && 2*threads <= maxThreadsPerBlock ) threads *= 2;
  // keep all the streaming multiprocessors busy
  blocks = nextPow2(multiProcessorCount);
  return 1;
}

template <class sobj, class Iterator>
__device__ void reduceBlock(volatile sobj *sdata, sobj mySum, const Iterator tid) {
  
  Iterator blockSize = blockDim.x;
  
  // cannot use overloaded operators for sobj as they are not volatile-qualified
  memcpy((void *)&sdata[tid], (void *)&mySum, sizeof(sobj));
  acceleratorSynchronise();
  
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
    acceleratorSynchronise();
  }
  acceleratorSynchroniseAll();
  
  if (threadIdx.x == 0) {
    beta  = Zero();
    for (Iterator i = 0; i < blockSize; i += VEC) {
      memcpy((void *)&temp, (void *)&sdata[i], sizeof(sobj));
      beta  += temp;
    }
    memcpy((void *)&sdata[0], (void *)&beta, sizeof(sobj));
  }
  acceleratorSynchroniseAll();
}


template <class vobj, class sobj, class Iterator>
__device__ void reduceBlocks(const vobj *g_idata, sobj *g_odata, Iterator n) 
{
  constexpr Iterator nsimd = vobj::Nsimd();
  
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
    sobj tmpD;
    tmpD=tmp;
    mySum   +=tmpD;
    
    if (i + blockSize < n) {
      lane = (i+blockSize) % nsimd;
      ss   = (i+blockSize) / nsimd;
      tmp = extractLane(lane,g_idata[ss]);
      tmpD = tmp;
      mySum += tmpD;
    }
    i += gridSize;
  }
  
  // copy mySum to shared memory and perform
  // reduction for all threads in this block
  reduceBlock(sdata, mySum, tid);
  if (tid == 0) g_odata[blockIdx.x] = sdata[0];
}

template <class vobj, class sobj,class Iterator>
__global__ void reduceKernel(const vobj *lat, sobj *buffer, Iterator n) {
  
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
    acceleratorFence();
    
    if (tid==0) {
      unsigned int ticket = atomicInc(&retirementCount, gridDim.x);
      // true if this block is the last block to be done
      amLast = (ticket == gridDim.x-1);
    }
    
    // each thread must read the correct value of amLast
    acceleratorSynchroniseAll();

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

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Possibly promote to double and sum
/////////////////////////////////////////////////////////////////////////////////////////////////////////

#undef GRID_REDUCTION_TIMING

template <class vobj>
inline typename vobj::scalar_objectD sumD_gpu_small(const vobj *lat, Integer osites)
{
  typedef typename vobj::scalar_objectD sobj;
  typedef decltype(lat) Iterator;

  Integer nsimd= vobj::Nsimd();
  Integer size = osites*nsimd;

  Integer numThreads, numBlocks;
  int ok = getNumBlocksAndThreads(size, sizeof(sobj), numThreads, numBlocks);
  GRID_ASSERT(ok);

  Integer smemSize = numThreads * sizeof(sobj);
  deviceVector<sobj> buffer(numBlocks);
  sobj *buffer_v = &buffer[0];
  sobj result;

#ifdef GRID_REDUCTION_TIMING
  RealD t_kernel = -usecond();
#endif
  reduceKernel<<< numBlocks, numThreads, smemSize, computeStream >>>(lat, buffer_v, size);
  accelerator_barrier();
#ifdef GRID_REDUCTION_TIMING
  t_kernel += usecond();
  RealD t_d2h = -usecond();
#endif
  acceleratorCopyFromDevice(buffer_v,&result,sizeof(result));
#ifdef GRID_REDUCTION_TIMING
  t_d2h += usecond();
  std::cout << GridLogDebug << "  sumD_gpu_small"
            << " sizeof(sobj)=" << sizeof(sobj)
            << " blocks=" << numBlocks << " threads=" << numThreads
            << " kernel+barrier=" << t_kernel << " us"
            << " D2H=" << t_d2h << " us" << std::endl;
#endif
  return result;
}

// Fused pack+reduce: reads R words of each vobj at word offset 'base',
// accumulates directly into iVector<iScalar<scalarD>,R> without staging
// through an intermediate bundle buffer.  One HBM pass instead of three.
template <int R, class vobj, class sobj, class Iterator>
__device__ void packReduceBlocks(
    const iScalar<typename vobj::vector_type> *idat,
    sobj *g_odata, Iterator osites, int base, int words)
{
  constexpr Iterator nsimd = vobj::Nsimd();
  Iterator blockSize = blockDim.x;

  extern __shared__ __align__(COALESCE_GRANULARITY) unsigned char shmem_pointer[];
  sobj *sdata = (sobj *)shmem_pointer;

  Iterator tid      = threadIdx.x;
  Iterator i        = blockIdx.x * (blockSize * 2) + threadIdx.x;
  Iterator gridSize = blockSize * 2 * gridDim.x;
  sobj mySum = Zero();

  while (i < osites * nsimd) {
    Iterator lane = i % nsimd;
    Iterator ss   = i / nsimd;
    sobj tmpD; zeroit(tmpD);
    for (int k = 0; k < R; k++) {
      auto w = extractLane(lane, idat[ss * words + base + k]);
      iScalar<typename vobj::scalar_typeD> wd; wd = w;
      tmpD._internal[k] = wd;
    }
    mySum += tmpD;

    if (i + blockSize < osites * nsimd) {
      lane = (i + blockSize) % nsimd;
      ss   = (i + blockSize) / nsimd;
      sobj tmpD2; zeroit(tmpD2);
      for (int k = 0; k < R; k++) {
        auto w = extractLane(lane, idat[ss * words + base + k]);
        iScalar<typename vobj::scalar_typeD> wd; wd = w;
        tmpD2._internal[k] = wd;
      }
      mySum += tmpD2;
    }
    i += gridSize;
  }

  reduceBlock(sdata, mySum, tid);
  if (tid == 0) g_odata[blockIdx.x] = sdata[0];
}

template <int R, class vobj, class sobj, class Iterator>
__global__ void packReduceKernel(
    const iScalar<typename vobj::vector_type> *idat,
    sobj *buffer, Iterator osites, int base, int words)
{
  Iterator blockSize = blockDim.x;

  packReduceBlocks<R, vobj, sobj>(idat, buffer, osites, base, words);

  if (gridDim.x > 1) {
    const Iterator tid = threadIdx.x;
    __shared__ bool amLast;
    extern __shared__ __align__(COALESCE_GRANULARITY) unsigned char shmem_pointer[];
    sobj *smem = (sobj *)shmem_pointer;

    acceleratorFence();

    if (tid == 0) {
      unsigned int ticket = atomicInc(&retirementCount, gridDim.x);
      amLast = (ticket == gridDim.x - 1);
    }
    acceleratorSynchroniseAll();

    if (amLast) {
      Iterator i = tid;
      sobj mySum = Zero();
      while (i < (Iterator)gridDim.x) {
        mySum += buffer[i];
        i += blockSize;
      }
      reduceBlock(smem, mySum, tid);
      if (tid == 0) {
        buffer[0] = smem[0];
        retirementCount = 0;
      }
    }
  }
}

template<int R, class vobj>
inline void sumD_gpu_reduce_words(const vobj *lat, Integer osites,
                                   typename vobj::scalar_typeD *ret_p, int base)
{
  typedef typename vobj::vector_type  vector;
  typedef typename vobj::scalar_typeD scalarD;
  using BundleScalarD = iVector<iScalar<scalarD>, R>;

  constexpr int Nsimd = vobj::Nsimd();
  const int words = sizeof(vobj) / sizeof(vector);
  const iScalar<vector> *idat = (const iScalar<vector> *)lat;

  Integer size = (Integer)osites * Nsimd;
  Integer numThreads, numBlocks;
  int ok = getNumBlocksAndThreads(size, sizeof(BundleScalarD), numThreads, numBlocks);
  GRID_ASSERT(ok);

  Integer smemSize = numThreads * sizeof(BundleScalarD);
  deviceVector<BundleScalarD> buffer(numBlocks);
  BundleScalarD *buffer_v = &buffer[0];
  BundleScalarD result;

#ifdef GRID_REDUCTION_TIMING
  RealD t_kernel = -usecond();
#endif
  packReduceKernel<R, vobj, BundleScalarD, Integer>
    <<<numBlocks, numThreads, smemSize, computeStream>>>
    (idat, buffer_v, osites, base, words);
  accelerator_barrier();
#ifdef GRID_REDUCTION_TIMING
  t_kernel += usecond();
  RealD t_d2h = -usecond();
#endif
  acceleratorCopyFromDevice(buffer_v, &result, sizeof(result));
#ifdef GRID_REDUCTION_TIMING
  t_d2h += usecond();
  std::cout << GridLogDebug << " sumD_gpu_reduce_words R=" << R
            << " base=" << base
            << " kernel=" << t_kernel << " D2H=" << t_d2h << " us" << std::endl;
#endif

  for (int k = 0; k < R; k++)
    ret_p[base + k] = TensorRemove(result._internal[k]);
}

template <class vobj>
inline typename vobj::scalar_objectD sumD_gpu_large(const vobj *lat, Integer osites)
{
  typedef typename vobj::vector_type    vector;
  typedef typename vobj::scalar_typeD   scalarD;
  typedef typename vobj::scalar_objectD sobjD;

  const int words = sizeof(vobj) / sizeof(vector);
  sobjD ret; zeroit(ret);
  scalarD *ret_p = (scalarD *)&ret;

#ifdef GRID_REDUCTION_TIMING
  RealD t_large = -usecond();
#endif
  int w = 0;
  while (w + 12 <= words) { sumD_gpu_reduce_words<12>(lat, osites, ret_p, w); w += 12; }
  while (w +  4 <= words) { sumD_gpu_reduce_words< 4>(lat, osites, ret_p, w); w +=  4; }
  while (w      <  words) { sumD_gpu_reduce_words< 1>(lat, osites, ret_p, w); w +=  1; }
#ifdef GRID_REDUCTION_TIMING
  t_large += usecond();
  std::cout << GridLogDebug << "sumD_gpu_large"
            << " sizeof(sobjD)=" << sizeof(sobjD)
            << " words=" << words << " total=" << t_large << " us" << std::endl;
#endif

  return ret;
}

template <class vobj>
inline typename vobj::scalar_objectD sumD_gpu(const vobj *lat, Integer osites)
{
  typedef typename vobj::scalar_objectD sobj;
  sobj ret;
  
  Integer nsimd= vobj::Nsimd();
  Integer size = osites*nsimd;
  Integer numThreads, numBlocks;
  int ok = getNumBlocksAndThreads(size, sizeof(sobj), numThreads, numBlocks);
  
  if ( ok ) {
    ret = sumD_gpu_small(lat,osites);
  } else {
    ret = sumD_gpu_large(lat,osites);
  }
  return ret;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Return as same precision as input performing reduction in double precision though
/////////////////////////////////////////////////////////////////////////////////////////////////////////
template <class vobj>
inline typename vobj::scalar_object sum_gpu(const vobj *lat, Integer osites) 
{
  typedef typename vobj::scalar_object sobj;
  sobj result;
  result = sumD_gpu(lat,osites);
  return result;
}

template <class vobj>
inline typename vobj::scalar_object sum_gpu_large(const vobj *lat, Integer osites)
{
  typedef typename vobj::scalar_object sobj;
  sobj result;
  result = sumD_gpu_large(lat,osites);
  return result;
}
template<class Word> Word checksum_gpu(Word *vec,uint64_t L)
{
  Word w;
  bzero(&w,sizeof(w));
  return w;
}

NAMESPACE_END(Grid);
