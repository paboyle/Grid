/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/Accelerator.h

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: paboyle <paboyle@ph.ed.ac.uk>

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
#pragma once
NAMESPACE_BEGIN(Grid);

//////////////////////////////////////////////////////////////////////////////////
// Accelerator primitives; fall back to threading if not CUDA or SYCL
//////////////////////////////////////////////////////////////////////////////////
//
// Function attributes
//
//    accelerator
//    accelerator_inline
//
// Parallel looping
// 
//    accelerator_for
//    accelerator_forNB 
//    uint32_t accelerator_barrier();         // device synchronise
//
// Parallelism control: Number of threads in thread block is acceleratorThreads*Nsimd
//
//    uint32_t acceleratorThreads(void);   
//    void     acceleratorThreads(uint32_t);
//
// Warp control and info:
//
//    void     acceleratorSynchronise(void); // synch warp etc..
//    int      acceleratorSIMTlane(int Nsimd);
//
// Memory management:
//
//    void *acceleratorAllocShared(size_t bytes);
//    void acceleratorFreeShared(void *ptr);
//
//    void *acceleratorAllocDevice(size_t bytes);
//    void acceleratorFreeDevice(void *ptr);
//
//    void *acceleratorCopyToDevice(void *from,void *to,size_t bytes);
//    void *acceleratorCopyFromDevice(void *from,void *to,size_t bytes);
//
//////////////////////////////////////////////////////////////////////////////////

uint32_t acceleratorThreads(void);   
void     acceleratorThreads(uint32_t);

//////////////////////////////////////////////
// CUDA acceleration
//////////////////////////////////////////////
#ifdef GRID_CUDA

#ifdef __CUDA_ARCH__
#define GRID_SIMT
#endif

#define accelerator        __host__ __device__
#define accelerator_inline __host__ __device__ inline

#define accelerator_barrier(dummy)					\
  {									\
    cudaDeviceSynchronize();						\
    cudaError err = cudaGetLastError();					\
    if ( cudaSuccess != err ) {						\
      printf("Cuda error %s \n", cudaGetErrorString( err ));		\
      puts(__FILE__);							\
      printf("Line %d\n",__LINE__);					\
      exit(0);								\
    }									\
  }

#define accelerator_forNB( iterator, num, nsimd, ... )			\
  {									\
    typedef uint64_t Iterator;						\
    auto lambda = [=] accelerator (Iterator lane,Iterator iterator) mutable { \
      __VA_ARGS__;							\
    };									\
    dim3 cu_threads(acceleratorThreads(),nsimd);			\
    dim3 cu_blocks ((num+acceleratorThreads()-1)/acceleratorThreads());			\
    LambdaApply<<<cu_blocks,cu_threads>>>(nsimd,num,lambda);	\
  }

#define accelerator_for( iterator, num, nsimd, ... )		\
  accelerator_forNB(iterator, num, nsimd, { __VA_ARGS__ } );	\
  accelerator_barrier(dummy);

inline void *acceleratorAllocShared(size_t bytes)
{
  void *ptr=NULL;
  auto err = cudaMallocManaged((void **)&ptr,bytes);
  if( err != cudaSuccess ) {
    ptr = (_Tp *) NULL;
    printf(" cudaMallocManaged failed for %d %s \n",bytes,cudaGetErrorString(err));
  }
  return ptr;
};
inline void *acceleratorAllocDevice(size_t bytes)
{
  void *ptr=NULL;
  auto err = cudaMalloc((void **)&ptr,bytes);
  if( err != cudaSuccess ) {
    ptr = (_Tp *) NULL;
    printf(" cudaMalloc failed for %d %s \n",bytes,cudaGetErrorString(err));
  }
  return ptr;
};
inline void acceleratorFreeShared(void *ptr){ cudaFree(ptr);};
inline void acceleratorFreeDevice(void *ptr){ cudaFree(ptr);};

template<typename lambda>  __global__
void LambdaApply(uint64_t Isites, uint64_t Osites, lambda Lambda)
{
  uint64_t isite = threadIdx.y;
  uint64_t osite = threadIdx.x+blockDim.x*blockIdx.x;
  if ( (osite <Osites) && (isite<Isites) ) {
    Lambda(isite,osite);
  }
}

#endif

//////////////////////////////////////////////
// SyCL acceleration
//////////////////////////////////////////////

#ifdef GRID_SYCL
NAMESPACE_END(Grid);
#include <CL/sycl.hpp>
#include <CL/sycl/usm.hpp>
NAMESPACE_BEGIN(Grid);

extern cl::sycl::queue *theGridAccelerator;

#ifdef __SYCL_DEVICE_ONLY__
#define GRID_SIMT
#endif

#define accelerator 
#define accelerator_inline strong_inline

#define accelerator_forNB(iterator,num,nsimd, ... )			\
  theGridAccelerator->submit([&](cl::sycl::handler &cgh) {		\
      cl::sycl::range<3> local {acceleratorThreads(),1,nsimd};			\
      cl::sycl::range<3> global{(unsigned long)num,1,(unsigned long)nsimd}; \
      cgh.parallel_for<class dslash>(					\
      cl::sycl::nd_range<3>(global,local),            \
      [=] (cl::sycl::nd_item<3> item) mutable {       \
      auto iterator = item.get_global_id(0);	      \
      auto lane     = item.get_global_id(2);	      \
      { __VA_ARGS__ };				      \
     });	   			              \
    });

#define accelerator_barrier(dummy) theGridAccelerator->wait();

#define accelerator_for( iterator, num, nsimd, ... )		\
  accelerator_forNB(iterator, num, nsimd, { __VA_ARGS__ } );	\
  accelerator_barrier(dummy);

inline void *acceleratorAllocShared(size_t bytes){ return malloc_shared(bytes,*theGridAccelerator);};
inline void *acceleratorAllocDevice(size_t bytes){ return malloc_device(bytes,*theGridAccelerator);};
inline void acceleratorFreeShared(void *ptr){free(ptr,*theGridAccelerator);};
inline void acceleratorFreeDevice(void *ptr){free(ptr,*theGridAccelerator);};

#endif

//////////////////////////////////////////////
// HIP acceleration
//////////////////////////////////////////////
#ifdef GRID_HIP

#ifdef __HIP_DEVICE_COMPILE__
#define GRID_SIMT
#endif

#define accelerator        __host__ __device__
#define accelerator_inline __host__ __device__ inline
#define accelerator_barrier(dummy)				\
  {								\
    hipDeviceSynchronize();					\
    auto err = hipGetLastError();				\
    if ( err != hipSuccess ) {					\
      printf("HIP error %s \n", hipGetErrorString( err )); \
      puts(__FILE__); \
      printf("Line %d\n",__LINE__);					\
      exit(0);							\
    }								\
  }

#define accelerator_forNB( iterator, num, nsimd, ... )			\
  {									\
    typedef uint64_t Iterator;						\
    auto lambda = [=] accelerator (Iterator lane,Iterator iterator) mutable { \
      __VA_ARGS__;							\
    };									\
    dim3 hip_threads(acceleratorThreads(),nsimd);				\
    dim3 hip_blocks ((num+acceleratorThreads()-1)/acceleratorThreads());			\
    hipLaunchKernelGGL(LambdaApply,hip_blocks,hip_threads,0,0,num,simd,lambda);\
  }

#define accelerator_for( iterator, num, nsimd, ... )		\
  accelerator_forNB(iterator, num, nsimd, { __VA_ARGS__ } );	\
  accelerator_barrier(dummy);

inline void *acceleratorAllocShared(size_t bytes)
{
  void *ptr=NULL;
  auto err = hipMallocManaged((void **)&ptr,bytes);
  if( err != hipSuccess ) {
    ptr = (_Tp *) NULL;
    printf(" hipMallocManaged failed for %d %s \n",bytes,hipGetErrorString(err));
  }
  return ptr;
};
inline void *acceleratorAllocDevice(size_t bytes)
{
  void *ptr=NULL;
  auto err = hipMalloc((void **)&ptr,bytes);
  if( err != hipSuccess ) {
    ptr = (_Tp *) NULL;
    printf(" hipMalloc failed for %d %s \n",bytes,hipGetErrorString(err));
  }
  return ptr;
};
inline void acceleratorFreeShared(void *ptr){ hipFree(ptr);};
inline void acceleratorFreeDevice(void *ptr){ hipFree(ptr);};

template<typename lambda>  __global__
void LambdaApply(uint64_t Isites, uint64_t Osites, lambda Lambda)
{
  uint64_t isite = hipThreadIdx_y;
  uint64_t osite = hipThreadIdx_x + hipBlockDim_x*hipBlockIdx_x;
  if ( (osite <Osites) && (isite<Isites) ) {
    Lambda(isite,osite);
  }
}

#endif

//////////////////////////////////////////////
// CPU Target - No accelerator just thread instead
//////////////////////////////////////////////
#if ( (!defined(GRID_SYCL)) && (!defined(GRID_CUDA)) && (!defined(GRID_HIP)) )

#undef GRID_SIMT

#define GRID_ALLOC_ALIGN (2*1024*1024) // 2MB aligned 

#define accelerator 
#define accelerator_inline strong_inline
#define accelerator_for(iterator,num,nsimd, ... )   thread_for(iterator, num, { __VA_ARGS__ });
#define accelerator_forNB(iterator,num,nsimd, ... ) thread_for(iterator, num, { __VA_ARGS__ });
#define accelerator_barrier(dummy) 

#ifdef HAVE_MALLOC_MALLOC_H
#include <malloc/malloc.h>
#endif
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#ifdef HAVE_MM_MALLOC_H
#include <mm_malloc.h>
#endif

#ifdef HAVE_MM_MALLOC_H
inline void *acceleratorAllocShared(size_t bytes){return _mm_malloc(bytes,GRID_ALLOC_ALIGN);};
inline void *acceleratorAllocDevice(size_t bytes){return _mm_malloc(bytes,GRID_ALLOC_ALIGN);};
inline void acceleratorFreeShared(void *ptr){_mm_free(ptr);};
inline void acceleratorFreeDevice(void *ptr){_mm_free(ptr);};
#else
inline void *acceleratorAllocShared(size_t bytes){return memalign(GRID_ALLOC_ALIGN,bytes);};
inline void *acceleratorAllocDevice(size_t bytes){return memalign(GRID_ALLOC_ALIGN,bytes);};
inline void acceleratorFreeShared(void *ptr){free(ptr);};
inline void acceleratorFreeDevice(void *ptr){free(ptr);};
#endif


#endif // CPU target

///////////////////////////////////////////////////
// Synchronise across local threads for divergence resynch
///////////////////////////////////////////////////
accelerator_inline void acceleratorSynchronise(void) 
{
#ifdef GRID_SIMT
#ifdef GRID_CUDA
  __syncwarp();
#endif
#ifdef GRID_SYCL
  // No barrier call on SYCL??  // Option get __spir:: stuff to do warp barrier
#endif
#ifdef GRID_HIP
  __syncthreads();
#endif
#endif
  return;
}

////////////////////////////////////////////////////
// Address subvectors on accelerators
////////////////////////////////////////////////////
#ifdef GRID_SIMT

#ifdef GRID_CUDA
accelerator_inline int acceleratorSIMTlane(int Nsimd) { return threadIdx.y; } // CUDA specific
#endif
#ifdef GRID_SYCL
accelerator_inline int acceleratorSIMTlane(int Nsimd) { return __spirv::initLocalInvocationId<3, cl::sycl::id<3>>()[2]; } // SYCL specific
#endif
#ifdef GRID_HIP
accelerator_inline int acceleratorSIMTlane(int Nsimd) { return hipThreadIdx_y; } // HIP specific
#endif

#else

accelerator_inline int acceleratorSIMTlane(int Nsimd) { return 0; } // CUDA specific

#endif

NAMESPACE_END(Grid);
