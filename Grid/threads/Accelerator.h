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

#include <string.h>

#ifdef HAVE_MALLOC_MALLOC_H
#include <malloc/malloc.h>
#endif
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#ifdef HAVE_MM_MALLOC_H
#include <mm_malloc.h>
#endif

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
//    acceleratorInit;
//    void     acceleratorSynchronise(void); // synch warp etc..
//    int      acceleratorSIMTlane(int Nsimd);
//
// Memory management:
//
//    int   acceleratorIsCommunicable(void *pointer);
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
void     acceleratorInit(void);

//////////////////////////////////////////////
// CUDA acceleration
//////////////////////////////////////////////

#ifdef GRID_CUDA
#include <cuda.h>

#ifdef __CUDA_ARCH__
#define GRID_SIMT
#endif

#define accelerator        __host__ __device__
#define accelerator_inline __host__ __device__ inline

extern int acceleratorAbortOnGpuError;

accelerator_inline int acceleratorSIMTlane(int Nsimd) {
#ifdef GRID_SIMT
  return threadIdx.x; 
#else
  return 0;
#endif
} // CUDA specific

#define accelerator_for2dNB( iter1, num1, iter2, num2, nsimd, ... )	\
  {									\
    int nt=acceleratorThreads();					\
    typedef uint64_t Iterator;						\
    auto lambda = [=] accelerator					\
      (Iterator iter1,Iterator iter2,Iterator lane) mutable {		\
      __VA_ARGS__;							\
    };									\
    dim3 cu_threads(nsimd,acceleratorThreads(),1);			\
    dim3 cu_blocks ((num1+nt-1)/nt,num2,1);				\
    LambdaApply<<<cu_blocks,cu_threads>>>(num1,num2,nsimd,lambda);	\
  }

#define accelerator_for6dNB(iter1, num1,				\
                            iter2, num2,				\
                            iter3, num3,				\
                            iter4, num4,				\
                            iter5, num5,				\
			    iter6, num6, ... )				\
  {									\
    typedef uint64_t Iterator;						\
    auto lambda = [=] accelerator					\
      (Iterator iter1,Iterator iter2,					\
       Iterator iter3,Iterator iter4,					\
       Iterator iter5,Iterator iter6) mutable {				\
      __VA_ARGS__;							\
    };									\
    dim3 cu_blocks (num1,num2,num3);					\
    dim3 cu_threads(num4,num5,num6);					\
    Lambda6Apply<<<cu_blocks,cu_threads>>>(num1,num2,num3,num4,num5,num6,lambda); \
  }

template<typename lambda>  __global__
void LambdaApply(uint64_t num1, uint64_t num2, uint64_t num3, lambda Lambda)
{
  // Weird permute is to make lane coalesce for large blocks
  uint64_t x = threadIdx.y + blockDim.y*blockIdx.x;
  uint64_t y = threadIdx.z + blockDim.z*blockIdx.y;
  uint64_t z = threadIdx.x;
  if ( (x < num1) && (y<num2) && (z<num3) ) {
    Lambda(x,y,z);
  }
}

template<typename lambda>  __global__
void Lambda6Apply(uint64_t num1, uint64_t num2, uint64_t num3,
		  uint64_t num4, uint64_t num5, uint64_t num6,
		  lambda Lambda)
{
  uint64_t iter1 = blockIdx.x;
  uint64_t iter2 = blockIdx.y;
  uint64_t iter3 = blockIdx.z;
  uint64_t iter4 = threadIdx.x;
  uint64_t iter5 = threadIdx.y;
  uint64_t iter6 = threadIdx.z;

  if ( (iter1 < num1) && (iter2<num2) && (iter3<num3)
    && (iter4 < num4) && (iter5<num5) && (iter6<num6) )
  {
    Lambda(iter1,iter2,iter3,iter4,iter5,iter6);
  }
}

#define accelerator_barrier(dummy)					\
  {									\
    cudaDeviceSynchronize();						\
    cudaError err = cudaGetLastError();					\
    if ( cudaSuccess != err ) {						\
      printf("accelerator_barrier(): Cuda error %s \n",			\
	     cudaGetErrorString( err ));				\
      printf("File %s Line %d\n",__FILE__,__LINE__);			\
      fflush(stdout);							\
      if (acceleratorAbortOnGpuError) assert(err==cudaSuccess);		\
    }									\
  }

inline void *acceleratorAllocShared(size_t bytes)
{
  void *ptr=NULL;
  auto err = cudaMallocManaged((void **)&ptr,bytes);
  if( err != cudaSuccess ) {
    ptr = (void *) NULL;
    printf(" cudaMallocManaged failed for %d %s \n",bytes,cudaGetErrorString(err));
  }
  return ptr;
};
inline void *acceleratorAllocDevice(size_t bytes)
{
  void *ptr=NULL;
  auto err = cudaMalloc((void **)&ptr,bytes);
  if( err != cudaSuccess ) {
    ptr = (void *) NULL;
    printf(" cudaMalloc failed for %d %s \n",bytes,cudaGetErrorString(err));
  }
  return ptr;
};
inline void acceleratorFreeShared(void *ptr){ cudaFree(ptr);};
inline void acceleratorFreeDevice(void *ptr){ cudaFree(ptr);};
inline void acceleratorCopyToDevice(void *from,void *to,size_t bytes)  { cudaMemcpy(to,from,bytes, cudaMemcpyHostToDevice);}
inline void acceleratorCopyDeviceToDevice(void *from,void *to,size_t bytes)  { cudaMemcpy(to,from,bytes, cudaMemcpyDeviceToDevice);}
inline void acceleratorCopyFromDevice(void *from,void *to,size_t bytes){ cudaMemcpy(to,from,bytes, cudaMemcpyDeviceToHost);}
inline void acceleratorMemSet(void *base,int value,size_t bytes) { cudaMemset(base,value,bytes);}
inline int  acceleratorIsCommunicable(void *ptr)
{
  //  int uvm=0;
  //  auto 
  //  cuerr = cuPointerGetAttribute( &uvm, CU_POINTER_ATTRIBUTE_IS_MANAGED, (CUdeviceptr) ptr);
  //  assert(cuerr == cudaSuccess );
  //  if(uvm) return 0;
  //  else    return 1;
    return 1;
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

accelerator_inline int acceleratorSIMTlane(int Nsimd) {
#ifdef GRID_SIMT
 return __spirv::initLocalInvocationId<3, cl::sycl::id<3>>()[2]; 
#else
 return 0;
#endif
} // SYCL specific

#define accelerator_for2dNB( iter1, num1, iter2, num2, nsimd, ... )	\
  theGridAccelerator->submit([&](cl::sycl::handler &cgh) {		\
      unsigned long nt=acceleratorThreads();				\
      unsigned long unum1 = num1;					\
      unsigned long unum2 = num2;					\
      if(nt < 8)nt=8;							\
      cl::sycl::range<3> local {nt,1,nsimd};				\
      cl::sycl::range<3> global{unum1,unum2,nsimd};			\
      cgh.parallel_for<class dslash>(					\
      cl::sycl::nd_range<3>(global,local), \
      [=] (cl::sycl::nd_item<3> item) /*mutable*/     \
      [[intel::reqd_sub_group_size(8)]]	      \
      {						      \
      auto iter1    = item.get_global_id(0);	      \
      auto iter2    = item.get_global_id(1);	      \
      auto lane     = item.get_global_id(2);	      \
      { __VA_ARGS__ };				      \
     });	   			              \
    });

#define accelerator_barrier(dummy) theGridAccelerator->wait();

inline void *acceleratorAllocShared(size_t bytes){ return malloc_shared(bytes,*theGridAccelerator);};
inline void *acceleratorAllocDevice(size_t bytes){ return malloc_device(bytes,*theGridAccelerator);};
inline void acceleratorFreeShared(void *ptr){free(ptr,*theGridAccelerator);};
inline void acceleratorFreeDevice(void *ptr){free(ptr,*theGridAccelerator);};
inline void acceleratorCopyDeviceToDevice(void *from,void *to,size_t bytes)  { theGridAccelerator->memcpy(to,from,bytes); theGridAccelerator->wait();}
inline void acceleratorCopyToDevice(void *from,void *to,size_t bytes)  { theGridAccelerator->memcpy(to,from,bytes); theGridAccelerator->wait();}
inline void acceleratorCopyFromDevice(void *from,void *to,size_t bytes){ theGridAccelerator->memcpy(to,from,bytes); theGridAccelerator->wait();}
inline void acceleratorMemSet(void *base,int value,size_t bytes) { theGridAccelerator->memset(base,value,bytes); theGridAccelerator->wait();}
inline int  acceleratorIsCommunicable(void *ptr)
{
#if 0
  auto uvm = cl::sycl::usm::get_pointer_type(ptr, theGridAccelerator->get_context());
  if ( uvm = cl::sycl::usm::alloc::shared ) return 1;
  else return 0;
#endif
  return 1;
}

#endif

//////////////////////////////////////////////
// HIP acceleration
//////////////////////////////////////////////
#ifdef GRID_HIP
NAMESPACE_END(Grid);
#include <hip/hip_runtime.h>
NAMESPACE_BEGIN(Grid);

#ifdef __HIP_DEVICE_COMPILE__
#define GRID_SIMT
#endif

#define accelerator        __host__ __device__
#define accelerator_inline __host__ __device__ inline

/*These routines define mapping from thread grid to loop & vector lane indexing */
accelerator_inline int acceleratorSIMTlane(int Nsimd) {
#ifdef GRID_SIMT
  return hipThreadIdx_z; 
#else
  return 0;
#endif
} // HIP specific

#define accelerator_for2dNB( iter1, num1, iter2, num2, nsimd, ... )	\
  {									\
    typedef uint64_t Iterator;						\
    auto lambda = [=] accelerator					\
      (Iterator iter1,Iterator iter2,Iterator lane ) mutable {		\
      { __VA_ARGS__;}							\
    };									\
    int nt=acceleratorThreads();					\
    dim3 hip_threads(nt,1,nsimd);					\
    dim3 hip_blocks ((num1+nt-1)/nt,num2,1);				\
    hipLaunchKernelGGL(LambdaApply,hip_blocks,hip_threads,		\
		       0,0,						\
		       num1,num2,nsimd,lambda);				\
  }

template<typename lambda>  __global__
void LambdaApply(uint64_t numx, uint64_t numy, uint64_t numz, lambda Lambda)
{
  uint64_t x = hipThreadIdx_x + hipBlockDim_x*hipBlockIdx_x;
  uint64_t y = hipThreadIdx_y + hipBlockDim_y*hipBlockIdx_y;
  uint64_t z = hipThreadIdx_z ;//+ hipBlockDim_z*hipBlockIdx_z;
  if ( (x < numx) && (y<numy) && (z<numz) ) {
    Lambda(x,y,z);
  }
}

#define accelerator_barrier(dummy)				\
  {								\
    hipDeviceSynchronize();					\
    auto err = hipGetLastError();				\
    if ( err != hipSuccess ) {					\
      printf("After hipDeviceSynchronize() : HIP error %s \n", hipGetErrorString( err )); \
      puts(__FILE__);							\
      printf("Line %d\n",__LINE__);				\
      exit(0);							\
    }								\
  }

inline void *acceleratorAllocShared(size_t bytes)
{
  void *ptr=NULL;
  auto err = hipMallocManaged((void **)&ptr,bytes);
  if( err != hipSuccess ) {
    ptr = (void *) NULL;
    printf(" hipMallocManaged failed for %ld %s \n",bytes,hipGetErrorString(err));
  }
  return ptr;
};
inline int  acceleratorIsCommunicable(void *ptr){ return 1; }

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

inline void acceleratorFreeShared(void *ptr){ hipFree(ptr);};
inline void acceleratorFreeDevice(void *ptr){ hipFree(ptr);};
inline void acceleratorCopyToDevice(void *from,void *to,size_t bytes)  { hipMemcpy(to,from,bytes, hipMemcpyHostToDevice);}
inline void acceleratorCopyFromDevice(void *from,void *to,size_t bytes){ hipMemcpy(to,from,bytes, hipMemcpyDeviceToHost);}
inline void acceleratorCopyDeviceToDevice(void *from,void *to,size_t bytes)  { hipMemcpy(to,from,bytes, hipMemcpyDeviceToDevice);}
inline void acceleratorMemSet(void *base,int value,size_t bytes) { hipMemset(base,value,bytes);}

#endif

//////////////////////////////////////////////
// Common on all GPU targets
//////////////////////////////////////////////
#if defined(GRID_SYCL) || defined(GRID_CUDA) || defined(GRID_HIP)
#define accelerator_forNB( iter1, num1, nsimd, ... ) accelerator_for2dNB( iter1, num1, iter2, 1, nsimd, {__VA_ARGS__} );

#define accelerator_for( iter, num, nsimd, ... )		\
  accelerator_forNB(iter, num, nsimd, { __VA_ARGS__ } );	\
  accelerator_barrier(dummy);

#define accelerator_for2d(iter1, num1, iter2, num2, nsimd, ... )	\
  accelerator_for2dNB(iter1, num1, iter2, num2, nsimd, { __VA_ARGS__ } ); \
  accelerator_barrier(dummy);

#endif

//////////////////////////////////////////////
// CPU Target - No accelerator just thread instead
//////////////////////////////////////////////

#if ( (!defined(GRID_SYCL)) && (!defined(GRID_CUDA)) && (!defined(GRID_HIP)) )

#undef GRID_SIMT

#define accelerator 
#define accelerator_inline strong_inline
#define accelerator_for(iterator,num,nsimd, ... )   thread_for(iterator, num, { __VA_ARGS__ });
#define accelerator_forNB(iterator,num,nsimd, ... ) thread_for(iterator, num, { __VA_ARGS__ });
#define accelerator_barrier(dummy) 
#define accelerator_for2d(iter1, num1, iter2, num2, nsimd, ... ) thread_for2d(iter1,num1,iter2,num2,{ __VA_ARGS__ });

accelerator_inline int acceleratorSIMTlane(int Nsimd) { return 0; } // CUDA specific
inline void acceleratorCopyToDevice(void *from,void *to,size_t bytes)  { memcpy(to,from,bytes);}
inline void acceleratorCopyFromDevice(void *from,void *to,size_t bytes){ memcpy(to,from,bytes);}
inline void acceleratorCopyDeviceToDevice(void *from,void *to,size_t bytes)  { memcpy(to,from,bytes);}

inline int  acceleratorIsCommunicable(void *ptr){ return 1; }
inline void acceleratorMemSet(void *base,int value,size_t bytes) { memset(base,value,bytes);}
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

#ifdef HAVE_MM_MALLOC_H
inline void *acceleratorAllocCpu(size_t bytes){return _mm_malloc(bytes,GRID_ALLOC_ALIGN);};
inline void acceleratorFreeCpu  (void *ptr){_mm_free(ptr);};
#else
inline void *acceleratorAllocCpu(size_t bytes){return memalign(GRID_ALLOC_ALIGN,bytes);};
inline void acceleratorFreeCpu  (void *ptr){free(ptr);};
#endif



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
  //cl::sycl::detail::workGroupBarrier();
#endif
#ifdef GRID_HIP
  __syncthreads();
#endif
#endif
  return;
}
accelerator_inline void acceleratorSynchroniseAll(void) 
{
#ifdef GRID_SIMT
#ifdef GRID_CUDA
  __syncthreads();
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
accelerator_inline void acceleratorFence(void) 
{
#ifdef GRID_SIMT
#ifdef GRID_CUDA
  __threadfence();
#endif
#ifdef GRID_SYCL
  // FIXMEE
#endif
#ifdef GRID_HIP
  __threadfence();
#endif
#endif
  return;
}

NAMESPACE_END(Grid);
