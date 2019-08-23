/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/Threads.h

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

#ifndef MAX
#define MAX(x,y) ((x)>(y)?(x):(y))
#define MIN(x,y) ((x)>(y)?(y):(x))
#endif

#define strong_inline     __attribute__((always_inline)) inline
#define UNROLL  _Pragma("unroll")

//////////////////////////////////////////////////////////////////////////////////
// New primitives; explicit host thread calls, and accelerator data parallel calls
//////////////////////////////////////////////////////////////////////////////////

#ifdef _OPENMP
#define GRID_OMP
#include <omp.h>
#endif

#ifdef GRID_OMP
#define DO_PRAGMA_(x) _Pragma (#x)
#define DO_PRAGMA(x) DO_PRAGMA_(x)
#define thread_num(a) omp_get_thread_num()
#define thread_max(a) omp_get_max_threads()
#else 
#define DO_PRAGMA_(x) 
#define DO_PRAGMA(x) 
#define thread_num(a) (0)
#define thread_max(a) (1)
#endif

#define thread_for( i, num, ... )                           DO_PRAGMA(omp parallel for schedule(static)) for ( uint64_t i=0;i<num;i++) { __VA_ARGS__ } ;
#define thread_foreach( i, container, ... )                 DO_PRAGMA(omp parallel for schedule(static)) for ( uint64_t i=container.begin();i<container.end();i++) { __VA_ARGS__ } ;
#define thread_for_in_region( i, num, ... )                 DO_PRAGMA(omp for schedule(static))          for ( uint64_t i=0;i<num;i++) { __VA_ARGS__ } ;
#define thread_for_collapse2( i, num, ... )                 DO_PRAGMA(omp parallel for collapse(2))      for ( uint64_t i=0;i<num;i++) { __VA_ARGS__ } ;
#define thread_for_collapse( N , i, num, ... )              DO_PRAGMA(omp parallel for collapse ( N ) )  for ( uint64_t i=0;i<num;i++) { __VA_ARGS__ } ;
#define thread_for_collapse_in_region( N , i, num, ... )    DO_PRAGMA(omp for collapse ( N ))            for ( uint64_t i=0;i<num;i++) { __VA_ARGS__ } ;
#define thread_region                                       DO_PRAGMA(omp parallel)
#define thread_critical                                     DO_PRAGMA(omp critical)


//////////////////////////////////////////////////////////////////////////////////
// Accelerator primitives; fall back to threading
//////////////////////////////////////////////////////////////////////////////////
#ifdef __NVCC__
#define GRID_NVCC
#endif

#ifdef GRID_NVCC

extern uint32_t gpu_threads;

#define accelerator        __host__ __device__
#define accelerator_inline __host__ __device__ inline

template<typename lambda>  __global__
void LambdaApplySIMT(uint64_t Isites, uint64_t Osites, lambda Lambda)
{
  uint64_t isite = threadIdx.y;
  uint64_t osite = threadIdx.x+blockDim.x*blockIdx.x;
  if ( (osite <Osites) && (isite<Isites) ) {
    Lambda(isite,osite);
  }
}

/////////////////////////////////////////////////////////////////
// Internal only really... but need to call when 
/////////////////////////////////////////////////////////////////
#define accelerator_barrier(dummy)				\
  {								\
    cudaDeviceSynchronize();					\
    cudaError err = cudaGetLastError();				\
    if ( cudaSuccess != err ) {					\
      printf("Cuda error %s \n", cudaGetErrorString( err )); \
      puts(__FILE__); \
      printf("Line %d\n",__LINE__);					\
      exit(0);							\
    }								\
  }

// Copy the for_each_n style ; Non-blocking variant
#define accelerator_forNB( iterator, num, nsimd, ... )			\
  {									\
    typedef uint64_t Iterator;						\
    auto lambda = [=] accelerator (Iterator lane,Iterator iterator) mutable { \
      __VA_ARGS__;							\
    };									\
    dim3 cu_threads(gpu_threads,nsimd);					\
    dim3 cu_blocks ((num+gpu_threads-1)/gpu_threads);			\
    LambdaApplySIMT<<<cu_blocks,cu_threads>>>(nsimd,num,lambda);	\
  }

// Copy the for_each_n style ; Non-blocking variant (default
#define accelerator_for( iterator, num, nsimd, ... )		\
  accelerator_forNB(iterator, num, nsimd, { __VA_ARGS__ } );	\
  accelerator_barrier(dummy);

#else

#define accelerator 
#define accelerator_inline strong_inline
#define accelerator_for(iterator,num,nsimd, ... )   thread_for(iterator, num, { __VA_ARGS__ });
#define accelerator_forNB(iterator,num,nsimd, ... ) thread_for(iterator, num, { __VA_ARGS__ });
#define accelerator_barrier(dummy) 

#endif
