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

#ifdef _OPENMP
#define GRID_OMP
#include <omp.h>
#endif

#ifdef __NVCC__
#define GRID_NVCC
#endif



//////////////////////////////////////////////////////////////////////////////////
// New primitives; explicit host thread calls, and accelerator data parallel calls
//////////////////////////////////////////////////////////////////////////////////
#ifdef GRID_OMP
#define thread_loop( range , ... )           _Pragma("omp parallel for schedule(static)") for range { __VA_ARGS__ ; };
#define thread_loop_in_region( range , ... ) _Pragma("omp for schedule(static)")          for range  { __VA_ARGS__ ; };
#define thread_loop_collapse( n, range , ... )  _Pragma("omp parallel for collapse(" #n ")")      for range  { __VA_ARGS__ };
#define thread_region                         _Pragma("omp parallel")
#define thread_critical                       _Pragma("omp critical")
#define thread_num(a) omp_get_thread_num()
#define thread_max(a) omp_get_max_threads()
#else
#define thread_loop( range , ... )            for range { __VA_ARGS__ ; };
#define thread_loop_in_region( range , ... )  for range { __VA_ARGS__ ; };
#define thread_loop_collapse( n, range , ... )   for range { __VA_ARGS__ ; };
#define thread_region                           
#define thread_critical                         
#define thread_num(a) (0)
#define thread_max(a) (1)
#endif


//////////////////////////////////////////////////////////////////////////////////
// Accelerator primitives; fall back to threading
//////////////////////////////////////////////////////////////////////////////////
#ifdef GRID_NVCC
#define accelerator __host__ __device__
#define accelerator_inline __host__ __device__ inline
// FIXME ; need to make this a CUDA kernel call
#define accelerator_loop( iterator, range, ... )		\
  typedef decltype(range.begin()) Iterator;			\
  auto lambda = [&] (Iterator iterator) {			\
    __VA_ARGS__;							\
  };								\
  for(auto it=range.begin();it<range.end();it++){		\
    lambda(it);							\
  }
#define cpu_loop( iterator, range, ... )   thread_loop( (auto iterator = range.begin();iterator<range.end();iterator++), { __VA_ARGS__ });
#else
#define accelerator 
#define accelerator_inline strong_inline
#define accelerator_loop( iterator, range, ... )			\
  thread_loop( (auto iterator = range.begin();iterator<range.end();iterator++), { __VA_ARGS__ });
#define cpu_loop( iterator, range, ... )				\
  thread_loop( (auto iterator = range.begin();iterator<range.end();iterator++), { __VA_ARGS__ });
#endif
