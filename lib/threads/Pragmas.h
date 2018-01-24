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

#define COMMA_SAFE(...) __VA_ARGS__

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
#define thread_loop( range , body )           _Pragma("omp parallel for schedule(static)") for range { body ; };
#define thread_loop_in_region( range , body ) _Pragma("omp for schedule(static)")          for range  { body ; };
#define thread_loop_collapse( range , body )  _Pragma("omp parallel for collapse(2)")      for range  { body };
#define thread_region                         _Pragma("omp parallel")
#define thread_critical                       _Pragma("omp critical")
#else
#define thread_loop( range , body )            for range { body ; };
#define thread_loop_in_region( range , body )  for range { body ; };
#define thread_loop_collapse( range , body )   for range { body ; };
#define thread_region                           
#define thread_critical                         
#endif

//////////////////////////////////////////////////////////////////////////////////
// Deprecated primitives; to eridicate when done delete and fix errors.
//////////////////////////////////////////////////////////////////////////////////

#ifdef GRID_OMP
#define parallel_for  _Pragma("omp parallel for schedule(static)") for
#define parallel_for_internal  _Pragma("omp for schedule(static)") for
#define parallel_region    thread_region
#define parallel_for_nest2 for
#else
#define parallel_for           for
#define parallel_for_internal  for
#define parallel_region    
#define parallel_for_nest2     for
#endif

//////////////////////////////////////////////////////////////////////////////////
// Accelerator primitives; fall back to threading
//////////////////////////////////////////////////////////////////////////////////
#ifdef GRID_NVCC
#define accelerator __host__ __device__
#define accelerator_inline __host__ __device__ inline
// FIXME ; need to make this a CUDA kernel call
#define accelerator_loop( iterator, range, body )		\
  typedef decltype(range.begin()) Iterator;			\
  auto lambda = [&] (Iterator iterator) {			\
    body;							\
  };								\
  for(auto it=range.begin();it<range.end();it++){		\
    lambda(it);							\
  }
#define cpu_loop( iterator, range, body )   thread_loop( (auto iterator = range.begin();iterator<range.end();iterator++), { body });
#else
#define accelerator 
#define accelerator_inline strong_inline
#define accelerator_loop( iterator, range, body )			\
  thread_loop( (auto iterator = range.begin();iterator<range.end();iterator++), { body });
#define cpu_loop( iterator, range, body )				\
  thread_loop( (auto iterator = range.begin();iterator<range.end();iterator++), { body });
#endif
