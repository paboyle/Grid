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
#define thread_loop_collapse2( range , ... )  _Pragma("omp parallel for collapse(2)")     for range  { __VA_ARGS__ };
#define thread_loop_collapse3( range , ... )  _Pragma("omp parallel for collapse(3)")     for range  { __VA_ARGS__ };
#define thread_loop_collapse4( range , ... )  _Pragma("omp parallel for collapse(4)")     for range  { __VA_ARGS__ };
#define thread_region                         _Pragma("omp parallel")
#define thread_critical                       _Pragma("omp critical")
#define thread_num(a) omp_get_thread_num()
#define thread_max(a) omp_get_max_threads()
#else
#define thread_loop( range , ... )            for range { __VA_ARGS__ ; };
#define thread_loop_in_region( range , ... )  for range { __VA_ARGS__ ; };
#define thread_loop_collapse2( range , ... )  for range { __VA_ARGS__ ; };
#define thread_loop_collapse3( range , ... )  for range { __VA_ARGS__ ; };
#define thread_loop_collapse4( range , ... )  for range { __VA_ARGS__ ; };
#define thread_region                           
#define thread_critical                         
#define thread_num(a) (0)
#define thread_max(a) (1)
#endif


//////////////////////////////////////////////////////////////////////////////////
// Accelerator primitives; fall back to threading
//////////////////////////////////////////////////////////////////////////////////
#ifdef GRID_NVCC

template<typename lambda>  __global__
void LambdaApply(int Num, lambda &Lambda)
{
  int ss = blockIdx.x;
  if ( ss < Num ) {
    Lambda(ss);
  }
}

#define accelerator_exec( ... )					\
  auto lambda = [=] accelerator (void) mutable {		\
    __VA_ARGS__;						\
  };								\
  lambda();							

#define accelerator __host__ __device__
#define accelerator_inline __host__ __device__ inline

#if 0 
#define accelerator_loop( iterator, range, ... )		\
  typedef decltype(range.begin()) Iterator;			\
  auto lambda = [=] accelerator (Iterator iterator) mutable {	\
    __VA_ARGS__;						\
  };								\
  for(auto it=range.begin();it<range.end();it++){		\
    lambda(it);							\
  }
#else
#define accelerator_loop( iterator, range, ... )		\
  typedef decltype(range.begin()) Iterator;			\
  auto lambda = [=] accelerator (Iterator iterator) mutable {	\
    __VA_ARGS__;						\
  };								\
  Iterator num = range.end();					\
  LambdaApply<<<num,1>>>(num,lambda);
#endif

#define cpu_loop( iterator, range, ... )   thread_loop( (auto iterator = range.begin();iterator<range.end();iterator++), { __VA_ARGS__ });

#define NVCC_DECLARE_VECTOR_ACCESSOR(Type) \
template <> accelerator_inline typename Grid::Vector<Type>::reference       Grid::Vector<Type>::operator[](Grid::Vector<Type>::size_type __n)       { return this->__begin_[__n]; } \
template <> accelerator_inline typename Grid::Vector<Type>::const_reference Grid::Vector<Type>::operator[](Grid::Vector<Type>::size_type __n) const { return this->__begin_[__n]; }

#define NVCC_DECLARE_STD_VECTOR_ACCESSOR(Type) \
template <> accelerator_inline typename std::vector<Type>::reference       std::vector<Type>::operator[](std::vector<Type>::size_type __n)       { return this->__begin_[__n]; } \
template <> accelerator_inline typename std::vector<Type>::const_reference std::vector<Type>::operator[](std::vector<Type>::size_type __n) const { return this->__begin_[__n]; }

#else
#define accelerator 
#define accelerator_inline strong_inline
#define accelerator_loop( iterator, range, ... )			\
  thread_loop( (auto iterator = range.begin();iterator<range.end();iterator++), { __VA_ARGS__ });
#define cpu_loop( iterator, range, ... )				\
  thread_loop( (auto iterator = range.begin();iterator<range.end();iterator++), { __VA_ARGS__ });

#define NVCC_DECLARE_STD_VECTOR_ACCESSOR(Type) 
#define NVCC_DECLARE_VECTOR_ACCESSOR(Type) 

#endif

NVCC_DECLARE_STD_VECTOR_ACCESSOR(int);
