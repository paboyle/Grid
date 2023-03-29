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
#define thread_for2d( i1, n1,i2,n2, ... )  \
  DO_PRAGMA(omp parallel for collapse(2))  \
  for ( uint64_t i1=0;i1<n1;i1++) {	   \
  for ( uint64_t i2=0;i2<n2;i2++) {	   \
  { __VA_ARGS__ } ;			   \
  }}
#define thread_foreach( i, container, ... )                 DO_PRAGMA(omp parallel for schedule(static)) for ( uint64_t i=container.begin();i<container.end();i++) { __VA_ARGS__ } ;
#define thread_for_in_region( i, num, ... )                 DO_PRAGMA(omp for schedule(static))          for ( uint64_t i=0;i<num;i++) { __VA_ARGS__ } ;
#define thread_for_collapse2( i, num, ... )                 DO_PRAGMA(omp parallel for collapse(2))      for ( uint64_t i=0;i<num;i++) { __VA_ARGS__ } ;
#define thread_for_collapse( N , i, num, ... )              DO_PRAGMA(omp parallel for collapse ( N ) )  for ( uint64_t i=0;i<num;i++) { __VA_ARGS__ } ;
#define thread_for_collapse_in_region( N , i, num, ... )    DO_PRAGMA(omp for collapse ( N ))            for ( uint64_t i=0;i<num;i++) { __VA_ARGS__ } ;
#define thread_region                                       DO_PRAGMA(omp parallel)
#define thread_critical                                     DO_PRAGMA(omp critical)

#ifdef GRID_OMP
inline void thread_bcopy(void *from, void *to,size_t bytes)
{
  uint64_t *ufrom = (uint64_t *)from;
  uint64_t *uto   = (uint64_t *)to;
  assert(bytes%8==0);
  uint64_t words=bytes/8;
  thread_for(w,words,{
      uto[w] = ufrom[w];
  });
}
#else
inline void thread_bcopy(void *from, void *to,size_t bytes)
{
  bcopy(from,to,bytes);
}
#endif
