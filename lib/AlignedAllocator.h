    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/AlignedAllocator.h

    Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>

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
#ifndef GRID_ALIGNED_ALLOCATOR_H
#define GRID_ALIGNED_ALLOCATOR_H

#ifdef HAVE_MALLOC_MALLOC_H
#include <malloc/malloc.h>
#endif
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif

#ifdef HAVE_MM_MALLOC_H
#include <mm_malloc.h>
#endif

#ifdef GRID_COMMS_SHMEM
extern "C" { 
#include <mpp/shmem.h>
extern void * shmem_align(size_t, size_t);
extern void  shmem_free(void *);
}
#endif

namespace Grid {

////////////////////////////////////////////////////////////////////
// A lattice of something, but assume the something is SIMDized.
////////////////////////////////////////////////////////////////////
template<typename _Tp>
class alignedAllocator {
public: 
  typedef std::size_t     size_type;
  typedef std::ptrdiff_t  difference_type;
  typedef _Tp*       pointer;
  typedef const _Tp* const_pointer;
  typedef _Tp&       reference;
  typedef const _Tp& const_reference;
  typedef _Tp        value_type;

  template<typename _Tp1>  struct rebind { typedef alignedAllocator<_Tp1> other; };

  alignedAllocator() throw() { }

  alignedAllocator(const alignedAllocator&) throw() { }

  template<typename _Tp1> alignedAllocator(const alignedAllocator<_Tp1>&) throw() { }

  ~alignedAllocator() throw() { }

  pointer       address(reference __x)       const { return &__x; }
  //  const_pointer address(const_reference __x) const { return &__x; }

  size_type  max_size() const throw() { return size_t(-1) / sizeof(_Tp); }

  pointer allocate(size_type __n, const void* _p= 0)
  { 
#ifdef GRID_COMMS_SHMEM

    _Tp *ptr = (_Tp *) shmem_align(__n*sizeof(_Tp),64);


#define PARANOID_SYMMETRIC_HEAP
#ifdef PARANOID_SYMMETRIC_HEAP
    static void * bcast;
    static long  psync[_SHMEM_REDUCE_SYNC_SIZE];

    bcast = (void *) ptr;
    shmem_broadcast32((void *)&bcast,(void *)&bcast,sizeof(void *)/4,0,0,0,shmem_n_pes(),psync);

    if ( bcast != ptr ) {
      std::printf("inconsistent alloc pe %d %lx %lx \n",shmem_my_pe(),bcast,ptr);std::fflush(stdout);
      BACKTRACEFILE();
      exit(0);
    }

    assert( bcast == (void *) ptr);

#endif 
#else

#ifdef HAVE_MM_MALLOC_H
    _Tp * ptr = (_Tp *) _mm_malloc(__n*sizeof(_Tp),128);
#else
    _Tp * ptr = (_Tp *) memalign(128,__n*sizeof(_Tp));
#endif

#endif
    _Tp tmp;
#undef FIRST_TOUCH_OPTIMISE
#ifdef FIRST_TOUCH_OPTIMISE
#pragma omp parallel for 
  for(int i=0;i<__n;i++){
    ptr[i]=tmp;
  }
#endif 
    return ptr;
  }

  void deallocate(pointer __p, size_type) { 
#ifdef GRID_COMMS_SHMEM
    shmem_free((void *)__p);
#else
#ifdef HAVE_MM_MALLOC_H
    _mm_free((void *)__p); 
#else
    free((void *)__p);
#endif
#endif
  }
  void construct(pointer __p, const _Tp& __val) { };
  void construct(pointer __p) { };

  void destroy(pointer __p) { };
};

template<typename _Tp>  inline bool
operator==(const alignedAllocator<_Tp>&, const alignedAllocator<_Tp>&){ return true; }

template<typename _Tp>  inline bool
operator!=(const alignedAllocator<_Tp>&, const alignedAllocator<_Tp>&){ return false; }
    
}; // namespace Grid
#endif
