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

namespace Grid {

  class PointerCache {
  private:

    static const int Ncache=8;
    static int victim;

    typedef struct { 
      void *address;
      size_t bytes;
      int valid;
    } PointerCacheEntry;
    
    static PointerCacheEntry Entries[Ncache];

  public:


    static void *Insert(void *ptr,size_t bytes) ;
    static void *Lookup(size_t bytes) ;

  };
  
  std::string sizeString(size_t bytes);

  struct MemoryStats
  {
    size_t totalAllocated{0}, maxAllocated{0}, 
           currentlyAllocated{0}, totalFreed{0};
  };
    
  class MemoryProfiler
  {
  public:
    static MemoryStats *stats;
    static bool        debug;
  };

  #define memString(bytes) std::to_string(bytes) + " (" + sizeString(bytes) + ")"
  #define profilerDebugPrint \
  if (MemoryProfiler::stats)\
  {\
    auto s = MemoryProfiler::stats;\
    std::cout << GridLogDebug << "[Memory debug] Stats " << MemoryProfiler::stats << std::endl;\
    std::cout << GridLogDebug << "[Memory debug] total  : " << memString(s->totalAllocated) \
              << std::endl;\
    std::cout << GridLogDebug << "[Memory debug] max    : " << memString(s->maxAllocated) \
              << std::endl;\
    std::cout << GridLogDebug << "[Memory debug] current: " << memString(s->currentlyAllocated) \
              << std::endl;\
    std::cout << GridLogDebug << "[Memory debug] freed  : " << memString(s->totalFreed) \
              << std::endl;\
  }

  #define profilerAllocate(bytes)\
  if (MemoryProfiler::stats)\
  {\
    auto s = MemoryProfiler::stats;\
    s->totalAllocated     += (bytes);\
    s->currentlyAllocated += (bytes);\
    s->maxAllocated        = std::max(s->maxAllocated, s->currentlyAllocated);\
  }\
  if (MemoryProfiler::debug)\
  {\
    std::cout << GridLogDebug << "[Memory debug] allocating " << memString(bytes) << std::endl;\
    profilerDebugPrint;\
  }

  #define profilerFree(bytes)\
  if (MemoryProfiler::stats)\
  {\
    auto s = MemoryProfiler::stats;\
    s->totalFreed         += (bytes);\
    s->currentlyAllocated -= (bytes);\
  }\
  if (MemoryProfiler::debug)\
  {\
    std::cout << GridLogDebug << "[Memory debug] freeing " << memString(bytes) << std::endl;\
    profilerDebugPrint;\
  }

  void check_huge_pages(void *Buf,uint64_t BYTES);

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
  size_type  max_size() const throw() { return size_t(-1) / sizeof(_Tp); }

  pointer allocate(size_type __n, const void* _p= 0)
  { 
    size_type bytes = __n*sizeof(_Tp);
    profilerAllocate(bytes);

    _Tp *ptr = (_Tp *) PointerCache::Lookup(bytes);
    //    if ( ptr != NULL ) 
    //      std::cout << "alignedAllocator "<<__n << " cache hit "<< std::hex << ptr <<std::dec <<std::endl;

    //////////////////
    // Hack 2MB align; could make option probably doesn't need configurability
    //////////////////
//define GRID_ALLOC_ALIGN (128)
#define GRID_ALLOC_ALIGN (2*1024*1024)
#ifdef HAVE_MM_MALLOC_H
    if ( ptr == (_Tp *) NULL ) ptr = (_Tp *) _mm_malloc(bytes,GRID_ALLOC_ALIGN);
#else
    if ( ptr == (_Tp *) NULL ) ptr = (_Tp *) memalign(GRID_ALLOC_ALIGN,bytes);
#endif
    //    std::cout << "alignedAllocator " << std::hex << ptr <<std::dec <<std::endl;
    // First touch optimise in threaded loop
    uint8_t *cp = (uint8_t *)ptr;
#ifdef GRID_OMP
#pragma omp parallel for
#endif
    for(size_type n=0;n<bytes;n+=4096){
      cp[n]=0;
    }
    return ptr;
  }

  void deallocate(pointer __p, size_type __n) { 
    size_type bytes = __n * sizeof(_Tp);

    profilerFree(bytes);

    pointer __freeme = (pointer)PointerCache::Insert((void *)__p,bytes);

#ifdef HAVE_MM_MALLOC_H
    if ( __freeme ) _mm_free((void *)__freeme); 
#else
    if ( __freeme ) free((void *)__freeme);
#endif
  }
  void construct(pointer __p, const _Tp& __val) { };
  void construct(pointer __p) { };
  void destroy(pointer __p) { };
};
template<typename _Tp>  inline bool operator==(const alignedAllocator<_Tp>&, const alignedAllocator<_Tp>&){ return true; }
template<typename _Tp>  inline bool operator!=(const alignedAllocator<_Tp>&, const alignedAllocator<_Tp>&){ return false; }

//////////////////////////////////////////////////////////////////////////////////////////
// MPI3 : comms must use shm region
// SHMEM: comms must use symmetric heap
//////////////////////////////////////////////////////////////////////////////////////////
#ifdef GRID_COMMS_SHMEM
extern "C" { 
#include <mpp/shmem.h>
extern void * shmem_align(size_t, size_t);
extern void  shmem_free(void *);
}
#define PARANOID_SYMMETRIC_HEAP
#endif

template<typename _Tp>
class commAllocator {
public: 
  typedef std::size_t     size_type;
  typedef std::ptrdiff_t  difference_type;
  typedef _Tp*       pointer;
  typedef const _Tp* const_pointer;
  typedef _Tp&       reference;
  typedef const _Tp& const_reference;
  typedef _Tp        value_type;

  template<typename _Tp1>  struct rebind { typedef commAllocator<_Tp1> other; };
  commAllocator() throw() { }
  commAllocator(const commAllocator&) throw() { }
  template<typename _Tp1> commAllocator(const commAllocator<_Tp1>&) throw() { }
  ~commAllocator() throw() { }
  pointer       address(reference __x)       const { return &__x; }
  size_type  max_size() const throw() { return size_t(-1) / sizeof(_Tp); }

#ifdef GRID_COMMS_SHMEM
  pointer allocate(size_type __n, const void* _p= 0)
  {
    size_type bytes = __n*sizeof(_Tp);

    profilerAllocate(bytes);
#ifdef CRAY
    _Tp *ptr = (_Tp *) shmem_align(bytes,64);
#else
    _Tp *ptr = (_Tp *) shmem_align(64,bytes);
#endif
#ifdef PARANOID_SYMMETRIC_HEAP
    static void * bcast;
    static long  psync[_SHMEM_REDUCE_SYNC_SIZE];

    bcast = (void *) ptr;
    shmem_broadcast32((void *)&bcast,(void *)&bcast,sizeof(void *)/4,0,0,0,shmem_n_pes(),psync);

    if ( bcast != ptr ) {
      std::printf("inconsistent alloc pe %d %lx %lx \n",shmem_my_pe(),bcast,ptr);std::fflush(stdout);
      //      BACKTRACEFILE();
      exit(0);
    }
    assert( bcast == (void *) ptr);
#endif 
    return ptr;
  }
  void deallocate(pointer __p, size_type __n) { 
    size_type bytes = __n*sizeof(_Tp);

    profilerFree(bytes);
    shmem_free((void *)__p);
  }
#else
  pointer allocate(size_type __n, const void* _p= 0) 
  {
    size_type bytes = __n*sizeof(_Tp);
    
    profilerAllocate(bytes);
#ifdef HAVE_MM_MALLOC_H
    _Tp * ptr = (_Tp *) _mm_malloc(bytes, GRID_ALLOC_ALIGN);
#else
    _Tp * ptr = (_Tp *) memalign(GRID_ALLOC_ALIGN, bytes);
#endif
    uint8_t *cp = (uint8_t *)ptr;
    if ( ptr ) { 
    // One touch per 4k page, static OMP loop to catch same loop order
#ifdef GRID_OMP
#pragma omp parallel for schedule(static)
#endif
      for(size_type n=0;n<bytes;n+=4096){
	cp[n]=0;
      }
    }
    return ptr;
  }
  void deallocate(pointer __p, size_type __n) {
    size_type bytes = __n*sizeof(_Tp);

    profilerFree(bytes);
#ifdef HAVE_MM_MALLOC_H
    _mm_free((void *)__p); 
#else
    free((void *)__p);
#endif
  }
#endif
  void construct(pointer __p, const _Tp& __val) { };
  void construct(pointer __p) { };
  void destroy(pointer __p) { };
};
template<typename _Tp>  inline bool operator==(const commAllocator<_Tp>&, const commAllocator<_Tp>&){ return true; }
template<typename _Tp>  inline bool operator!=(const commAllocator<_Tp>&, const commAllocator<_Tp>&){ return false; }

////////////////////////////////////////////////////////////////////////////////
// Template typedefs
////////////////////////////////////////////////////////////////////////////////
template<class T> using Vector     = std::vector<T,alignedAllocator<T> >;           
template<class T> using commVector = std::vector<T,commAllocator<T> >;              
template<class T> using Matrix     = std::vector<std::vector<T,alignedAllocator<T> > >;
    
}; // namespace Grid
#endif
