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

#define POINTER_CACHE
#define GRID_ALLOC_ALIGN (2*1024*1024)

NAMESPACE_BEGIN(Grid);

// Move control to configure.ac and Config.h?
#ifdef POINTER_CACHE
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
#endif  

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
#define profilerDebugPrint						\
  if (MemoryProfiler::stats)						\
    {									\
      auto s = MemoryProfiler::stats;					\
      std::cout << GridLogDebug << "[Memory debug] Stats " << MemoryProfiler::stats << std::endl; \
      std::cout << GridLogDebug << "[Memory debug] total  : " << memString(s->totalAllocated) \
		<< std::endl;						\
      std::cout << GridLogDebug << "[Memory debug] max    : " << memString(s->maxAllocated) \
		<< std::endl;						\
      std::cout << GridLogDebug << "[Memory debug] current: " << memString(s->currentlyAllocated) \
		<< std::endl;						\
      std::cout << GridLogDebug << "[Memory debug] freed  : " << memString(s->totalFreed) \
		<< std::endl;						\
    }

#define profilerAllocate(bytes)						\
  if (MemoryProfiler::stats)						\
    {									\
      auto s = MemoryProfiler::stats;					\
      s->totalAllocated     += (bytes);					\
      s->currentlyAllocated += (bytes);					\
      s->maxAllocated        = std::max(s->maxAllocated, s->currentlyAllocated); \
    }									\
  if (MemoryProfiler::debug)						\
    {									\
      std::cout << GridLogDebug << "[Memory debug] allocating " << memString(bytes) << std::endl; \
      profilerDebugPrint;						\
    }

#define profilerFree(bytes)						\
  if (MemoryProfiler::stats)						\
    {									\
      auto s = MemoryProfiler::stats;					\
      s->totalFreed         += (bytes);					\
      s->currentlyAllocated -= (bytes);					\
    }									\
  if (MemoryProfiler::debug)						\
    {									\
      std::cout << GridLogDebug << "[Memory debug] freeing " << memString(bytes) << std::endl; \
      profilerDebugPrint;						\
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


#ifdef POINTER_CACHE
    _Tp *ptr = (_Tp *) PointerCache::Lookup(bytes);
#else
    pointer ptr = nullptr;
#endif

#ifdef GRID_NVCC
    ////////////////////////////////////
    // Unified (managed) memory
    ////////////////////////////////////
    if ( ptr == (_Tp *) NULL ) {
      auto err = cudaMallocManaged((void **)&ptr,bytes);
      if( err != cudaSuccess ) {
	ptr = (_Tp *) NULL;
	std::cerr << " cudaMallocManaged failed for " << bytes<<" bytes " <<cudaGetErrorString(err)<< std::endl;
	assert(0);
      }
    }
#else 
    //////////////////////////////////////////////////////////////////////////////////////////
    // 2MB align; could make option probably doesn't need configurability
    //////////////////////////////////////////////////////////////////////////////////////////
  #ifdef HAVE_MM_MALLOC_H
    if ( ptr == (_Tp *) NULL ) ptr = (_Tp *) _mm_malloc(bytes,GRID_ALLOC_ALIGN);
  #else
    if ( ptr == (_Tp *) NULL ) ptr = (_Tp *) memalign(GRID_ALLOC_ALIGN,bytes);
  #endif
#endif
    assert( ptr != (_Tp *)NULL);

    /////////////////////////////////////////
    // First touch optimise in threaded loop
    /////////////////////////////////////////
    uint8_t *cp = (uint8_t *)ptr;
    thread_loop( (size_type n=0;n<bytes;n+=4096) , {
      cp[n]=0;
    });
    return ptr;
  }

  void deallocate(pointer __p, size_type __n) { 
    size_type bytes = __n * sizeof(_Tp);

    profilerFree(bytes);

#ifdef POINTER_CACHE
    pointer __freeme = (pointer)PointerCache::Insert((void *)__p,bytes);
#else 
    pointer __freeme = __p;
#endif

#ifdef GRID_NVCC
    if ( __freeme ) cudaFree((void *)__freeme);
#else 
  #ifdef HAVE_MM_MALLOC_H
    if ( __freeme ) _mm_free((void *)__freeme); 
  #else
    if ( __freeme ) free((void *)__freeme);
  #endif
#endif
  }
  void construct(pointer __p, const _Tp& __val) { };
  void construct(pointer __p) { };
  void destroy(pointer __p) { };
};
template<typename _Tp>  inline bool operator==(const alignedAllocator<_Tp>&, const alignedAllocator<_Tp>&){ return true; }
template<typename _Tp>  inline bool operator!=(const alignedAllocator<_Tp>&, const alignedAllocator<_Tp>&){ return false; }

////////////////////////////////////////////////////////////////////////////////
// Template typedefs
////////////////////////////////////////////////////////////////////////////////
template<class T> using commAllocator = alignedAllocator<T>;
template<class T> using Vector     = std::vector<T,alignedAllocator<T> >;           
template<class T> using commVector = std::vector<T,alignedAllocator<T> >;
template<class T> using Matrix     = std::vector<std::vector<T,alignedAllocator<T> > >;

NAMESPACE_END(Grid);

#endif
