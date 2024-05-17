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
#pragma once

NAMESPACE_BEGIN(Grid);

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
    _Tp *ptr = (_Tp*) MemoryManager::CpuAllocate(bytes);
    assert( ( (_Tp*)ptr != (_Tp *)NULL ) );
    return ptr;
  }

  void deallocate(pointer __p, size_type __n) 
  { 
    size_type bytes = __n * sizeof(_Tp);
    profilerFree(bytes);
    MemoryManager::CpuFree((void *)__p,bytes);
  }

  // FIXME: hack for the copy constructor: it must be avoided to avoid single thread loop
  void construct(pointer __p, const _Tp& __val) { assert(0);};
  void construct(pointer __p) { };
  void destroy(pointer __p) { };
};
template<typename _Tp>  inline bool operator==(const alignedAllocator<_Tp>&, const alignedAllocator<_Tp>&){ return true; }
template<typename _Tp>  inline bool operator!=(const alignedAllocator<_Tp>&, const alignedAllocator<_Tp>&){ return false; }

//////////////////////////////////////////////////////////////////////////////////////
// Unified virtual memory
//////////////////////////////////////////////////////////////////////////////////////
template<typename _Tp>
class uvmAllocator {
public: 
  typedef std::size_t     size_type;
  typedef std::ptrdiff_t  difference_type;
  typedef _Tp*       pointer;
  typedef const _Tp* const_pointer;
  typedef _Tp&       reference;
  typedef const _Tp& const_reference;
  typedef _Tp        value_type;

  template<typename _Tp1>  struct rebind { typedef uvmAllocator<_Tp1> other; };
  uvmAllocator() throw() { }
  uvmAllocator(const uvmAllocator&) throw() { }
  template<typename _Tp1> uvmAllocator(const uvmAllocator<_Tp1>&) throw() { }
  ~uvmAllocator() throw() { }
  pointer       address(reference __x)       const { return &__x; }
  size_type  max_size() const throw() { return size_t(-1) / sizeof(_Tp); }

  pointer allocate(size_type __n, const void* _p= 0)
  { 
    size_type bytes = __n*sizeof(_Tp);
    profilerAllocate(bytes);
    _Tp *ptr = (_Tp*) MemoryManager::SharedAllocate(bytes);
    assert( ( (_Tp*)ptr != (_Tp *)NULL ) );
    return ptr;
  }

  void deallocate(pointer __p, size_type __n) 
  { 
    size_type bytes = __n * sizeof(_Tp);
    profilerFree(bytes);
    MemoryManager::SharedFree((void *)__p,bytes);
  }

  void construct(pointer __p, const _Tp& __val) { new((void *)__p) _Tp(__val); };
  void construct(pointer __p) { };
  void destroy(pointer __p) { };
};
template<typename _Tp>  inline bool operator==(const uvmAllocator<_Tp>&, const uvmAllocator<_Tp>&){ return true; }
template<typename _Tp>  inline bool operator!=(const uvmAllocator<_Tp>&, const uvmAllocator<_Tp>&){ return false; }

////////////////////////////////////////////////////////////////////////////////
// Device memory
////////////////////////////////////////////////////////////////////////////////
template<typename _Tp>
class devAllocator {
public: 
  typedef std::size_t     size_type;
  typedef std::ptrdiff_t  difference_type;
  typedef _Tp*       pointer;
  typedef const _Tp* const_pointer;
  typedef _Tp&       reference;
  typedef const _Tp& const_reference;
  typedef _Tp        value_type;

  template<typename _Tp1>  struct rebind { typedef devAllocator<_Tp1> other; };
  devAllocator() throw() { }
  devAllocator(const devAllocator&) throw() { }
  template<typename _Tp1> devAllocator(const devAllocator<_Tp1>&) throw() { }
  ~devAllocator() throw() { }
  pointer       address(reference __x)       const { return &__x; }
  size_type  max_size() const throw() { return size_t(-1) / sizeof(_Tp); }

  pointer allocate(size_type __n, const void* _p= 0)
  { 
    size_type bytes = __n*sizeof(_Tp);
    profilerAllocate(bytes);
    _Tp *ptr = (_Tp*) MemoryManager::AcceleratorAllocate(bytes);
    assert( ( (_Tp*)ptr != (_Tp *)NULL ) );
    return ptr;
  }

  void deallocate(pointer __p, size_type __n) 
  { 
    size_type bytes = __n * sizeof(_Tp);
    profilerFree(bytes);
    MemoryManager::AcceleratorFree((void *)__p,bytes);
  }
  void construct(pointer __p, const _Tp& __val) { };
  void construct(pointer __p) { };
  void destroy(pointer __p) { };
};
template<typename _Tp>  inline bool operator==(const devAllocator<_Tp>&, const devAllocator<_Tp>&){ return true; }
template<typename _Tp>  inline bool operator!=(const devAllocator<_Tp>&, const devAllocator<_Tp>&){ return false; }

////////////////////////////////////////////////////////////////////////////////
// Template typedefs
////////////////////////////////////////////////////////////////////////////////
#ifdef ACCELERATOR_CSHIFT
// Cshift on device
template<class T> using cshiftAllocator = devAllocator<T>;
#else
// Cshift on host
template<class T> using cshiftAllocator = std::allocator<T>;
#endif

template<class T> using Vector        = std::vector<T,uvmAllocator<T> >;           
template<class T> using stencilVector = std::vector<T,alignedAllocator<T> >;           
template<class T> using commVector    = std::vector<T,devAllocator<T> >;
template<class T> using deviceVector  = std::vector<T,devAllocator<T> >;
template<class T> using cshiftVector  = std::vector<T,cshiftAllocator<T> >;

/*
template<class T> class vecView
{
 protected:
  T * data;
  uint64_t size;
  ViewMode mode;
  void * cpu_ptr;
 public:
  accelerator_inline T & operator[](size_t i) const { return this->data[i]; };
  vecView(std::vector<T> &refer_to_me,ViewMode _mode)
  {
    cpu_ptr = &refer_to_me[0];
    size = refer_to_me.size();
    mode = _mode;
    data =(T *) MemoryManager::ViewOpen(cpu_ptr,
					size*sizeof(T),
					mode,
					AdviseDefault);
  }
  void ViewClose(void)
  { // Inform the manager
    MemoryManager::ViewClose(this->cpu_ptr,this->mode);    
  }
};

template<class T> vecView<T> VectorView(std::vector<T> &vec,ViewMode _mode)
{
  vecView<T> ret(vec,_mode); // does the open
  return ret;                // must be closed
}

// Little autoscope assister
template<class View> 
class VectorViewCloser
{
  View v;  // Take a copy of view and call view close when I go out of scope automatically
 public:
  VectorViewCloser(View &_v) : v(_v) {};
  ~VectorViewCloser() { auto ptr = v.cpu_ptr; v.ViewClose();  MemoryManager::NotifyDeletion(ptr);}
};

#define autoVecView(v_v,v,mode)					\
  auto v_v = VectorView(v,mode);				\
  ViewCloser<decltype(v_v)> _autoView##v_v(v_v);
*/

NAMESPACE_END(Grid);


