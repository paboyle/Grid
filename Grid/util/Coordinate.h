/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/cartesian/Coordinate.h

    Copyright (C) 2018

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

////////////////////////////////////////////////////////////////////////////
// Provide a stack resident container for coordinates, or extents
////////////////////////////////////////////////////////////////////////////
NAMESPACE_BEGIN(Grid);
template<class _T,int MaxEntries>
class AcceleratorVector { 
public:

  typedef _T  value;
  typedef int size_type;
  typedef value & reference;
  typedef const value & const_reference;
  typedef value * pointer;
  typedef const value * const_pointer;

private:
  value _data[MaxEntries];
  size_type _size;
  
public:
  accelerator_inline reference       operator[](size_type __n)       { return _data[__n];}
  accelerator_inline const_reference operator[](size_type __n) const { return _data[__n];}
  accelerator_inline size_type size(void) const { return _size; };
  accelerator_inline void  clear(void) { resize(0);}
  accelerator_inline void  resize(size_type sz) {
#ifndef GRID_HIP
    assert(sz>=0);
    assert(sz<=MaxEntries);
#endif
    _size = sz;
  }
  accelerator_inline void  resize(size_type sz,const value &val) {
    resize(sz);
    for(int s=0;s<sz;s++) _data[s]=val;
  }
  accelerator_inline pointer begin(void)                   { return &_data[0]; } 
  accelerator_inline const_pointer begin(void) const       { return &_data[0]; } 
  accelerator_inline pointer end  (void)                   { return &_data[_size]; } 
  accelerator_inline const_pointer end  (void) const       { return &_data[_size]; } 
  accelerator_inline void push_back(const value &val)      { resize(_size+1); _data[_size-1] = val;}
  accelerator_inline AcceleratorVector()                   { resize(0); }
  accelerator_inline AcceleratorVector(size_type sz)           { resize(sz); }
  accelerator_inline AcceleratorVector(size_type sz,const value &val) { resize(sz,val); }
  AcceleratorVector(const std::vector<value> &copyme) { 
    resize(copyme.size());
    for(int s=0;s<_size;s++){
      _data[s] = copyme[s];
    }
  }
  std::vector<value> toVector(void) const { 
    auto clone =*this;
    std::vector<value> ret;
    std::copy(clone.begin(),clone.end(),std::back_inserter(ret));
    return ret;
  }
};

////////////////////////////////////////////////////////////////
// Coordinate class, maxdims = 8 for now.
////////////////////////////////////////////////////////////////
#define GRID_MAX_LATTICE_DIMENSION (8)
#define GRID_MAX_SIMD              (32)

static constexpr int MaxDims = GRID_MAX_LATTICE_DIMENSION;

typedef AcceleratorVector<int,MaxDims> Coordinate;

template<class T,int _ndim>
inline bool operator==(const AcceleratorVector<T,_ndim> &v,const AcceleratorVector<T,_ndim> &w)
{
  if (v.size()!=w.size()) return false;
  for(int i=0;i<v.size();i++) if ( v[i]!=w[i] ) return false;
  return true;
}
template<class T,int _ndim>
inline std::ostream & operator<<(std::ostream &os, const AcceleratorVector<T,_ndim> &v)
{
  os << "[";
  for(int s=0;s<v.size();s++) {
    os << v[s];
    if( s < (v.size()-1) ){
      os << " ";
    }
  }
  os << "]";
  return os;
}

NAMESPACE_END(Grid);
