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
    assert(sz>=0);
    assert(sz<MaxEntries);
    _size = sz;
  }
  accelerator_inline void  resize(size_type sz,const value &val) {
    assert(sz>=0);
    assert(sz<=MaxEntries);
    _size = sz;
    for(int s=0;s<sz;s++) _data[s]=val;
  }
  accelerator_inline pointer begin(void)                   { return _data; } 
  accelerator_inline pointer end  (void)                   { return &_data[_size]; } 
  accelerator_inline void push_back(const value &val)      { resize(_size+1); _data[_size-1] = val;}
  accelerator_inline AcceleratorVector()                   { _size = 0; }
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
#define GRID_MAX_SIMD              (16)

static constexpr int MaxDims = GRID_MAX_LATTICE_DIMENSION;

typedef AcceleratorVector<int,MaxDims> Coordinate;

inline std::ostream & operator<<(std::ostream &os, const Coordinate &v)
{
  os << "[";
  for(int s=0;s<v.size();s++) {
    os << v[s] << " ";
  }
  if (v.size() > 0) {
    os << "\b";
  }
  os << "]";
  return os;
}

NAMESPACE_END(Grid);
