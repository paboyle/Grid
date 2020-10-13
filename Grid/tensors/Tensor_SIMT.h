/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/tensors/Tensor_SIMT.h

    Copyright (C) 2015

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

#include <string.h>

NAMESPACE_BEGIN(Grid);


#ifndef GRID_SIMT
//////////////////////////////////////////
// Trivial mapping of vectors on host
//////////////////////////////////////////
template<class vobj> accelerator_inline
vobj coalescedRead(const vobj & __restrict__ vec,int lane=0)
{
  return vec;
}
template<class vobj> accelerator_inline
vobj coalescedReadPermute(const vobj & __restrict__ vec,int ptype,int doperm,int lane=0)
{
  if ( doperm ) {
    vobj ret;
    permute(ret,vec, ptype);
    return ret;
  } else { 
    return vec;
  }
}
template<class vobj> accelerator_inline
void coalescedWrite(vobj & __restrict__ vec,const vobj & __restrict__ extracted,int lane=0)
{
  vec = extracted;
}
template<class vobj> accelerator_inline
void coalescedWriteNonTemporal(vobj & __restrict__ vec,const vobj & __restrict__ extracted,int lane=0)
{
  vstream(vec, extracted);
}
#else

//////////////////////////////////////////
// Extract and insert slices on the GPU
//////////////////////////////////////////
template<class vobj> accelerator_inline
typename vobj::scalar_object coalescedRead(const vobj & __restrict__ vec,int lane=acceleratorSIMTlane(vobj::Nsimd()))
{
  return extractLane(lane,vec);
}
template<class vobj> accelerator_inline
typename vobj::scalar_object coalescedReadPermute(const vobj & __restrict__ vec,int ptype,int doperm,int lane=acceleratorSIMTlane(vobj::Nsimd()))
{
  int mask = vobj::Nsimd() >> (ptype + 1);		
  int plane= doperm ? lane ^ mask : lane;
  return extractLane(plane,vec);
}
template<class vobj> accelerator_inline
void coalescedWrite(vobj & __restrict__ vec,const typename vobj::scalar_object & __restrict__ extracted,int lane=acceleratorSIMTlane(vobj::Nsimd()))
{
  insertLane(lane,vec,extracted);
}
template<class vobj> accelerator_inline
void coalescedWriteNonTemporal(vobj & __restrict__ vec,const vobj & __restrict__ extracted,int lane=0)
{
  insertLane(lane,vec,extracted);
}
#endif


NAMESPACE_END(Grid);

