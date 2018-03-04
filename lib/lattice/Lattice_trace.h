/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/lattice/Lattice_trace.h

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
#ifndef GRID_LATTICE_TRACE_H
#define GRID_LATTICE_TRACE_H

///////////////////////////////////////////////
// Tracing, transposing, peeking, poking
///////////////////////////////////////////////

NAMESPACE_BEGIN(Grid);

////////////////////////////////////////////////////////////////////////////////////////////////////
// Trace
////////////////////////////////////////////////////////////////////////////////////////////////////
template<class vobj>
inline auto trace(const Lattice<vobj> &lhs)  -> Lattice<decltype(trace(vobj()))>
{
  Lattice<decltype(trace(vobj()))> ret(lhs.Grid());
  auto ret_v = ret.View();
  auto lhs_v = lhs.View();
  accelerator_loop( ss, lhs_v, {
    ret_v[ss] = trace(lhs_v[ss]);
  });
  return ret;
};
    
////////////////////////////////////////////////////////////////////////////////////////////////////
// Trace Index level dependent operation
////////////////////////////////////////////////////////////////////////////////////////////////////
template<int Index,class vobj>
inline auto TraceIndex(const Lattice<vobj> &lhs) -> Lattice<decltype(traceIndex<Index>(vobj()))>
{
  Lattice<decltype(traceIndex<Index>(vobj()))> ret(lhs.Grid());
  auto ret_v = ret.View();
  auto lhs_v = lhs.View();
  accelerator_loop( ss, lhs_v, {
    ret_v[ss] = traceIndex<Index>(lhs_v[ss]);
  });
  return ret;
};

NAMESPACE_END(Grid);
#endif

