/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/lattice/Lattice_local.h

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
#ifndef GRID_LATTICE_LOCALREDUCTION_H
#define GRID_LATTICE_LOCALREDUCTION_H

///////////////////////////////////////////////
// localInner, localNorm, outerProduct
///////////////////////////////////////////////

NAMESPACE_BEGIN(Grid);

/////////////////////////////////////////////////////
// Non site, reduced locally reduced routines
/////////////////////////////////////////////////////
  
// localNorm2,
template<class vobj>
inline auto localNorm2 (const Lattice<vobj> &rhs)-> Lattice<typename vobj::tensor_reduced>
{
  Lattice<typename vobj::tensor_reduced> ret(rhs.Grid());
  auto rhs_v = rhs.View();
  auto ret_v = ret.View();
  accelerator_loop(ss,rhs_v,{
    ret_v[ss]=innerProduct(rhs_v[ss],rhs_v[ss]);
  });
  return ret;
}
  
// localInnerProduct
template<class vobj>
inline auto localInnerProduct (const Lattice<vobj> &lhs,const Lattice<vobj> &rhs) -> Lattice<typename vobj::tensor_reduced>
{
  Lattice<typename vobj::tensor_reduced> ret(rhs.Grid());
  auto lhs_v = lhs.View();
  auto rhs_v = rhs.View();
  auto ret_v = ret.View();
  accelerator_loop(ss,rhs_v,{
    ret_v[ss]=innerProduct(lhs_v[ss],rhs_v[ss]);
  });
  return ret;
}
  
// outerProduct Scalar x Scalar -> Scalar
//              Vector x Vector -> Matrix
template<class ll,class rr>
inline auto outerProduct (const Lattice<ll> &lhs,const Lattice<rr> &rhs) -> Lattice<decltype(outerProduct(ll(),rr()))>
{
  Lattice<decltype(outerProduct(ll(),rr()))> ret(rhs.Grid());
  auto lhs_v = lhs.View();
  auto rhs_v = rhs.View();
  auto ret_v = ret.View();
  accelerator_loop(ss,rhs_v,{
    ret_v[ss]=outerProduct(lhs_v[ss],rhs_v[ss]);
  });
  return ret;
}
NAMESPACE_END(Grid);
#endif
