    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/tensors/Tensor_logical.h

    Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>

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
#ifndef GRID_TENSOR_LOGICAL_H
#define GRID_TENSOR_LOGICAL_H

namespace Grid {

#define LOGICAL_BINOP(Op)\
template<class v> strong_inline iScalar<v> operator Op (const iScalar<v>& lhs,const iScalar<v>& rhs) \
{\
  iScalar<v> ret;\
  ret._internal = lhs._internal Op rhs._internal ;\
  return ret;\
}\
template<class l> strong_inline iScalar<l> operator Op (const iScalar<l>& lhs,Integer rhs) \
{\
  typename iScalar<l>::scalar_type t; t=rhs;\
  typename iScalar<l>::tensor_reduced srhs; srhs=t;\
  return lhs Op srhs;\
}\
template<class l> strong_inline iScalar<l> operator Op (Integer lhs,const iScalar<l>& rhs) \
{\
  typename iScalar<l>::scalar_type t;t=lhs;\
  typename iScalar<l>::tensor_reduced slhs;slhs=t;\
  return slhs Op rhs;\
}

LOGICAL_BINOP(|);
LOGICAL_BINOP(&);
LOGICAL_BINOP(||);
LOGICAL_BINOP(&&);

}
#endif
