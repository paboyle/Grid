    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/tensors/Tensor_arith_scalar.h

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
#ifndef GRID_MATH_ARITH_SCALAR_H
#define GRID_MATH_ARITH_SCALAR_H

namespace Grid {


//////////////////////////////////////////////////////////////////////////////////////////
// Must support native C++ types Integer, Complex, Real
//////////////////////////////////////////////////////////////////////////////////////////

// multiplication by fundamental scalar type
template<class l, class r> strong_inline typename
std::enable_if<std::is_constructible<typename iScalar<l>::scalar_type,
r>::value, iScalar<l>>::type operator * (const iScalar<l>& lhs, const r& rhs)
{
  typename iScalar<l>::scalar_type t;t=rhs;
  typename iScalar<l>::tensor_reduced srhs; srhs=t;
  return lhs*srhs;
}
template<class l, class r> strong_inline typename
std::enable_if<std::is_constructible<typename iScalar<r>::scalar_type,
l>::value, iScalar<r>>::type operator * (const l& lhs, const iScalar<r>& rhs) {  return rhs*lhs; }

template<class l, class r, int N> strong_inline typename
std::enable_if<std::is_constructible<typename iScalar<l>::scalar_type,
r>::value, iVector<l,N>>::type operator * (const iVector<l,N>& lhs,const r& rhs)
{
  typename iScalar<l>::scalar_type t;t=rhs;
  typename iVector<l,N>::tensor_reduced srhs; srhs=t;
  return lhs*srhs;
}
template<class l, class r, int N> strong_inline typename
std::enable_if<std::is_constructible<typename iScalar<r>::scalar_type,
l>::value, iVector<r,N>>::type operator * (const l& lhs,const iVector<r,N>& rhs) {  return rhs*lhs; }

template<class l, class r, int N> strong_inline typename
std::enable_if<std::is_constructible<typename iScalar<l>::scalar_type,
r>::value, iMatrix<l,N>>::type operator * (const iMatrix<l,N>& lhs,const r& rhs)
{
  typename iScalar<l>::scalar_type t;t=rhs;
  typename iMatrix<l,N>::tensor_reduced srhs; srhs=t;
  return lhs*srhs;
}
template<class l, class r, int N> strong_inline typename
std::enable_if<std::is_constructible<typename iScalar<r>::scalar_type,
l>::value, iMatrix<r,N>>::type operator * (const l& lhs,const iMatrix<r,N>& rhs) {  return rhs*lhs; }

///////////////////////////////////////////////////////////////////////////////////////////////
// addition by fundamental scalar type applies to matrix(down diag) and scalar
///////////////////////////////////////////////////////////////////////////////////////////////
template<class l, class r> strong_inline typename
std::enable_if<std::is_constructible<typename iScalar<l>::scalar_type,
r>::value, iScalar<l>>::type operator + (const iScalar<l>& lhs, const r& rhs)
{
  typename iScalar<l>::scalar_type t;t=rhs;
  typename iScalar<l>::tensor_reduced srhs; srhs=t;
  return lhs+srhs;
}
template<class l, class r> strong_inline typename
std::enable_if<std::is_constructible<typename iScalar<r>::scalar_type,
l>::value, iScalar<r>>::type operator + (const l& lhs, const iScalar<r>& rhs) {  return rhs+lhs; }

template<class l, class r, int N> strong_inline typename
std::enable_if<std::is_constructible<typename iScalar<l>::scalar_type,
r>::value, iMatrix<l,N>>::type operator + (const iMatrix<l,N>& lhs,const r& rhs)
{
  typename iScalar<l>::scalar_type t;t=rhs;
  typename iMatrix<l,N>::tensor_reduced srhs; srhs=t;
  return lhs+srhs;
}
template<class l, class r, int N> strong_inline typename
std::enable_if<std::is_constructible<typename iScalar<r>::scalar_type,
l>::value, iMatrix<r,N>>::type operator + (const l& lhs,const iMatrix<r,N>& rhs) {  return rhs+lhs; }

///////////////////////////////////////////////////////////////////////////////////////////////
// subtraction of fundamental scalar type applies to matrix(down diag) and scalar
///////////////////////////////////////////////////////////////////////////////////////////////
template<class l, class r> strong_inline typename
std::enable_if<std::is_constructible<typename iScalar<l>::scalar_type,
r>::value, iScalar<l>>::type operator - (const iScalar<l>& lhs, const r& rhs)
{
  typename iScalar<l>::scalar_type t;t=rhs;
  typename iScalar<l>::tensor_reduced srhs; srhs=t;
  return lhs-srhs;
}
template<class l, class r> strong_inline typename
std::enable_if<std::is_constructible<typename iScalar<r>::scalar_type,
l>::value, iScalar<r>>::type operator - (const l& lhs, const iScalar<r>& rhs)
{
  typename iScalar<r>::scalar_type t;t=lhs;
  typename iScalar<r>::tensor_reduced slhs; slhs=t;
  return slhs-rhs;
}

template<class l, class r, int N> strong_inline typename
std::enable_if<std::is_constructible<typename iScalar<l>::scalar_type,
r>::value, iMatrix<l,N>>::type operator - (const iMatrix<l,N>& lhs,const r& rhs)
{
  typename iScalar<l>::scalar_type t;t=rhs;
  typename iMatrix<l,N>::tensor_reduced srhs; srhs=t;
  return lhs-srhs;
}
template<class l, class r, int N> strong_inline typename
std::enable_if<std::is_constructible<typename iScalar<r>::scalar_type,
l>::value, iMatrix<r,N>>::type operator - (const l& lhs,const iMatrix<r,N>& rhs)
{
  typename iScalar<r>::scalar_type t;t=lhs;
  typename iMatrix<r,N>::tensor_reduced slhs; slhs=t;
  return slhs-rhs;
}


}
#endif
