    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/tensors/Tensor_outer.h

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
#ifndef GRID_MATH_OUTER_H
#define GRID_MATH_OUTER_H
namespace Grid {
    ///////////////////////////////////////////////////////////////////////////////////////
    // outerProduct Scalar x Scalar -> Scalar
    //              Vector x Vector -> Matrix
    ///////////////////////////////////////////////////////////////////////////////////////

template<class l,class r,int N> inline
auto outerProduct (const iVector<l,N>& lhs,const iVector<r,N>& rhs) -> iMatrix<decltype(outerProduct(lhs._internal[0],rhs._internal[0])),N>
{
    typedef decltype(outerProduct(lhs._internal[0],rhs._internal[0])) ret_t;
    iMatrix<ret_t,N> ret;
    for(int c1=0;c1<N;c1++){
    for(int c2=0;c2<N;c2++){
        ret._internal[c1][c2] = outerProduct(lhs._internal[c1],rhs._internal[c2]);
    }}
    return ret;
}
template<class l,class r> inline
auto outerProduct (const iScalar<l>& lhs,const iScalar<r>& rhs) -> iScalar<decltype(outerProduct(lhs._internal,rhs._internal))>
{
    typedef decltype(outerProduct(lhs._internal,rhs._internal)) ret_t;
    iScalar<ret_t> ret;
    ret._internal = outerProduct(lhs._internal,rhs._internal);
    return ret;
}

inline ComplexF outerProduct(const ComplexF &l, const ComplexF& r)
{
  std::cout << "outer product taking conj "<<r<<" "<<conj(r)<<std::endl;
  return l*conj(r);
}
inline ComplexD outerProduct(const ComplexD &l, const ComplexD& r)
{
  std::cout << "outer product taking conj "<<r<<" "<<conj(r)<<std::endl;
  return l*conj(r);
}
inline RealF outerProduct(const RealF &l, const RealF& r)
{
  return l*r;
}
inline RealD outerProduct(const RealD &l, const RealD& r)
{
  return l*r;
}

}
#endif
