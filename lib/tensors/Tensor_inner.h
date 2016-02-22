    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/tensors/Tensor_inner.h

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
#ifndef GRID_MATH_INNER_H
#define GRID_MATH_INNER_H
namespace Grid {
    ///////////////////////////////////////////////////////////////////////////////////////
    // innerProduct Scalar x Scalar -> Scalar
    // innerProduct Vector x Vector -> Scalar
    // innerProduct Matrix x Matrix -> Scalar
    ///////////////////////////////////////////////////////////////////////////////////////
    template<class sobj> inline RealD norm2(const sobj &arg){
      typedef typename sobj::scalar_type scalar;
      decltype(innerProduct(arg,arg)) nrm;
      nrm = innerProduct(arg,arg);
      RealD ret = real(nrm);
      return ret;
    }

    template<class l,class r,int N> inline
    auto innerProduct (const iVector<l,N>& lhs,const iVector<r,N>& rhs) -> iScalar<decltype(innerProduct(lhs._internal[0],rhs._internal[0]))>
    {
        typedef decltype(innerProduct(lhs._internal[0],rhs._internal[0])) ret_t;
        iScalar<ret_t> ret;
	ret=zero;
        for(int c1=0;c1<N;c1++){
            ret._internal += innerProduct(lhs._internal[c1],rhs._internal[c1]);
        }
        return ret;
    }
    template<class l,class r,int N> inline
    auto innerProduct (const iMatrix<l,N>& lhs,const iMatrix<r,N>& rhs) -> iScalar<decltype(innerProduct(lhs._internal[0][0],rhs._internal[0][0]))>
    {
        typedef decltype(innerProduct(lhs._internal[0][0],rhs._internal[0][0])) ret_t;
        iScalar<ret_t> ret;
        iScalar<ret_t> tmp;
	ret=zero;
        for(int c1=0;c1<N;c1++){
        for(int c2=0;c2<N;c2++){
	  ret._internal+=innerProduct(lhs._internal[c1][c2],rhs._internal[c1][c2]);
        }}
        return ret;
    }
    template<class l,class r> inline
    auto innerProduct (const iScalar<l>& lhs,const iScalar<r>& rhs) -> iScalar<decltype(innerProduct(lhs._internal,rhs._internal))>
    {
        typedef decltype(innerProduct(lhs._internal,rhs._internal)) ret_t;
        iScalar<ret_t> ret;
        ret._internal = innerProduct(lhs._internal,rhs._internal);
        return ret;
    }

}
#endif
