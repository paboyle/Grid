    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/tensors/Tensor_unary.h

    Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: paboyle <paboyle@ph.ed.ac.uk>

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
#ifndef GRID_TENSOR_UNARY_H
#define GRID_TENSOR_UNARY_H
namespace Grid {

#define UNARY(func)\
template<class obj> inline auto func(const iScalar<obj> &z) -> iScalar<obj>\
{\
    iScalar<obj> ret;\
    ret._internal = func( (z._internal));\
    return ret;\
}\
template<class obj,int N> inline auto func(const iVector<obj,N> &z) -> iVector<obj,N>\
{\
    iVector<obj,N> ret;\
    for(int c1=0;c1<N;c1++){\
      ret._internal[c1] = func( (z._internal[c1]));\
    }\
    return ret;\
}\
template<class obj,int N> inline auto func(const iMatrix<obj,N> &z) -> iMatrix<obj,N>\
{\
    iMatrix<obj,N> ret;\
    for(int c1=0;c1<N;c1++){\
    for(int c2=0;c2<N;c2++){\
      ret._internal[c1][c2] = func( (z._internal[c1][c2]));\
    }}\
    return ret;\
}


#define BINARY_RSCALAR(func,scal)					\
template<class obj> inline iScalar<obj> func(const iScalar<obj> &z,scal y)	\
{\
    iScalar<obj> ret;\
    ret._internal = func(z._internal,y);	\
    return ret;\
}\
 template<class obj,int N> inline iVector<obj,N> func(const iVector<obj,N> &z,scal y) \
{\
    iVector<obj,N> ret;\
    for(int c1=0;c1<N;c1++){\
      ret._internal[c1] = func(z._internal[c1],y);	\
    }\
    return ret;\
}\
 template<class obj,int N> inline  iMatrix<obj,N> func(const iMatrix<obj,N> &z, scal y) \
{\
    iMatrix<obj,N> ret;\
    for(int c1=0;c1<N;c1++){\
    for(int c2=0;c2<N;c2++){\
      ret._internal[c1][c2] = func(z._internal[c1][c2],y);	\
    }}\
    return ret;\
}

UNARY(sqrt);
UNARY(rsqrt);
UNARY(sin);
UNARY(cos);
UNARY(log);
UNARY(exp);
UNARY(abs);
UNARY(Not);


template<class obj> inline auto toReal(const iScalar<obj> &z) -> typename iScalar<obj>::Realified
{
  typename iScalar<obj>::Realified ret;
  ret._internal = toReal(z._internal);
  return ret;
}
 template<class obj,int N> inline auto toReal(const iVector<obj,N> &z) -> typename iVector<obj,N>::Realified
{
  typename iVector<obj,N>::Realified ret;
  for(int c1=0;c1<N;c1++){  
    ret._internal[c1] = toReal(z._internal[c1]); 
  }
  return ret;
}
template<class obj,int N> inline auto toReal(const iMatrix<obj,N> &z) -> typename iMatrix<obj,N>::Realified
{
  typename iMatrix<obj,N>::Realified ret;
  for(int c1=0;c1<N;c1++){
  for(int c2=0;c2<N;c2++){
    ret._internal[c1][c2] = toReal(z._internal[c1][c2]);
  }}
  return ret;
}

template<class obj> inline auto toComplex(const iScalar<obj> &z) -> typename iScalar<obj>::Complexified
{
  typename iScalar<obj>::Complexified ret;
  ret._internal = toComplex(z._internal);
  return ret;
}
 template<class obj,int N> inline auto toComplex(const iVector<obj,N> &z) -> typename iVector<obj,N>::Complexified
{
  typename iVector<obj,N>::Complexified ret;
  for(int c1=0;c1<N;c1++){  
    ret._internal[c1] = toComplex(z._internal[c1]); 
  }
  return ret;
}
template<class obj,int N> inline auto toComplex(const iMatrix<obj,N> &z) -> typename iMatrix<obj,N>::Complexified
{
  typename iMatrix<obj,N>::Complexified ret;
  for(int c1=0;c1<N;c1++){
  for(int c2=0;c2<N;c2++){
    ret._internal[c1][c2] = toComplex(z._internal[c1][c2]);
  }}
  return ret;
}

BINARY_RSCALAR(div,Integer);
BINARY_RSCALAR(mod,Integer);
BINARY_RSCALAR(pow,RealD);

#undef UNARY
#undef BINARY_RSCALAR


}
#endif
