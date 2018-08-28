    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/tensors/Tensor_reality.h

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: neo <cossu@post.kek.jp>

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
#ifndef GRID_MATH_REALITY_H
#define GRID_MATH_REALITY_H
namespace Grid {

/////////////////////////////////////////////// 
// multiply by I; make recursive.
/////////////////////////////////////////////// 
template<class vtype> inline iScalar<vtype> timesI(const iScalar<vtype>&r) 
{
    iScalar<vtype> ret;
    timesI(ret._internal,r._internal);
    return ret;
}
template<class vtype,int N> inline iVector<vtype,N> timesI(const iVector<vtype,N>&r) 
{
  iVector<vtype,N> ret;
  for(int i=0;i<N;i++){
    timesI(ret._internal[i],r._internal[i]);
  }
  return ret;
}
template<class vtype,int N> inline iMatrix<vtype,N> timesI(const iMatrix<vtype,N>&r)
{
  iMatrix<vtype,N> ret;
  for(int i=0;i<N;i++){
  for(int j=0;j<N;j++){
    timesI(ret._internal[i][j],r._internal[i][j]);
  }}
  return ret;
}

template<class vtype> inline void timesI(iScalar<vtype> &ret,const iScalar<vtype>&r) 
{
  timesI(ret._internal,r._internal);
}
template<class vtype,int N> inline void timesI(iVector<vtype,N> &ret,const iVector<vtype,N>&r) 
{
  for(int i=0;i<N;i++){
    timesI(ret._internal[i],r._internal[i]);
  }
}
template<class vtype,int N> inline void  timesI(iMatrix<vtype,N> &ret,const iMatrix<vtype,N>&r)
{
  for(int i=0;i<N;i++){
  for(int j=0;j<N;j++){
    timesI(ret._internal[i][j],r._internal[i][j]);
  }}
}


template<class vtype> inline iScalar<vtype> timesMinusI(const iScalar<vtype>&r) 
{
    iScalar<vtype> ret;
    timesMinusI(ret._internal,r._internal);
    return ret;
}
template<class vtype,int N> inline iVector<vtype,N> timesMinusI(const iVector<vtype,N>&r) 
{
  iVector<vtype,N> ret;
  for(int i=0;i<N;i++){
    timesMinusI(ret._internal[i],r._internal[i]);
  }
  return ret;
}
template<class vtype,int N> inline iMatrix<vtype,N> timesMinusI(const iMatrix<vtype,N>&r)
{
  iMatrix<vtype,N> ret;
  for(int i=0;i<N;i++){
  for(int j=0;j<N;j++){
    timesMinusI(ret._internal[i][j],r._internal[i][j]);
  }}
  return ret;
}

template<class vtype>  inline void timesMinusI(iScalar<vtype> &ret,const iScalar<vtype>&r) 
{
  timesMinusI(ret._internal,r._internal);
}
template<class vtype,int N> inline void timesMinusI(iVector<vtype,N> &ret,const iVector<vtype,N>&r) 
{
  for(int i=0;i<N;i++){
    timesMinusI(ret._internal[i],r._internal[i]);
  }
}
template<class vtype,int N> inline void  timesMinusI(iMatrix<vtype,N> &ret,const iMatrix<vtype,N>&r)
{
  for(int i=0;i<N;i++){
  for(int j=0;j<N;j++){
    timesMinusI(ret._internal[i][j],r._internal[i][j]);
  }}
}


/////////////////////////////////////////////// 
// Conj function for scalar, vector, matrix
/////////////////////////////////////////////// 
template<class vtype> inline iScalar<vtype> conjugate(const iScalar<vtype>&r)
{
    iScalar<vtype> ret;
    ret._internal = conjugate(r._internal);
    return ret;
}
template<class vtype,int N> inline iVector<vtype,N> conjugate(const iVector<vtype,N>&r)
{
  iVector<vtype,N> ret;
  for(int i=0;i<N;i++){
    ret._internal[i] = conjugate(r._internal[i]);
  }
  return ret;
}
template<class vtype,int N> inline iMatrix<vtype,N> conjugate(const iMatrix<vtype,N>&r)
{
  iMatrix<vtype,N> ret;
  for(int i=0;i<N;i++){
  for(int j=0;j<N;j++){
    ret._internal[i][j] = conjugate(r._internal[i][j]);
  }}
  return ret;
}

/////////////////////////////////////////////// 
// Adj function for scalar, vector, matrix
/////////////////////////////////////////////// 
template<class vtype> inline iScalar<vtype> adj(const iScalar<vtype>&r)
{
    iScalar<vtype> ret;
    ret._internal = adj(r._internal);
    return ret;
}
template<class vtype,int N> inline iVector<vtype,N> adj(const iVector<vtype,N>&r)
{
    iVector<vtype,N> ret;
    for(int i=0;i<N;i++){
        ret._internal[i] = adj(r._internal[i]);
    }
    return ret;
}
template<class vtype,int N> inline iMatrix<vtype,N> adj(const iMatrix<vtype,N> &arg)
{
    iMatrix<vtype,N> ret;
    for(int c1=0;c1<N;c1++){
    for(int c2=0;c2<N;c2++){
        ret._internal[c1][c2]=adj(arg._internal[c2][c1]);
    }}
    return ret;
}






/////////////////////////////////////////////////////////////////
// Can only take the real/imag part of scalar objects, since
// lattice objects of different complex nature are non-conformable.
/////////////////////////////////////////////////////////////////
template<class itype> inline auto real(const iScalar<itype> &z) -> iScalar<decltype(real(z._internal))>
{
    iScalar<decltype(real(z._internal))> ret;
    ret._internal = real(z._internal);
    return ret;
}
template<class itype,int N> inline auto real(const iMatrix<itype,N> &z) -> iMatrix<decltype(real(z._internal[0][0])),N>
{
    iMatrix<decltype(real(z._internal[0][0])),N> ret;
    for(int c1=0;c1<N;c1++){
    for(int c2=0;c2<N;c2++){
        ret._internal[c1][c2] = real(z._internal[c1][c2]);
    }}
    return ret;
}
template<class itype,int N> inline auto real(const iVector<itype,N> &z) -> iVector<decltype(real(z._internal[0])),N>
{
    iVector<decltype(real(z._internal[0])),N> ret;
    for(int c1=0;c1<N;c1++){
        ret._internal[c1] = real(z._internal[c1]);
    }
    return ret;
}
    
template<class itype> inline auto imag(const iScalar<itype> &z) -> iScalar<decltype(imag(z._internal))>
{
    iScalar<decltype(imag(z._internal))> ret;
    ret._internal = imag(z._internal);
    return ret;
}
template<class itype,int N> inline auto imag(const iMatrix<itype,N> &z) -> iMatrix<decltype(imag(z._internal[0][0])),N>
{
    iMatrix<decltype(imag(z._internal[0][0])),N> ret;
    for(int c1=0;c1<N;c1++){
    for(int c2=0;c2<N;c2++){
        ret._internal[c1][c2] = imag(z._internal[c1][c2]);
    }}
    return ret;
}
template<class itype,int N> inline auto imag(const iVector<itype,N> &z) -> iVector<decltype(imag(z._internal[0])),N>
{
    iVector<decltype(imag(z._internal[0])),N> ret;
    for(int c1=0;c1<N;c1++){
        ret._internal[c1] = imag(z._internal[c1]);
    }
    return ret;
}


}
#endif
