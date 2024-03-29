/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/tensors/Tensor_trace.h

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
#ifndef GRID_MATH_TRACE_H
#define GRID_MATH_TRACE_H

NAMESPACE_BEGIN(Grid);

//////////////////////////////////////////////////////////////////
// Traces: both all indices and a specific index. Indices must be
// either scalar or matrix
/////////////////////////////////////////////////////////////////

accelerator_inline ComplexF trace( const ComplexF &arg){    return arg;}
accelerator_inline ComplexD trace( const ComplexD &arg){    return arg;}
accelerator_inline RealF trace( const RealF &arg){    return arg;}
accelerator_inline RealD trace( const RealD &arg){    return arg;}

template<class vtype,int N>
accelerator_inline auto trace(const iMatrix<vtype,N> &arg) -> iScalar<decltype(trace(arg._internal[0][0]))>
{
  iScalar<decltype( trace(arg._internal[0][0] )) > ret;
  zeroit(ret._internal);
  for(int i=0;i<N;i++){
    ret._internal=ret._internal+trace(arg._internal[i][i]);
  }
  return ret;
}

template<class vtype>
accelerator_inline auto trace(const iScalar<vtype> &arg) -> iScalar<decltype(trace(arg._internal))>
{
  iScalar<decltype(trace(arg._internal))> ret;
  ret._internal=trace(arg._internal);
  return ret;
}

template<class vtype,int N>
accelerator_inline auto trace(const iVector<vtype,N> &arg) -> iVector<decltype(trace(arg._internal[0])),N>
{
  iVector<decltype(trace(arg._internal[0])),N> ret;
  for(int i=0;i<N;i++){
    ret._internal[i]=trace(arg._internal[i]);
  }
  return ret;
}
////////////////////////////
// Fast path traceProduct
////////////////////////////
template<class S1 , class S2, IfNotGridTensor<S1> = 0, IfNotGridTensor<S2> = 0>
accelerator_inline auto traceProduct( const S1 &arg1,const S2 &arg2)
  -> decltype(arg1*arg2)
{
  return arg1*arg2;
}

template<class vtype,class rtype,int N >
accelerator_inline auto traceProduct(const iMatrix<vtype,N> &arg1,const iMatrix<rtype,N> &arg2) -> iScalar<decltype(trace(arg1._internal[0][0]*arg2._internal[0][0]))>
{
  iScalar<decltype( trace(arg1._internal[0][0]*arg2._internal[0][0] )) > ret;
  zeroit(ret._internal);
  for(int i=0;i<N;i++){
  for(int j=0;j<N;j++){
    ret._internal=ret._internal+traceProduct(arg1._internal[i][j],arg2._internal[j][i]);
  }}
  return ret;
}

template<class vtype,class rtype >
accelerator_inline auto traceProduct(const iScalar<vtype> &arg1,const iScalar<rtype> &arg2) -> iScalar<decltype(trace(arg1._internal*arg2._internal))>
{
  iScalar<decltype(trace(arg1._internal*arg2._internal))> ret;
  ret._internal=traceProduct(arg1._internal,arg2._internal);
  return ret;
}

NAMESPACE_END(Grid);

#endif
