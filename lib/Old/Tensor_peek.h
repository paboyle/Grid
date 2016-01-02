    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/Old/Tensor_peek.h

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
#ifndef GRID_MATH_PEEK_H
#define GRID_MATH_PEEK_H
namespace Grid {

//////////////////////////////////////////////////////////////////////////////
// Peek on a specific index; returns a scalar in that index, tensor inherits rest
//////////////////////////////////////////////////////////////////////////////
// If we hit the right index, return scalar with no further recursion

//template<int Level> inline ComplexF peekIndex(const ComplexF arg) { return arg;}
//template<int Level> inline ComplexD peekIndex(const ComplexD arg) { return arg;}
//template<int Level> inline RealF peekIndex(const RealF arg) { return arg;}
//template<int Level> inline RealD peekIndex(const RealD arg) { return arg;}
#if 0
// Scalar peek, no indices
template<int Level,class vtype,typename std::enable_if< iScalar<vtype>::TensorLevel == Level >::type * =nullptr> inline 
  auto peekIndex(const iScalar<vtype> &arg) ->  iScalar<vtype> 
{
  return arg;
}
// Vector peek, one index
template<int Level,class vtype,int N,typename std::enable_if< iScalar<vtype>::TensorLevel == Level >::type * =nullptr> inline 
  auto peekIndex(const iVector<vtype,N> &arg,int i) -> iScalar<vtype> // Index matches
{
  iScalar<vtype> ret;                              // return scalar
  ret._internal = arg._internal[i];
  return ret;
}
// Matrix peek, two indices
template<int Level,class vtype,int N,typename std::enable_if< iScalar<vtype>::TensorLevel == Level >::type * =nullptr> inline 
  auto peekIndex(const iMatrix<vtype,N> &arg,int i,int j) ->  iScalar<vtype>
{
  iScalar<vtype> ret;                              // return scalar
  ret._internal = arg._internal[i][j];
  return ret;
}

/////////////
// No match peek for scalar,vector,matrix must forward on either 0,1,2 args. Must have 9 routines with notvalue
/////////////
// scalar
template<int Level,class vtype,typename std::enable_if< iScalar<vtype>::TensorLevel != Level >::type * =nullptr> inline 
  auto peekIndex(const iScalar<vtype> &arg) -> iScalar<decltype(peekIndex<Level>(arg._internal))>
{
  iScalar<decltype(peekIndex<Level>(arg._internal))> ret;
  ret._internal= peekIndex<Level>(arg._internal);
  return ret;
}
template<int Level,class vtype, typename std::enable_if< iScalar<vtype>::TensorLevel != Level >::type * =nullptr> inline 
  auto peekIndex(const iScalar<vtype> &arg,int i) ->  iScalar<decltype(peekIndex<Level>(arg._internal,i))> 
{
  iScalar<decltype(peekIndex<Level>(arg._internal,i))> ret;
  ret._internal=peekIndex<Level>(arg._internal,i);
  return ret;
}
template<int Level,class vtype, typename std::enable_if< iScalar<vtype>::TensorLevel != Level >::type * =nullptr> inline 
  auto peekIndex(const iScalar<vtype> &arg,int i,int j) ->  iScalar<decltype(peekIndex<Level>(arg._internal,i,j))>
{
  iScalar<decltype(peekIndex<Level>(arg._internal,i,j))> ret;
  ret._internal=peekIndex<Level>(arg._internal,i,j);
  return ret;
}
// vector
template<int Level,class vtype,int N, typename std::enable_if< iScalar<vtype>::TensorLevel != Level >::type * =nullptr> inline 
auto peekIndex(const iVector<vtype,N> &arg) ->   iVector<decltype(peekIndex<Level>(arg._internal[0])),N>
{
  iVector<decltype(peekIndex<Level>(arg._internal[0])),N> ret;
  for(int ii=0;ii<N;ii++){
    ret._internal[ii]=peekIndex<Level>(arg._internal[ii]);
  }
  return ret;
}
template<int Level,class vtype,int N, typename std::enable_if< iScalar<vtype>::TensorLevel != Level >::type * =nullptr> inline 
  auto peekIndex(const iVector<vtype,N> &arg,int i) ->  iVector<decltype(peekIndex<Level>(arg._internal[0],i)),N>
{
  iVector<decltype(peekIndex<Level>(arg._internal[0],i)),N> ret;
  for(int ii=0;ii<N;ii++){
    ret._internal[ii]=peekIndex<Level>(arg._internal[ii],i);
  }
  return ret;
}
template<int Level,class vtype,int N, typename std::enable_if< iScalar<vtype>::TensorLevel != Level >::type * =nullptr> inline 
  auto peekIndex(const iVector<vtype,N> &arg,int i,int j) ->  iVector<decltype(peekIndex<Level>(arg._internal[0],i,j)),N> 
{
  iVector<decltype(peekIndex<Level>(arg._internal[0],i,j)),N> ret;
  for(int ii=0;ii<N;ii++){
    ret._internal[ii]=peekIndex<Level>(arg._internal[ii],i,j);
  }
  return ret;
}

// matrix
template<int Level,class vtype,int N, typename std::enable_if< iScalar<vtype>::TensorLevel != Level >::type * =nullptr> inline 
auto peekIndex(const iMatrix<vtype,N> &arg) ->   iMatrix<decltype(peekIndex<Level>(arg._internal[0][0])),N> 
{
  iMatrix<decltype(peekIndex<Level>(arg._internal[0][0])),N> ret;
  for(int ii=0;ii<N;ii++){
  for(int jj=0;jj<N;jj++){
    ret._internal[ii][jj]=peekIndex<Level>(arg._internal[ii][jj]);// Could avoid this because peeking a scalar is dumb
  }}
  return ret;
}
template<int Level,class vtype,int N, typename std::enable_if< iScalar<vtype>::TensorLevel != Level >::type * =nullptr> inline 
  auto peekIndex(const iMatrix<vtype,N> &arg,int i) ->   iMatrix<decltype(peekIndex<Level>(arg._internal[0][0],i)),N>
{
  iMatrix<decltype(peekIndex<Level>(arg._internal[0][0],i)),N> ret;
  for(int ii=0;ii<N;ii++){
  for(int jj=0;jj<N;jj++){
    ret._internal[ii][jj]=peekIndex<Level>(arg._internal[ii][jj],i);
  }}
  return ret;
}
template<int Level,class vtype,int N, typename std::enable_if< iScalar<vtype>::TensorLevel != Level >::type * =nullptr> inline 
  auto peekIndex(const iMatrix<vtype,N> &arg,int i,int j) ->   iMatrix<decltype(peekIndex<Level>(arg._internal[0][0],i,j)),N>
{
  iMatrix<decltype(peekIndex<Level>(arg._internal[0][0],i,j)),N> ret;
  for(int ii=0;ii<N;ii++){
  for(int jj=0;jj<N;jj++){
    ret._internal[ii][jj]=peekIndex<Level>(arg._internal[ii][jj],i,j);
  }}
  return ret;
}
#endif


}
#endif
