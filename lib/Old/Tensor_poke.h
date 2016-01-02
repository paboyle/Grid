    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/Old/Tensor_poke.h

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
#ifndef GRID_MATH_POKE_H
#define GRID_MATH_POKE_H
namespace Grid {

//////////////////////////////////////////////////////////////////////////////
// Poke a specific index; 
//////////////////////////////////////////////////////////////////////////////
#if 0
// Scalar poke
template<int Level,class vtype,typename std::enable_if< iScalar<vtype>::TensorLevel == Level >::type * =nullptr> inline 
  void pokeIndex(iScalar<vtype> &ret, const iScalar<vtype> &arg)
{
  ret._internal = arg._internal;
}
// Vector poke, one index
template<int Level,class vtype,int N,typename std::enable_if< iScalar<vtype>::TensorLevel == Level >::type * =nullptr> inline 
  void pokeIndex(iVector<vtype,N> &ret, const iScalar<vtype> &arg,int i)
{
  ret._internal[i] = arg._internal;
}
//Matrix poke, two indices
template<int Level,class vtype,int N,typename std::enable_if< iScalar<vtype>::TensorLevel == Level >::type * =nullptr> inline 
  void pokeIndex(iMatrix<vtype,N> &ret, const iScalar<vtype> &arg,int i,int j)
{
  ret._internal[i][j] = arg._internal;
}

/////////////
// No match poke for scalar,vector,matrix must forward on either 0,1,2 args. Must have 9 routines with notvalue
/////////////
// scalar
template<int Level,class vtype,typename std::enable_if< iScalar<vtype>::TensorLevel != Level >::type * =nullptr> inline 
void pokeIndex(iScalar<vtype> &ret, const iScalar<decltype(peekIndex<Level>(ret._internal))>  &arg)
{
  pokeIndex<Level>(ret._internal,arg._internal);
}
template<int Level,class vtype,typename std::enable_if< iScalar<vtype>::TensorLevel != Level >::type * =nullptr> inline 
  void pokeIndex(iScalar<vtype> &ret, const iScalar<decltype(peekIndex<Level>(ret._internal,0))> &arg, int i)
		 
{
  pokeIndex<Level>(ret._internal,arg._internal,i);
}
template<int Level,class vtype,typename std::enable_if< iScalar<vtype>::TensorLevel != Level >::type * =nullptr> inline 
  void pokeIndex(iScalar<vtype> &ret, const iScalar<decltype(peekIndex<Level>(ret._internal,0,0))> &arg,int i,int j)
{
  pokeIndex<Level>(ret._internal,arg._internal,i,j);
}

// Vector
template<int Level,class vtype,int N,typename std::enable_if< iScalar<vtype>::TensorLevel != Level >::type * =nullptr> inline 
  void pokeIndex(iVector<vtype,N> &ret, iVector<decltype(peekIndex<Level>(ret._internal)),N>  &arg)
{
  for(int ii=0;ii<N;ii++){
    pokeIndex<Level>(ret._internal[ii],arg._internal[ii]);
  }
}
template<int Level,class vtype,int N,typename std::enable_if< iScalar<vtype>::TensorLevel != Level >::type * =nullptr> inline 
  void pokeIndex(iVector<vtype,N> &ret, const iVector<decltype(peekIndex<Level>(ret._internal,0)),N> &arg,int i)
{
  for(int ii=0;ii<N;ii++){
    pokeIndex<Level>(ret._internal[ii],arg._internal[ii],i);
  }
}
template<int Level,class vtype,int N,typename std::enable_if< iScalar<vtype>::TensorLevel != Level >::type * =nullptr> inline 
  void pokeIndex(iVector<vtype,N> &ret, const iVector<decltype(peekIndex<Level>(ret._internal,0,0)),N> &arg,int i,int j)
{
  for(int ii=0;ii<N;ii++){
    pokeIndex<Level>(ret._internal[ii],arg._internal[ii],i,j);
  }
}

// Matrix
template<int Level,class vtype,int N,typename std::enable_if< iScalar<vtype>::TensorLevel != Level >::type * =nullptr> inline 
  void pokeIndex(iMatrix<vtype,N> &ret, const iMatrix<decltype(peekIndex<Level>(ret._internal)),N> &arg)		 
{
  for(int ii=0;ii<N;ii++){
  for(int jj=0;jj<N;jj++){
    pokeIndex<Level>(ret._internal[ii][jj],arg._internal[ii][jj]);
  }}
}
template<int Level,class vtype,int N,typename std::enable_if< iScalar<vtype>::TensorLevel != Level >::type * =nullptr> inline 
  void pokeIndex(iMatrix<vtype,N> &ret, const iMatrix<decltype(peekIndex<Level>(ret._internal,0)),N> &arg,int i)
{
  for(int ii=0;ii<N;ii++){
  for(int jj=0;jj<N;jj++){
    pokeIndex<Level>(ret._internal[ii][jj],arg._internal[ii][jj],i);
  }}
}
template<int Level,class vtype,int N,typename std::enable_if< iScalar<vtype>::TensorLevel != Level >::type * =nullptr> inline 
  void pokeIndex(iMatrix<vtype,N> &ret, const iMatrix<decltype(peekIndex<Level>(ret._internal,0,0)),N> &arg, int i,int j)
{
  for(int ii=0;ii<N;ii++){
  for(int jj=0;jj<N;jj++){
    pokeIndex<Level>(ret._internal[ii][jj],arg._internal[ii][jj],i,j);
  }}
}
#endif

}
#endif
