    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/tensors/Tensor_transpose.h

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
#ifndef GRID_MATH_TRANSPOSE_H
#define GRID_MATH_TRANSPOSE_H
namespace Grid {



/////////////////////////////////////////////////////////////////
// Transpose all indices
/////////////////////////////////////////////////////////////////

inline ComplexD transpose(ComplexD &rhs){  return rhs;}
inline ComplexF transpose(ComplexF &rhs){  return rhs;}
inline RealD transpose(RealD &rhs){  return rhs;}
inline RealF transpose(RealF &rhs){  return rhs;}

template<class vtype,int N>
  inline typename std::enable_if<isGridTensor<vtype>::value, iMatrix<vtype,N> >::type 
  transpose(iMatrix<vtype,N> arg)
  {
    iMatrix<vtype,N> ret;
    for(int i=0;i<N;i++){
      for(int j=0;j<N;j++){
	ret._internal[i][j] = transpose(arg._internal[j][i]); // NB recurses
      }}
    return ret;
  }
template<class vtype,int N>
  inline typename std::enable_if<isGridTensor<vtype>::notvalue, iMatrix<vtype,N> >::type 
  transpose(iMatrix<vtype,N> arg)
  {
    iMatrix<vtype,N> ret;
    for(int i=0;i<N;i++){
      for(int j=0;j<N;j++){
	ret._internal[i][j] = arg._internal[j][i]; // Stop recursion if not a tensor type
      }}
    return ret;
  }

template<class vtype>
  inline typename std::enable_if<isGridTensor<vtype>::value, iScalar<vtype> >::type 
  transpose(iScalar<vtype> arg)
  {
    iScalar<vtype> ret;
    ret._internal = transpose(arg._internal); // NB recurses
    return ret;
  }

template<class vtype>
  inline typename std::enable_if<isGridTensor<vtype>::notvalue, iScalar<vtype> >::type 
  transpose(iScalar<vtype> arg)
  {
    iScalar<vtype> ret;
    ret._internal = arg._internal; // NB recursion stops
    return ret;
  }


////////////////////////////////////////////////////////////////////////////////////////////
// Transpose a specific index; instructive to compare this style of recursion termination
// to that of adj; which is easiers?
////////////////////////////////////////////////////////////////////////////////////////////
#if 0
template<int Level,class vtype,int N> inline 
  typename std::enable_if<matchGridTensorIndex<iMatrix<vtype,N>,Level>::value, iMatrix<vtype,N> >::type 
transposeIndex (const iMatrix<vtype,N> &arg)
{
  iMatrix<vtype,N> ret;
  for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){
      ret._internal[i][j] = arg._internal[j][i]; 
  }}
  return ret;
}
// or not
template<int Level,class vtype,int N> inline 
typename std::enable_if<matchGridTensorIndex<iMatrix<vtype,N>,Level>::notvalue, iMatrix<vtype,N> >::type 
transposeIndex (const iMatrix<vtype,N> &arg)
{
  iMatrix<vtype,N> ret;
  for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){
      ret._internal[i][j] = transposeIndex<Level>(arg._internal[i][j]); 
  }}
  return ret;
}
template<int Level,class vtype> inline 
typename std::enable_if<matchGridTensorIndex<iScalar<vtype>,Level>::notvalue, iScalar<vtype> >::type 
transposeIndex (const iScalar<vtype> &arg)
{
  iScalar<vtype> ret;
  ret._internal=transposeIndex<Level>(arg._internal);
  return ret;
}
template<int Level,class vtype> inline 
typename std::enable_if<matchGridTensorIndex<iScalar<vtype>,Level>::value, iScalar<vtype> >::type 
transposeIndex (const iScalar<vtype> &arg)
{
  return arg;
}
#endif

}
#endif
