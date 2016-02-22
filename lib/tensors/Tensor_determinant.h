    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/tensors/Tensor_determinant.h

    Copyright (C) 2015

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
#ifndef GRID_MATH_DET_H
#define GRID_MATH_DET_H
namespace Grid {
  /////////////////////////////////////////////// 
  // Determinant function for scalar, vector, matrix
  /////////////////////////////////////////////// 
  inline ComplexF Determinant( const ComplexF &arg){    return arg;}
  inline ComplexD Determinant( const ComplexD &arg){    return arg;}
  inline RealF Determinant( const RealF &arg){    return arg;}
  inline RealD Determinant( const RealD &arg){    return arg;}

  template<class vtype> inline auto Determinant(const iScalar<vtype>&r) -> iScalar<decltype(Determinant(r._internal))>
    {
      iScalar<decltype(Determinant(r._internal))> ret;
      ret._internal = Determinant(r._internal);
      return ret;
    }

  template<class vtype,int N, typename std::enable_if< GridTypeMapper<vtype>::TensorLevel == 0 >::type * =nullptr> 
    inline iScalar<vtype> Determinant(const iMatrix<vtype,N> &arg)
    {
      iMatrix<vtype,N> ret(arg);
      iScalar<vtype> det = vtype(1.0);
      /* Conversion of matrix to upper triangular */
      for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
	  if(j>i){
	    vtype ratio = ret._internal[j][i]/ret._internal[i][i];
	    for(int k = 0; k < N; k++){
	      ret._internal[j][k] -= ratio * ret._internal[i][k];
	    }
	  }
        }
      }      

      for(int i = 0; i < N; i++)
	det *= ret._internal[i][i];   

      return det;
    }



}
#endif
