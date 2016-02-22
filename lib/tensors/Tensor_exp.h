    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/tensors/Tensor_exp.h

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
#ifndef GRID_MATH_EXP_H
#define GRID_MATH_EXP_H

#define DEFAULT_MAT_EXP 12

namespace Grid {

  /////////////////////////////////////////////// 
  // Exponentiate function for scalar, vector, matrix
  /////////////////////////////////////////////// 


  template<class vtype> inline iScalar<vtype> Exponentiate(const iScalar<vtype>&r, ComplexD alpha ,  Integer Nexp = DEFAULT_MAT_EXP)
    {
      iScalar<vtype> ret;
      ret._internal = Exponentiate(r._internal, alpha, Nexp);
      return ret;
    }


  template<class vtype,int N, typename std::enable_if< GridTypeMapper<vtype>::TensorLevel == 0 >::type * =nullptr> 
    inline iMatrix<vtype,N> Exponentiate(const iMatrix<vtype,N> &arg, ComplexD alpha  , Integer Nexp = DEFAULT_MAT_EXP )
    {
      iMatrix<vtype,N> unit(1.0);
      iMatrix<vtype,N> temp(unit);
      
      for(int i=Nexp; i>=1;--i){
	temp *= alpha/ComplexD(i);
	temp = unit + temp*arg;
      }
      
      return ProjectOnGroup(temp);//maybe not strictly necessary
      
    }



}
#endif
