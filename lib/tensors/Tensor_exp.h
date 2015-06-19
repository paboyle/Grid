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
