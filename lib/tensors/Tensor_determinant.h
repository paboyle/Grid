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
