#ifndef GRID_MATH_TA_H
#define GRID_MATH_TA_H
namespace Grid {

  /////////////////////////////////////////////// 
  // Ta function for scalar, vector, matrix
  /////////////////////////////////////////////// 
  inline ComplexF Ta( const ComplexF &arg){    return arg;}
  inline ComplexD Ta( const ComplexD &arg){    return arg;}
  inline RealF Ta( const RealF &arg){    return arg;}
  inline RealD Ta( const RealD &arg){    return arg;}


  template<class vtype> inline iScalar<vtype> Ta(const iScalar<vtype>&r)
    {
      iScalar<vtype> ret;
      ret._internal = Ta(r._internal);
      return ret;
    }
  template<class vtype,int N> inline iVector<vtype,N> Ta(const iVector<vtype,N>&r)
    {
      iVector<vtype,N> ret;
      for(int i=0;i<N;i++){
        ret._internal[i] = Ta(r._internal[i]);
      }
      return ret;
    }
  template<class vtype,int N> inline iMatrix<vtype,N> Ta(const iMatrix<vtype,N> &arg)
    {
      iMatrix<vtype,N> ret(arg);
      double factor = (1/(double)N);
      ret = (ret - adj(arg))*0.5;
      ret -= trace(ret)*factor;
      return ret;
    }


  /////////////////////////////////////////////// 
  // ProjectOnGroup function for scalar, vector, matrix
  /////////////////////////////////////////////// 


  template<class vtype> inline iScalar<vtype> ProjectOnGroup(const iScalar<vtype>&r)
    {
      iScalar<vtype> ret;
      ret._internal = ProjectOnGroup(r._internal);
      return ret;
    }
  template<class vtype,int N> inline iVector<vtype,N> ProjectOnGroup(const iVector<vtype,N>&r)
    {
      iVector<vtype,N> ret;
      for(int i=0;i<N;i++){
        ret._internal[i] = ProjectOnGroup(r._internal[i]);
      }
      return ret;
    }
  template<class vtype,int N, typename std::enable_if< GridTypeMapper<vtype>::TensorLevel == 0 >::type * =nullptr> 
    inline iMatrix<vtype,N> ProjectOnGroup(const iMatrix<vtype,N> &arg)
    {
      // need a check for the group type?
      iMatrix<vtype,N> ret(arg);
      double nrm;
      for(int c1=0;c1<N;c1++){
	nrm = 0.0; 
	for(int c2=0;c2<N;c2++)
	  nrm = real(innerProduct(ret._internal[c1][c2],ret._internal[c1][c2]));
	nrm = 1.0/sqrt(nrm);
	for(int c2=0;c2<N;c2++)
	  ret._internal[c1][c2]*= nrm;
      
	for (int b=c1+1; b<N; ++b){
	  decltype(ret._internal[b][b]*ret._internal[b][b]) pr = 0.0;
	  for(int c=0; c<N; ++c)
	    pr += ret._internal[c1][c]*ret._internal[b][c];
	  
	  for(int c=0; c<N; ++c){
	    ret._internal[b][c] -= pr * ret._internal[c1][c];
	  }
	}
	  
      }
      return ret;
    }


  /////////////////////////////////////////////// 
  // Exponentiate function for scalar, vector, matrix
  /////////////////////////////////////////////// 


  template<class vtype> inline iScalar<vtype> Exponentiate(const iScalar<vtype>&r, double alpha,  int Nexp)
    {
      iScalar<vtype> ret;
      ret._internal = Exponentiate(r._internal, alpha, Nexp);
      return ret;
    }


  template<class vtype,int N, typename std::enable_if< GridTypeMapper<vtype>::TensorLevel == 0 >::type * =nullptr> 
    inline iMatrix<vtype,N> Exponentiate(const iMatrix<vtype,N> &arg, double alpha, int Nexp)
    {
      iMatrix<vtype,N> unit(1.0);
      iMatrix<vtype,N> temp(unit);
      
      for(int i=Nexp; i>=1;--i){
	temp *= alpha/double(i);
	temp = unit + temp*arg;
      }
      return ProjectOnGroup(temp);
    }



}
#endif
