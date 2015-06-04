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
      for(int c1=0;c1<N;c1++){
	for(int c2=0;c2<N;c2++){
	  ret._internal[c1][c2]= (ret._internal[c1][c2] - adj(arg._internal[c2][c1]));
	  ret._internal[c1][c2] *= 0.5;
	}}
      //ret = (ret - adj(arg))*0.5;
      ret -= trace(ret)*factor;
      return ret;
    }

}
#endif
