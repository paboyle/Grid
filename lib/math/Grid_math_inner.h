#ifndef GRID_MATH_INNER_H
#define GRID_MATH_INNER_H
namespace Grid {
    ///////////////////////////////////////////////////////////////////////////////////////
    // innerProduct Scalar x Scalar -> Scalar
    // innerProduct Vector x Vector -> Scalar
    // innerProduct Matrix x Matrix -> Scalar
    ///////////////////////////////////////////////////////////////////////////////////////
    template<class sobj> inline RealD norm2(const sobj &arg){
      typedef typename sobj::scalar_type scalar;
      decltype(innerProduct(arg,arg)) nrm;
      nrm = innerProduct(arg,arg);
      return real(nrm);
    }

    template<class l,class r,int N> inline
    auto innerProduct (const iVector<l,N>& lhs,const iVector<r,N>& rhs) -> iScalar<decltype(innerProduct(lhs._internal[0],rhs._internal[0]))>
    {
        typedef decltype(innerProduct(lhs._internal[0],rhs._internal[0])) ret_t;
        iScalar<ret_t> ret;
	ret=zero;
        for(int c1=0;c1<N;c1++){
            ret._internal += innerProduct(lhs._internal[c1],rhs._internal[c1]);
        }
        return ret;
    }
    template<class l,class r,int N> inline
    auto innerProduct (const iMatrix<l,N>& lhs,const iMatrix<r,N>& rhs) -> iScalar<decltype(innerProduct(lhs._internal[0][0],rhs._internal[0][0]))>
    {
        typedef decltype(innerProduct(lhs._internal[0][0],rhs._internal[0][0])) ret_t;
        iScalar<ret_t> ret;
        iScalar<ret_t> tmp;
	ret=zero;
        for(int c1=0;c1<N;c1++){
        for(int c2=0;c2<N;c2++){
	  ret._internal+=innerProduct(lhs._internal[c1][c2],rhs._internal[c1][c2]);
        }}
        return ret;
    }
    template<class l,class r> inline
    auto innerProduct (const iScalar<l>& lhs,const iScalar<r>& rhs) -> iScalar<decltype(innerProduct(lhs._internal,rhs._internal))>
    {
        typedef decltype(innerProduct(lhs._internal,rhs._internal)) ret_t;
        iScalar<ret_t> ret;
        ret._internal = innerProduct(lhs._internal,rhs._internal);
        return ret;
    }

}
#endif
