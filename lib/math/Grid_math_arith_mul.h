#ifndef GRID_MATH_ARITH_MUL_H
#define GRID_MATH_ARITH_MUL_H

namespace Grid {


    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////// MUL         ///////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    
template<class rtype,class vtype,class mtype>
inline void mult(iScalar<rtype> * __restrict__ ret,const iScalar<mtype> * __restrict__ lhs,const iScalar<vtype> * __restrict__ rhs){
    mult(&ret->_internal,&lhs->_internal,&rhs->_internal);
}

template<class rrtype,class ltype,class rtype,int N>
inline void mult(iMatrix<rrtype,N> * __restrict__ ret,const iMatrix<ltype,N> * __restrict__ lhs,const iMatrix<rtype,N> * __restrict__ rhs){
    for(int c2=0;c2<N;c2++){
    for(int c1=0;c1<N;c1++){
        mult(&ret->_internal[c1][c2],&lhs->_internal[c1][0],&rhs->_internal[0][c2]);
        for(int c3=1;c3<N;c3++){
            mac(&ret->_internal[c1][c2],&lhs->_internal[c1][c3],&rhs->_internal[c3][c2]);
        }
    }}
    return;
}
template<class rrtype,class ltype,class rtype,int N>
inline void mult(iMatrix<rrtype,N> * __restrict__ ret,const iMatrix<ltype,N> * __restrict__ lhs,const iScalar<rtype> * __restrict__ rhs){
    for(int c2=0;c2<N;c2++){
    for(int c1=0;c1<N;c1++){
        mult(&ret->_internal[c1][c2],&lhs->_internal[c1][c2],&rhs->_internal);
    }}
    return;
}

template<class rrtype,class ltype,class rtype, int N>
inline void mult(iMatrix<rrtype,N> * __restrict__ ret,const iScalar<ltype>   * __restrict__ lhs,const iMatrix<rtype,N> * __restrict__ rhs){
    for(int c2=0;c2<N;c2++){
    for(int c1=0;c1<N;c1++){
        mult(&ret->_internal[c1][c2],&lhs->_internal,&rhs->_internal[c1][c2]);
    }}
    return;
}
// Matrix left multiplies vector
template<class rtype,class vtype,class mtype,int N>
inline void mult(iVector<rtype,N> * __restrict__ ret,const iMatrix<mtype,N> * __restrict__ lhs,const iVector<vtype,N> * __restrict__ rhs)
{
    for(int c1=0;c1<N;c1++){
        mult(&ret->_internal[c1],&lhs->_internal[c1][0],&rhs->_internal[0]);
        for(int c2=1;c2<N;c2++){
            mac(&ret->_internal[c1],&lhs->_internal[c1][c2],&rhs->_internal[c2]);
        }
    }
    return;
}
template<class rtype,class vtype,class mtype,int N>
inline void mult(iVector<rtype,N> * __restrict__ ret,
                 const iScalar<mtype>   * __restrict__ lhs,
                 const iVector<vtype,N> * __restrict__ rhs){
    for(int c1=0;c1<N;c1++){
        mult(&ret->_internal[c1],&lhs->_internal,&rhs->_internal[c1]);
    }
}
template<class rtype,class vtype,class mtype,int N>
inline void mult(iVector<rtype,N> * __restrict__ ret,
                 const iVector<vtype,N> * __restrict__ rhs,
                 const iScalar<mtype> * __restrict__ lhs){
    mult(ret,lhs,rhs);
}
    


template<class rtype,class vtype,class mtype,int N> inline
iVector<rtype,N> operator * (const iMatrix<mtype,N>& lhs,const iVector<vtype,N>& rhs)
{
    iVector<rtype,N> ret;
    mult(&ret,&lhs,&rhs);
    return ret;
}

template<class rtype,class vtype,class mtype,int N> inline
iVector<rtype,N> operator * (const iScalar<mtype>& lhs,const iVector<vtype,N>& rhs)
{
    iVector<rtype,N> ret;
    mult(&ret,&lhs,&rhs);
    return ret;
}

template<class rtype,class vtype,class mtype,int N> inline
iVector<rtype,N> operator * (const iVector<mtype,N>& lhs,const iScalar<vtype>& rhs)
{
    iVector<rtype,N> ret;
    mult(&ret,&lhs,&rhs);
    return ret;
}
    
    //////////////////////////////////////////////////////////////////
    // Glue operators to mult routines. Must resolve return type cleverly from typeof(internal)
    // since nesting matrix<scalar> x matrix<matrix>-> matrix<matrix>
    // while         matrix<scalar> x matrix<scalar>-> matrix<scalar>
    // so return type depends on argument types in nasty way.
    //////////////////////////////////////////////////////////////////
    // scal x scal = scal
    // mat x  mat  = mat
    // mat  x scal = mat
    // scal x mat  = mat
    // mat  x vec  = vec
    // vec  x scal = vec
    // scal x vec  = vec
    //
    // We can special case scalar_type ??
template<class l,class r>
inline auto operator * (const iScalar<l>& lhs,const iScalar<r>& rhs) -> iScalar<decltype(lhs._internal * rhs._internal)>
{
    typedef iScalar<decltype(lhs._internal*rhs._internal)> ret_t;
    ret_t ret;
    mult(&ret,&lhs,&rhs);
    return ret;
}
template<class l,class r,int N> inline
auto operator * (const iMatrix<l,N>& lhs,const iMatrix<r,N>& rhs) -> iMatrix<decltype(lhs._internal[0][0]*rhs._internal[0][0]),N>
{
    typedef decltype(lhs._internal[0][0]*rhs._internal[0][0]) ret_t;
    iMatrix<ret_t,N> ret;
    mult(&ret,&lhs,&rhs);
    return ret;
}
template<class l,class r, int N> inline
auto operator * (const iMatrix<r,N>& lhs,const iScalar<l>& rhs) -> iMatrix<decltype(lhs._internal[0][0]*rhs._internal),N>
{
    typedef decltype(lhs._internal[0][0]*rhs._internal) ret_t;
        
    iMatrix<ret_t,N> ret;
    for(int c1=0;c1<N;c1++){
    for(int c2=0;c2<N;c2++){
        mult(&ret._internal[c1][c2],&lhs._internal[c1][c2],&rhs._internal);
    }}
    return ret;
}
template<class l,class r,int N> inline
auto operator * (const iScalar<l>& lhs,const iMatrix<r,N>& rhs) -> iMatrix<decltype(lhs._internal*rhs._internal[0][0]),N>
{
    typedef decltype(lhs._internal*rhs._internal[0][0]) ret_t;
    iMatrix<ret_t,N> ret;
    for(int c1=0;c1<N;c1++){
    for(int c2=0;c2<N;c2++){
        mult(&ret._internal[c1][c2],&lhs._internal,&rhs._internal[c1][c2]);
    }}
    return ret;
}
template<class l,class r,int N> inline
auto operator * (const iMatrix<l,N>& lhs,const iVector<r,N>& rhs) -> iVector<decltype(lhs._internal[0][0]*rhs._internal[0]),N>
{
    typedef decltype(lhs._internal[0][0]*rhs._internal[0]) ret_t;
    iVector<ret_t,N> ret;
    for(int c1=0;c1<N;c1++){
        mult(&ret._internal[c1],&lhs._internal[c1][0],&rhs._internal[0]);
        for(int c2=1;c2<N;c2++){
            mac(&ret._internal[c1],&lhs._internal[c1][c2],&rhs._internal[c2]);
        }
    }
    return ret;
}
template<class l,class r,int N> inline
auto operator * (const iScalar<l>& lhs,const iVector<r,N>& rhs) -> iVector<decltype(lhs._internal*rhs._internal[0]),N>
{
    typedef decltype(lhs._internal*rhs._internal[0]) ret_t;
    iVector<ret_t,N> ret;
    for(int c1=0;c1<N;c1++){
        mult(&ret._internal[c1],&lhs._internal,&rhs._internal[c1]);
    }
    return ret;
}
template<class l,class r,int N> inline
auto operator * (const iVector<l,N>& lhs,const iScalar<r>& rhs) -> iVector<decltype(lhs._internal[0]*rhs._internal),N>
{
    typedef decltype(lhs._internal[0]*rhs._internal) ret_t;
    iVector<ret_t,N> ret;
    for(int c1=0;c1<N;c1++){
        mult(&ret._internal[c1],&lhs._internal[c1],&rhs._internal);
    }
    return ret;
}


}

#endif
