#ifndef GRID_MATH_ARITH_MAC_H
#define GRID_MATH_ARITH_MAC_H

namespace Grid {


///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////// MAC         ///////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////

    ///////////////////////////
    // Legal multiplication table
    ///////////////////////////
    // scal x scal = scal
    // mat x  mat  = mat
    // mat  x scal = mat
    // scal x mat  = mat
    // mat  x vec  = vec
    // vec  x scal = vec
    // scal x vec  = vec
    ///////////////////////////
template<class rtype,class vtype,class mtype>
inline  void mac(iScalar<rtype> * __restrict__ ret,const iScalar<vtype> * __restrict__ lhs,const iScalar<mtype> * __restrict__ rhs)
{
    mac(&ret->_internal,&lhs->_internal,&rhs->_internal);
}
template<class rrtype,class ltype,class rtype,int N>
inline void mac(iMatrix<rrtype,N> * __restrict__ ret,const iMatrix<ltype,N> * __restrict__ lhs,const iMatrix<rtype,N> * __restrict__ rhs){
    for(int c2=0;c2<N;c2++){
    for(int c1=0;c1<N;c1++){
    for(int c3=0;c3<N;c3++){
        mac(&ret->_internal[c1][c2],&lhs->_internal[c1][c3],&rhs->_internal[c3][c2]);
    }}}
    return;
}
template<class rrtype,class ltype,class rtype,int N>
inline void mac(iMatrix<rrtype,N> * __restrict__ ret,const iMatrix<ltype,N> * __restrict__ lhs,const iScalar<rtype> * __restrict__ rhs){
    for(int c1=0;c1<N;c1++){
    for(int c2=0;c2<N;c2++){
        mac(&ret->_internal[c1][c2],&lhs->_internal[c1][c2],&rhs->_internal);
    }}
    return;
}
template<class rrtype,class ltype,class rtype,int N>
inline void mac(iMatrix<rrtype,N> * __restrict__ ret,const iScalar<ltype> * __restrict__ lhs,const iMatrix<rtype,N> * __restrict__ rhs){
    for(int c1=0;c1<N;c1++){
    for(int c2=0;c2<N;c2++){
        mac(&ret->_internal[c1][c2],&lhs->_internal,&rhs->_internal[c1][c2]);
    }}
    return;
}
template<class rrtype,class ltype,class rtype,int N>
inline void mac(iVector<rrtype,N> * __restrict__ ret,const iMatrix<ltype,N> * __restrict__ lhs,const iVector<rtype,N> * __restrict__ rhs)
{
    for(int c1=0;c1<N;c1++){
    for(int c2=0;c2<N;c2++){
        mac(&ret->_internal[c1],&lhs->_internal[c1][c2],&rhs->_internal[c2]);
    }}
    return;
}
template<class rrtype,class ltype,class rtype,int N>
inline void mac(iVector<rrtype,N> * __restrict__ ret,const iScalar<ltype> * __restrict__ lhs,const iVector<rtype,N> * __restrict__ rhs)
{
    for(int c1=0;c1<N;c1++){
        mac(&ret->_internal[c1],&lhs->_internal,&rhs->_internal[c1]);
    }
    return;
}
template<class rrtype,class ltype,class rtype,int N>
inline void mac(iVector<rrtype,N> * __restrict__ ret,const iVector<ltype,N> * __restrict__ lhs,const iScalar<rtype> * __restrict__ rhs)
{
    for(int c1=0;c1<N;c1++){
        mac(&ret->_internal[c1],&lhs->_internal[c1],&rhs->_internal);
    }
    return;
}
}

#endif
