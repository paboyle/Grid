//
//  Grid.cpp
//  simd
//
//  Created by Peter Boyle on 09/05/2014.
//  Copyright (c) 2014 University of Edinburgh. All rights reserved.
//


#ifndef GRID_V3_H
#define GRID_V3_H

#include <stdio.h>
#include <complex>
#include <vector>
#include <iostream>
#include <cassert>
#include <random>
#include <functional>
#include <stdlib.h>

#ifdef OMP
#include <omp.h>
#endif

#ifdef MAC
#include <malloc/malloc.h>
#else
#include <malloc.h>
#endif

#ifndef POOH
#define ALIGN_DIRECTIVE(A) __attribute__ ((aligned(A)))
#else
#define ALIGN_DIRECTIVE(A) __declspec(align(A))
#endif

#if defined(AVX1) || defined (AVX2)
#include <immintrin.h>
#define SIMDalign ALIGN_DIRECTIVE(32)
#endif
#ifdef SSE2
#include <pmmintrin.h>
#define SIMDalign ALIGN_DIRECTIVE(16)
#endif
#ifdef AVX512
#include <immintrin.h>
#define SIMDalign ALIGN_DIRECTIVE(64)
#endif

#include <sys/time.h>
#include <stdio.h>

namespace dpo {

  void Grid_init(void);


inline double usecond(void)
{
    struct timeval tv;
    gettimeofday(&tv,NULL);
    return 1.0*tv.tv_usec + 1.0e6*tv.tv_sec;
}
    
    typedef  float  RealF;
    typedef  double RealD;
    typedef  RealF  Real;
    
    typedef std::complex<RealF> ComplexF;
    typedef std::complex<RealD> ComplexD;
    typedef std::complex<Real>  Complex;
    

    class Zero{};
    static Zero zero;
    template<class itype> inline void ZeroIt(itype &arg){ arg=zero;};
    template<>            inline void ZeroIt(ComplexF &arg){ arg=0; };
    template<>            inline void ZeroIt(ComplexD &arg){ arg=0; };
    template<>            inline void ZeroIt(RealF &arg){ arg=0; };
    template<>            inline void ZeroIt(RealD &arg){ arg=0; };

    // TODO
    //
    // Base class to share common code between vRealF, VComplexF etc...
    //
    // lattice Broad cast assignment
    //
    // where() support
    // implement with masks, and/or? Type of the mask & boolean support?
    //
    // Unary functions
    // cos,sin, tan, acos, asin, cosh, acosh, tanh, sinh, // Scalar<vReal> only arg
    // exp, log, sqrt, fabs
    //
    // transposeColor, transposeSpin,
    // adjColor, adjSpin,
    // traceColor, traceSpin.
    // peekColor, peekSpin + pokeColor PokeSpin
    //
    // copyMask.
    //
    // localMaxAbs
    //
    // norm2,
    // sumMulti equivalent.
    // Fourier transform equivalent.
    //
    
    ////////////////////////////////////////////////////////////////////////////////
    //Provide support functions for basic real and complex data types required by dpo
    //Single and double precision versions. Should be able to template this once only.
    ////////////////////////////////////////////////////////////////////////////////
    
    inline void mac (ComplexD * __restrict__ y,const ComplexD * __restrict__ a,const ComplexD *__restrict__ x){ *y = (*a) * (*x)+(*y); };
    inline void mult(ComplexD * __restrict__ y,const ComplexD * __restrict__ l,const ComplexD *__restrict__ r){ *y = (*l) * (*r);}
    inline void sub (ComplexD * __restrict__ y,const ComplexD * __restrict__ l,const ComplexD *__restrict__ r){ *y = (*l) - (*r);}
    inline void add (ComplexD * __restrict__ y,const ComplexD * __restrict__ l,const ComplexD *__restrict__ r){ *y = (*l) + (*r);}
    inline ComplexD adj(const ComplexD& r){ return(conj(r)); }
    // conj already supported for complex
    
    inline void mac (ComplexF * __restrict__ y,const ComplexF * __restrict__ a,const ComplexF *__restrict__ x){ *y = (*a) * (*x)+(*y); }
    inline void mult(ComplexF * __restrict__ y,const ComplexF * __restrict__ l,const ComplexF *__restrict__ r){ *y = (*l) * (*r); }
    inline void sub (ComplexF * __restrict__ y,const ComplexF * __restrict__ l,const ComplexF *__restrict__ r){ *y = (*l) - (*r); }
    inline void add (ComplexF * __restrict__ y,const ComplexF * __restrict__ l,const ComplexF *__restrict__ r){ *y = (*l) + (*r); }
    inline Complex  adj(const Complex& r ){ return(conj(r)); }
    //conj already supported for complex
    
    inline void mac (RealD * __restrict__ y,const RealD * __restrict__ a,const RealD *__restrict__ x){  *y = (*a) * (*x)+(*y);}
    inline void mult(RealD * __restrict__ y,const RealD * __restrict__ l,const RealD *__restrict__ r){ *y = (*l) * (*r);}
    inline void sub (RealD * __restrict__ y,const RealD * __restrict__ l,const RealD *__restrict__ r){ *y = (*l) - (*r);}
    inline void add (RealD * __restrict__ y,const RealD * __restrict__ l,const RealD *__restrict__ r){ *y = (*l) + (*r);}
    inline RealD adj(const RealD & r){ return r; }  // No-op for real
    inline RealD conj(const RealD & r){ return r; }
    
    inline void mac (RealF * __restrict__ y,const RealF * __restrict__ a,const RealF *__restrict__ x){  *y = (*a) * (*x)+(*y); }
    inline void mult(RealF * __restrict__ y,const RealF * __restrict__ l,const RealF *__restrict__ r){ *y = (*l) * (*r); }
    inline void sub (RealF * __restrict__ y,const RealF * __restrict__ l,const RealF *__restrict__ r){ *y = (*l) - (*r); }
    inline void add (RealF * __restrict__ y,const RealF * __restrict__ l,const RealF *__restrict__ r){ *y = (*l) + (*r); }
    inline RealF adj(const RealF  & r){ return r; }
    inline RealF conj(const RealF  & r){ return r; }
    
    ////////////////////////////////////////////////////////////////////////
    //  Vector types are arch dependent/////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
#if defined (SSE2)
    typedef __m128 fvec;
    typedef __m128d dvec;
    typedef __m128 cvec;
    typedef __m128d zvec;
#endif
#if defined (AVX1) || defined (AVX2)
    typedef __m256 fvec;
    typedef __m256d dvec;
    typedef __m256 cvec;
    typedef __m256d zvec;
#endif
#if defined (AVX512)
    typedef __m512  fvec;
    typedef __m512d dvec;
    typedef __m512  cvec;
    typedef __m512d zvec;
#endif
#if defined (QPX)
    typedef float  fvec __attribute__ ((vector_size (16))); // QPX has same SIMD width irrespective of precision
    typedef float  cvec __attribute__ ((vector_size (16)));
    
    typedef vector4double dvec;
    typedef vector4double zvec;
#endif

#if defined (AVX1) || defined (AVX2) || defined (AVX512)
    inline void v_prefetch0(int size, const char *ptr){
          for(int i=0;i<size;i+=64){
            _mm_prefetch(ptr+i+4096,_MM_HINT_T1);
            _mm_prefetch(ptr+i+512,_MM_HINT_T0);
          }
    }
#endif
/*
      typedef  vComplexF vFComplex;
      typedef  vComplexD vDComplex;
     typedef  vComplexF vComplex;
    
    void zeroit(vRealF &z){ vzero(z);}
    void zeroit(vRealD &z){ vzero(z);}
    void zeroit(vComplexF &z){ vzero(z);}
    void zeroit(vComplexD &z){ vzero(z);}
    inline void zeroit(float &z){ z=0;}
    inline void zeroit(double &z){ z=0;}
    inline void zeroit(ComplexF &z){ z=0;}
    inline void zeroit(ComplexD &z){ z=0;}
*/

///////////////////////////////////////////////////
// Scalar, Vector, Matrix objects.
// These can be composed to form tensor products of internal indices.
///////////////////////////////////////////////////
    
template<class vtype> class iScalar
{
public:
  SIMDalign vtype _internal;
    iScalar(){};
    iScalar(Zero &z){ *this = zero; };
    iScalar<vtype> & operator= (const Zero &hero){
        zeroit(*this);
        return *this;
    }
    friend void zeroit(iScalar<vtype> &that){
        zeroit(that._internal);
    }
    // Unary negation
    friend inline iScalar<vtype> operator -(const iScalar<vtype> &r) {
        iScalar<vtype> ret;
        ret._internal= -r._internal;
        return ret;
    }
    // *=,+=,-= operators
    inline iScalar<vtype> &operator *=(const iScalar<vtype> &r) {
        *this = (*this)*r;
        return *this;
    }
    inline iScalar<vtype> &operator -=(const iScalar<vtype> &r) {
        *this = (*this)-r;
        return *this;
    }
    inline iScalar<vtype> &operator +=(const iScalar<vtype> &r) {
        *this = (*this)+r;
        return *this;
    }
    

};
    
template<class vtype,int N> class iVector
{
public:
  SIMDalign vtype _internal[N];
    iVector(Zero &z){ *this = zero; };
    iVector() {};
    iVector<vtype,N> & operator= (Zero &hero){
        zeroit(*this);
        return *this;
    }
    friend void zeroit(iVector<vtype,N> &that){
        for(int i=0;i<N;i++){
            zeroit(that._internal[i]);
        }
    }
    // Unary negation
    friend inline iVector<vtype,N> operator -(const iVector<vtype,N> &r) {
        iVector<vtype,N> ret;
        for(int i=0;i<N;i++) ret._internal[i]= -r._internal[i];
        return ret;
    }
    // *=,+=,-= operators
    inline iVector<vtype,N> &operator *=(const iScalar<vtype> &r) {
        *this = (*this)*r;
        return *this;
    }
    inline iVector<vtype,N> &operator -=(const iVector<vtype,N> &r) {
        *this = (*this)-r;
        return *this;
    }
    inline iVector<vtype,N> &operator +=(const iVector<vtype,N> &r) {
        *this = (*this)+r;
        return *this;
    }

};
    
    
template<class vtype,int N> class iMatrix
{
public:
  SIMDalign    vtype _internal[N][N];
    iMatrix(Zero &z){ *this = zero; };
    iMatrix() {};
    iMatrix<vtype,N> & operator= (Zero &hero){
        zeroit(*this);
        return *this;
    }
    friend void zeroit(iMatrix<vtype,N> &that){
        for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
                zeroit(that._internal[i][j]);
        }}
    }
    // Unary negation
    friend inline iMatrix<vtype,N> operator -(const iMatrix<vtype,N> &r) {
        iMatrix<vtype,N> ret;
        for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            ret._internal[i][j]= -r._internal[i][j];
        }}
        return ret;
    }
    // *=,+=,-= operators
    template<class T>
    inline iMatrix<vtype,N> &operator *=(const T &r) {
        *this = (*this)*r;
        return *this;
    }
    template<class T>
    inline iMatrix<vtype,N> &operator -=(const T &r) {
        *this = (*this)-r;
        return *this;
    }
    template<class T>
    inline iMatrix<vtype,N> &operator +=(const T &r) {
        *this = (*this)+r;
        return *this;
    }

};
/*
    inline vComplexD localInnerProduct(const vComplexD & l, const vComplexD & r) { return conj(l)*r; }
    inline vComplexF localInnerProduct(const vComplexF & l, const vComplexF & r) { return conj(l)*r; }
    inline vRealD localInnerProduct(const vRealD & l, const vRealD & r) { return conj(l)*r; }
    inline vRealF localInnerProduct(const vRealF & l, const vRealF & r) { return conj(l)*r; }
*/
    inline ComplexD localInnerProduct(const ComplexD & l, const ComplexD & r) { return conj(l)*r; }
    inline ComplexF localInnerProduct(const ComplexF & l, const ComplexF & r) { return conj(l)*r; }
    inline RealD localInnerProduct(const RealD & l, const RealD & r) { return conj(l)*r; }
    inline RealF localInnerProduct(const RealF & l, const RealF & r) { return conj(l)*r; }

    
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////// ADD         ///////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    

// ADD is simple for now; cannot mix types and straightforward template
// Scalar +/- Scalar
// Vector +/- Vector
// Matrix +/- Matrix
template<class vtype,class ltype,class rtype> inline void add(iScalar<vtype> * __restrict__ ret,
                                                              const iScalar<ltype> * __restrict__ lhs,
                                                              const iScalar<rtype> * __restrict__ rhs)
{
    add(&ret->_internal,&lhs->_internal,&rhs->_internal);
}
template<class vtype,class ltype,class rtype,int N> inline void add(iVector<vtype,N> * __restrict__ ret,
                                                                    const iVector<ltype,N> * __restrict__ lhs,
                                                                    const iVector<rtype,N> * __restrict__ rhs)
{
    for(int c=0;c<N;c++){
        ret->_internal[c]=lhs->_internal[c]+rhs->_internal[c];
    }
    return;
}
template<class vtype,class ltype,class rtype, int N> inline  void add(iMatrix<vtype,N> * __restrict__ ret,
                                                                      const iMatrix<ltype,N> * __restrict__ lhs,
                                                                      const iMatrix<rtype,N> * __restrict__ rhs)
{
    for(int c2=0;c2<N;c2++){
    for(int c1=0;c1<N;c1++){
        add(&ret->_internal[c1][c2],&lhs->_internal[c1][c2],&rhs->_internal[c1][c2]);
    }}
    return;
}
template<class vtype,class ltype,class rtype, int N> inline  void add(iMatrix<vtype,N> * __restrict__ ret,
                                                                      const iScalar<ltype>   * __restrict__ lhs,
                                                                      const iMatrix<rtype,N> * __restrict__ rhs)
{
    for(int c2=0;c2<N;c2++){
    for(int c1=0;c1<N;c1++){
        add(&ret->_internal[c1][c2],&lhs->_internal,&rhs->_internal[c1][c2]);
    }}
    return;
}
template<class vtype,class ltype,class rtype, int N> inline  void add(iMatrix<vtype,N> * __restrict__ ret,
                                                                      const iMatrix<ltype,N> * __restrict__ lhs,
                                                                      const iScalar<rtype>   * __restrict__ rhs)
{
    for(int c2=0;c2<N;c2++){
    for(int c1=0;c1<N;c1++){
        if ( c1==c2)
            add(&ret->_internal[c1][c2],&lhs->_internal[c1][c2],&rhs->_internal);
        else
            ret->_internal[c1][c2]=lhs->_internal[c1][c2];
    }}
    return;
}
// Need to figure multi-precision.
template<class Mytype>  Mytype timesI(Mytype &r)
{
    iScalar<Complex> i;
    i._internal = Complex(0,1);
    return r*i;
}

                // + operator for scalar, vector, matrix
template<class ltype,class rtype>
//inline auto operator + (iScalar<ltype>& lhs,iScalar<rtype>&& rhs) -> iScalar<decltype(lhs._internal + rhs._internal)>
inline auto operator + (const iScalar<ltype>& lhs,const iScalar<rtype>& rhs) -> iScalar<decltype(lhs._internal + rhs._internal)>
{
    typedef iScalar<decltype(lhs._internal+rhs._internal)> ret_t;
    ret_t ret;
    add(&ret,&lhs,&rhs);
    return ret;
}
template<class ltype,class rtype,int N>
inline auto operator + (const iVector<ltype,N>& lhs,const iVector<rtype,N>& rhs) ->iVector<decltype(lhs._internal[0]+rhs._internal[0]),N>
{
    typedef iVector<decltype(lhs._internal[0]+rhs._internal[0]),N> ret_t;
    ret_t ret;
    add(&ret,&lhs,&rhs);
    return ret;
}
template<class ltype,class rtype,int N>
inline auto operator + (const iMatrix<ltype,N>& lhs,const iMatrix<rtype,N>& rhs) ->iMatrix<decltype(lhs._internal[0][0]+rhs._internal[0][0]),N>
{
    typedef iMatrix<decltype(lhs._internal[0][0]+rhs._internal[0][0]),N> ret_t;
    ret_t ret;
    add(&ret,&lhs,&rhs);
    return ret;
}
template<class ltype,class rtype,int N>
inline auto operator + (const iScalar<ltype>& lhs,const iMatrix<rtype,N>& rhs)->iMatrix<decltype(lhs._internal+rhs._internal[0][0]),N>
{
    typedef iMatrix<decltype(lhs._internal+rhs._internal[0][0]),N> ret_t;
    ret_t ret;
    add(&ret,&lhs,&rhs);
    return ret;
}
template<class ltype,class rtype,int N>
inline auto operator + (const iMatrix<ltype,N>& lhs,const iScalar<rtype>& rhs)->iMatrix<decltype(lhs._internal[0][0]+rhs._internal),N>
{
    typedef iMatrix<decltype(lhs._internal[0][0]+rhs._internal),N> ret_t;
    ret_t ret;
    add(&ret,&lhs,&rhs);
    return ret;
}


    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////// SUB         ///////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    

// SUB is simple for now; cannot mix types and straightforward template
// Scalar +/- Scalar
// Vector +/- Vector
// Matrix +/- Matrix
// Matrix /- scalar
template<class vtype,class ltype,class rtype> inline void sub(iScalar<vtype> * __restrict__ ret,
                                                              const iScalar<ltype> * __restrict__ lhs,
                                                              const iScalar<rtype> * __restrict__ rhs)
{
    sub(&ret->_internal,&lhs->_internal,&rhs->_internal);
}

template<class vtype,class ltype,class rtype,int N> inline void sub(iVector<vtype,N> * __restrict__ ret,
                                                                    const iVector<ltype,N> * __restrict__ lhs,
                                                                    const iVector<rtype,N> * __restrict__ rhs)
{
    for(int c=0;c<N;c++){
        ret->_internal[c]=lhs->_internal[c]-rhs->_internal[c];
    }
    return;
}
template<class vtype,class ltype,class rtype, int N> inline void sub(iMatrix<vtype,N> * __restrict__ ret,
                                                                     const iMatrix<ltype,N> * __restrict__ lhs,
                                                                     const iMatrix<rtype,N> * __restrict__ rhs){
    for(int c2=0;c2<N;c2++){
    for(int c1=0;c1<N;c1++){
        sub(&ret->_internal[c1][c2],&lhs->_internal[c1][c2],&rhs->_internal[c1][c2]);
    }}
    return;
}
template<class vtype,class ltype,class rtype, int N> inline void sub(iMatrix<vtype,N> * __restrict__ ret,
                                                                     const iScalar<ltype> * __restrict__ lhs,
                                                                     const iMatrix<rtype,N> * __restrict__ rhs){
    for(int c2=0;c2<N;c2++){
    for(int c1=0;c1<N;c1++){
        if ( c1!=c2) {
            sub(&ret->_internal[c1][c2],&lhs->_internal,&rhs->_internal[c1][c2]);
        } else {
            // Fails -- need unary minus. Catalogue other unops?
            ret->_internal[c1][c2]=zero;
            ret->_internal[c1][c2]=ret->_internal[c1][c2]-rhs->_internal[c1][c2];

        }
    }}
    return;
}
template<class vtype,class ltype,class rtype, int N> inline void sub(iMatrix<vtype,N> * __restrict__ ret,
                                                                     const iMatrix<ltype,N> * __restrict__ lhs,
                                                                     const iScalar<rtype> * __restrict__ rhs){
    for(int c2=0;c2<N;c2++){
    for(int c1=0;c1<N;c1++){
        if ( c1!=c2)
            sub(&ret->_internal[c1][c2],&lhs->_internal[c1][c2],&rhs->_internal);
        else
            ret->_internal[c1][c2]=lhs->_internal[c1][c2];
    }}
    return;
}

template<class v> void vprefetch(const iScalar<v> &vv)
{
  vprefetch(vv._internal);
}
template<class v,int N> void vprefetch(const iVector<v,N> &vv)
{
  for(int i=0;i<N;i++){
    vprefetch(vv._internal[i]);
  }
}
template<class v,int N> void vprefetch(const iMatrix<v,N> &vv)
{
  for(int i=0;i<N;i++){
  for(int j=0;j<N;j++){
    vprefetch(vv._internal[i][j]);
  }}
}

    // - operator for scalar, vector, matrix
template<class ltype,class rtype> inline auto
operator - (const iScalar<ltype>& lhs, const iScalar<rtype>& rhs) -> iScalar<decltype(lhs._internal - rhs._internal)>
{
    typedef iScalar<decltype(lhs._internal-rhs._internal)> ret_t;
    ret_t ret;
    sub(&ret,&lhs,&rhs);
    return ret;
}
template<class ltype,class rtype,int N>
inline auto operator - (const iVector<ltype,N>& lhs,const iVector<rtype,N>& rhs) ->iVector<decltype(lhs._internal[0]-rhs._internal[0]),N>
{
    typedef iVector<decltype(lhs._internal[0]-rhs._internal[0]),N> ret_t;
    ret_t ret;
    sub(&ret,&lhs,&rhs);
    return ret;
}
template<class ltype,class rtype,int N>
inline auto operator - (const iMatrix<ltype,N>& lhs,const iMatrix<rtype,N>& rhs) ->iMatrix<decltype(lhs._internal[0][0]-rhs._internal[0][0]),N>
{
    typedef iMatrix<decltype(lhs._internal[0][0]-rhs._internal[0][0]),N> ret_t;
    ret_t ret;
    sub(&ret,&lhs,&rhs);
    return ret;
}
template<class ltype,class rtype,int N>
inline auto operator - (const iScalar<ltype>& lhs,const iMatrix<rtype,N>& rhs)->iMatrix<decltype(lhs._internal-rhs._internal[0][0]),N>
{
    typedef iMatrix<decltype(lhs._internal-rhs._internal[0][0]),N> ret_t;
    ret_t ret;
    sub(&ret,&lhs,&rhs);
    return ret;
}
template<class ltype,class rtype,int N>
inline auto operator - (const iMatrix<ltype,N>& lhs,const iScalar<rtype>& rhs)->iMatrix<decltype(lhs._internal[0][0]-rhs._internal),N>
{
    typedef iMatrix<decltype(lhs._internal[0][0]-rhs._internal),N> ret_t;
    ret_t ret;
    sub(&ret,&lhs,&rhs);
    return ret;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////// MAC         ///////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

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
    ///////////////////////////////////////////////////////////////////////////////////////
    // localInnerProduct Scalar x Scalar -> Scalar
    // localInnerProduct Vector x Vector -> Scalar
    // localInnerProduct Matrix x Matrix -> Scalar
    ///////////////////////////////////////////////////////////////////////////////////////
    template<class l,class r,int N> inline
    auto localInnerProduct (const iVector<l,N>& lhs,const iVector<r,N>& rhs) -> iScalar<decltype(localInnerProduct(lhs._internal[0],rhs._internal[0]))>
    {
        typedef decltype(localInnerProduct(lhs._internal[0],rhs._internal[0])) ret_t;
        iScalar<ret_t> ret=zero;
        for(int c1=0;c1<N;c1++){
            ret._internal += localInnerProduct(lhs._internal[c1],rhs._internal[c1]);
        }
        return ret;
    }
    template<class l,class r,int N> inline
    auto localInnerProduct (const iMatrix<l,N>& lhs,const iMatrix<r,N>& rhs) -> iScalar<decltype(localInnerProduct(lhs._internal[0][0],rhs._internal[0][0]))>
    {
        typedef decltype(localInnerProduct(lhs._internal[0][0],rhs._internal[0][0])) ret_t;
        iScalar<ret_t> ret=zero;
        for(int c1=0;c1<N;c1++){
        for(int c2=0;c2<N;c2++){
            ret._internal += localInnerProduct(lhs._internal[c1][c2],rhs._internal[c1][c2]);
        }}
        return ret;
    }
    template<class l,class r> inline
    auto localInnerProduct (const iScalar<l>& lhs,const iScalar<r>& rhs) -> iScalar<decltype(localInnerProduct(lhs._internal,rhs._internal))>
    {
        typedef decltype(localInnerProduct(lhs._internal,rhs._internal)) ret_t;
        iScalar<ret_t> ret;
        ret._internal = localInnerProduct(lhs._internal,rhs._internal);
        return ret;
    }

    ///////////////////////////////////////////////////////////////////////////////////////
    // outerProduct Scalar x Scalar -> Scalar
    //              Vector x Vector -> Matrix
    ///////////////////////////////////////////////////////////////////////////////////////

template<class l,class r,int N> inline
auto outerProduct (const iVector<l,N>& lhs,const iVector<r,N>& rhs) -> iMatrix<decltype(outerProduct(lhs._internal[0],rhs._internal[0])),N>
{
    typedef decltype(outerProduct(lhs._internal[0],rhs._internal[0])) ret_t;
    iMatrix<ret_t,N> ret;
    for(int c1=0;c1<N;c1++){
    for(int c2=0;c2<N;c2++){
        ret._internal[c1][c2] = outerProduct(lhs._internal[c1],rhs._internal[c2]);
    }}
    return ret;
}
template<class l,class r> inline
auto outerProduct (const iScalar<l>& lhs,const iScalar<r>& rhs) -> iScalar<decltype(outerProduct(lhs._internal,rhs._internal))>
{
    typedef decltype(outerProduct(lhs._internal,rhs._internal)) ret_t;
    iScalar<ret_t> ret;
    ret._internal = outerProduct(lhs._internal,rhs._internal);
    return ret;
}
/*
    inline vComplexF outerProduct(const vComplexF &l, const vComplexF& r)
    {
        return l*r;
    }
    inline vComplexD outerProduct(const vComplexD &l, const vComplexD& r)
    {
        return l*r;
    }
    inline vRealF outerProduct(const vRealF &l, const vRealF& r)
    {
        return l*r;
    }
    inline vRealD outerProduct(const vRealD &l, const vRealD& r)
    {
        return l*r;
    }
*/
    inline ComplexF outerProduct(const ComplexF &l, const ComplexF& r)
    {
        return l*r;
    }
    inline ComplexD outerProduct(const ComplexD &l, const ComplexD& r)
    {
        return l*r;
    }
    inline RealF outerProduct(const RealF &l, const RealF& r)
    {
        return l*r;
    }
    inline RealD outerProduct(const RealD &l, const RealD& r)
    {
        return l*r;
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////// CONJ         ///////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
 
// Conj function for scalar, vector, matrix
template<class vtype> inline iScalar<vtype> conj(const iScalar<vtype>&r)
{
    iScalar<vtype> ret;
    ret._internal = conj(r._internal);
    return ret;
}

// Adj function for scalar, vector, matrix
template<class vtype> inline iScalar<vtype> adj(const iScalar<vtype>&r)
{
    iScalar<vtype> ret;
    ret._internal = adj(r._internal);
    return ret;
}
template<class vtype,int N> inline iVector<vtype,N> adj(const iVector<vtype,N>&r)
{
    iVector<vtype,N> ret;
    for(int i=0;i<N;i++){
        ret._internal[i] = adj(r._internal[i]);
    }
    return ret;
}
template<class vtype,int N> inline iMatrix<vtype,N> adj(const iMatrix<vtype,N> &arg)
{
    iMatrix<vtype,N> ret;
    for(int c1=0;c1<N;c1++){
    for(int c2=0;c2<N;c2++){
        ret._internal[c1][c2]=adj(arg._internal[c2][c1]);
    }}
    return ret;
}

/////////////////////////////////////////////////////////////////
// Can only take the real/imag part of scalar objects, since
// lattice objects of different complexity are non-conformable.
/////////////////////////////////////////////////////////////////
template<class itype> inline auto real(const iScalar<itype> &z) -> iScalar<decltype(real(z._internal))>
{
    iScalar<decltype(real(z._internal))> ret;
    ret._internal = real(z._internal);
    return ret;
}
template<class itype,int N> inline auto real(const iMatrix<itype,N> &z) -> iMatrix<decltype(real(z._internal[0][0])),N>
{
    iMatrix<decltype(real(z._internal[0][0])),N> ret;
    for(int c1=0;c1<N;c1++){
    for(int c2=0;c2<N;c2++){
        ret._internal[c1][c2] = real(z._internal[c1][c2]);
    }}
    return ret;
}
template<class itype,int N> inline auto real(const iVector<itype,N> &z) -> iVector<decltype(real(z._internal[0])),N>
{
    iVector<decltype(real(z._internal[0])),N> ret;
    for(int c1=0;c1<N;c1++){
        ret._internal[c1] = real(z._internal[c1]);
    }
    return ret;
}
    
template<class itype> inline auto imag(const iScalar<itype> &z) -> iScalar<decltype(imag(z._internal))>
{
    iScalar<decltype(imag(z._internal))> ret;
    ret._internal = imag(z._internal);
    return ret;
}
template<class itype,int N> inline auto imag(const iMatrix<itype,N> &z) -> iMatrix<decltype(imag(z._internal[0][0])),N>
{
    iMatrix<decltype(imag(z._internal[0][0])),N> ret;
    for(int c1=0;c1<N;c1++){
    for(int c2=0;c2<N;c2++){
        ret._internal[c1][c2] = imag(z._internal[c1][c2]);
    }}
    return ret;
}
template<class itype,int N> inline auto imag(const iVector<itype,N> &z) -> iVector<decltype(imag(z._internal[0])),N>
{
    iVector<decltype(imag(z._internal[0])),N> ret;
    for(int c1=0;c1<N;c1++){
        ret._internal[c1] = imag(z._internal[c1]);
    }
    return ret;
}

    /////////////////////////////////
    // Trace of scalar and matrix
    /////////////////////////////////

inline Complex trace( const Complex &arg){
    return arg;
}
//inline vComplex trace(const vComplex &arg){
//    return arg;
//}
template<class vtype,int N>
inline auto trace(const iMatrix<vtype,N> &arg) -> iScalar<decltype(trace(arg._internal[0][0]))>
{
    iScalar<decltype( trace(arg._internal[0][0] )) > ret;
    ZeroIt(ret._internal);
    for(int i=0;i<N;i++){
        ret._internal=ret._internal+trace(arg._internal[i][i]);
    }
    return ret;
}
template<class vtype>
inline auto trace(const iScalar<vtype> &arg) -> iScalar<decltype(trace(arg._internal))>
{
    iScalar<decltype(trace(arg._internal))> ret;
    ret._internal=trace(arg._internal);
    return ret;
}
    
/////////////////////////////////////////////////////////////////////////
// Generic routine to promote object<complex> -> object<vcomplex>
// Supports the array reordering transformation that gives me SIMD utilisation
/////////////////////////////////////////////////////////////////////////
/*
template<template<class> class object>
inline object<vComplex> splat(object<Complex >s){
    object<vComplex> ret;
    vComplex * v_ptr = (vComplex *)& ret;
    Complex * s_ptr = (Complex *) &s;
    for(int i=0;i<sizeof(ret);i+=sizeof(vComplex)){
        vsplat(*(v_ptr++),*(s_ptr++));
    }
    return ret;
}
*/
    
    /////////////////////////////////////////////////////////////////////////////////////////
    // Grid Support. Following will go into Grid.h.
    /////////////////////////////////////////////////////////////////////////////////////////
    // Cartesian grids
    // dpo::Grid
    // dpo::GridCartesian
    // dpo::GridCartesianRedBlack
    
class Grid {
public:
    // Give Lattice access
    template<class object> friend class Lattice;

        
//protected:
        
    // Lattice wide random support. not yet fully implemented. Need seed strategy
    // and one generator per site.
    //std::default_random_engine generator;
  //    static std::mt19937  generator( 9 );

        
    // Grid information.
    unsigned long _ndimension;
    std::vector<int> _layout;     // Which dimensions get relayed out over simd lanes.
    std::vector<int> _dimensions; // Dimensions of array
    std::vector<int> _rdimensions;// Reduced dimensions with simd lane images removed
    std::vector<int> _ostride;    // Outer stride for each dimension
    std::vector<int> _istride;    // Inner stride i.e. within simd lane
    int _osites;                  // _isites*_osites = product(dimensions).
    int _isites;
        
    // subslice information
    std::vector<int> _slice_block;
    std::vector<int> _slice_stride;
    std::vector<int> _slice_nblock;
public:
    
    // These routines are key. Subdivide the linearised cartesian index into
    //      "inner" index identifying which simd lane of object<vFcomplex> is associated with coord
    //      "outer" index identifying which element of _odata in class "Lattice" is associated with coord.
    // Compared to, say, Blitz++ we simply need to store BOTH an inner stride and an outer
    // stride per dimension. The cost of evaluating the indexing information is doubled for an n-dimensional
    // coordinate. Note, however, for data parallel operations the "inner" indexing cost is not paid and all
    // lanes are operated upon simultaneously.
    
    inline int oIndexReduced(std::vector<int> &rcoor)
    {
        int idx=0;
        for(int d=0;d<_ndimension;d++) idx+=_ostride[d]*rcoor[d];
        return idx;
    }
    virtual int oIndex(std::vector<int> &coor)
    {
        int idx=0;
        for(int d=0;d<_ndimension;d++) idx+=_ostride[d]*(coor[d]%_rdimensions[d]);
        return idx;
    }
    inline int iIndex(std::vector<int> &rcoor)
    {
        int idx=0;
        for(int d=0;d<_ndimension;d++) idx+=_istride[d]*(rcoor[d]/_rdimensions[d]);
        return idx;
    }
        
    inline int oSites(void) { return _osites; };
    inline int iSites(void) { return _isites; };
    virtual int CheckerBoard(std::vector<int> site)=0;
    virtual int CheckerBoardDestination(int source_cb,int shift)=0;
    virtual int CheckerBoardShift(int source_cb,int dim,int shift)=0;
};

////////////////////////////////////////////////////////////////////
// A lattice of something, but assume the something is SIMDized.
////////////////////////////////////////////////////////////////////
template<typename _Tp>
class myallocator {
public: 
  typedef std::size_t     size_type;
  typedef std::ptrdiff_t  difference_type;
  typedef _Tp*       pointer;
  typedef const _Tp* const_pointer;
  typedef _Tp&       reference;
  typedef const _Tp& const_reference;
  typedef _Tp        value_type;

  template<typename _Tp1>  struct rebind { typedef myallocator<_Tp1> other; };
  myallocator() throw() { }
  myallocator(const myallocator&) throw() { }
  template<typename _Tp1> myallocator(const myallocator<_Tp1>&) throw() { }
  ~myallocator() throw() { }
  pointer address(reference __x) const { return &__x; }
  const_pointer address(const_reference __x) const { return &__x; }
  size_type  max_size() const throw() { return size_t(-1) / sizeof(_Tp); }
  // Should override allocate and deallocate
  pointer allocate(size_type __n, const void* = 0)
  { 
    //_Tp * ptr = (_Tp *) memalign(sizeof(_Tp),__n*sizeof(_Tp));
    // _Tp * ptr = (_Tp *) memalign(128,__n*sizeof(_Tp));
#ifdef AVX512
    _Tp * ptr = (_Tp *) memalign(128,__n*sizeof(_Tp));
#else
    _Tp * ptr = (_Tp *) _mm_malloc(__n*sizeof(_Tp),128);
#endif

    return ptr;
  }
  void deallocate(pointer __p, size_type) { 
    free(__p); 
  }
  void construct(pointer __p, const _Tp& __val) { };
  void construct(pointer __p) { };
  void destroy(pointer __p) { };
};

template<typename _Tp>  inline bool
operator==(const myallocator<_Tp>&, const myallocator<_Tp>&){ return true; }

template<typename _Tp>  inline bool
operator!=(const myallocator<_Tp>&, const myallocator<_Tp>&){ return false; }

    
}; // namespace dpo
    
////////////////////////////////////////////////////////////////////////////////////////
// Test code
////////////////////////////////////////////////////////////////////////////////////////
/*
using namespace std;
using namespace dpo;
using namespace dpo::QCD;
*/


#endif
