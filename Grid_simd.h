#ifndef GRID_SIMD_H
#define GRID_SIMD_H

////////////////////////////////////////////////////////////////////////
// Define scalar and vector floating point types
//
// Scalar:   RealF, RealD, ComplexF, ComplexD
//
// Vector:  vRealF, vRealD, vComplexF, vComplexD
//
// Vector types are arch dependent
////////////////////////////////////////////////////////////////////////

namespace dpo {

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


  ////////////////////////////////////////////////////////////
  // SIMD Alignment controls
  ////////////////////////////////////////////////////////////
#ifdef HAVE_VAR_ATTRIBUTE_ALIGNED
#define ALIGN_DIRECTIVE(A) __attribute__ ((aligned(A)))
#else
#define ALIGN_DIRECTIVE(A) __declspec(align(A))
#endif

#ifdef SSE2
#include <pmmintrin.h>
#define SIMDalign ALIGN_DIRECTIVE(16)
#endif

#if defined(AVX1) || defined (AVX2)
#include <immintrin.h>
#define SIMDalign ALIGN_DIRECTIVE(32)
#endif

#ifdef AVX512
#include <immintrin.h>
#define SIMDalign ALIGN_DIRECTIVE(64)
#endif

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
          for(int i=0;i<size;i+=64){ //  Define L1 linesize above// What about SSE?
            _mm_prefetch(ptr+i+4096,_MM_HINT_T1);
            _mm_prefetch(ptr+i+512,_MM_HINT_T0);
          }
    }
#else 
    inline void v_prefetch0(int size, const char *ptr){};
#endif

};

#include <Grid_vRealF.h>
#include <Grid_vRealD.h>
#include <Grid_vComplexF.h>
#include <Grid_vComplexD.h>


#endif
