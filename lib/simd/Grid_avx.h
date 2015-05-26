//----------------------------------------------------------------------
/*! @file Grid_avx.h
  @brief Optimization libraries for AVX1/2 instructions set

  Using intrinsics
*/
// Time-stamp: <2015-05-22 18:58:27 neo>
//----------------------------------------------------------------------

#include <immintrin.h>
// _mm256_set_m128i(hi,lo); // not defined in all versions of immintrin.h
#ifndef _mm256_set_m128i
#define _mm256_set_m128i(hi,lo) _mm256_insertf128_si256(_mm256_castsi128_si256(lo),(hi),1)
#endif

namespace Optimization {
  
  struct Vsplat{
    //Complex float
    inline __m256 operator()(float a, float b){
      return _mm256_set_ps(b,a,b,a,b,a,b,a);
    }
    // Real float
    inline __m256 operator()(float a){
      return _mm256_set_ps(a,a,a,a,a,a,a,a);
    }
    //Complex double
    inline __m256d operator()(double a, double b){
      return _mm256_set_pd(b,a,b,a);
    }
    //Real double
    inline __m256d operator()(double a){
      return _mm256_set_pd(a,a,a,a);
    }
    //Integer
    inline __m256i operator()(Integer a){
      return _mm256_set1_epi32(a);
    }
  };

  struct Vstore{
    //Float 
    inline void operator()(__m256 a, float* F){
      _mm256_store_ps(F,a);
    }
    //Double
    inline void operator()(__m256d a, double* D){
      _mm256_store_pd(D,a);
    }
    //Integer
    inline void operator()(__m256i a, Integer* I){
      _mm256_store_si256((__m256i*)I,a);
    }

  };


  struct Vstream{
    //Float
    inline void operator()(float * a, __m256 b){
      _mm256_stream_ps(a,b);
    }
    //Double
    inline void operator()(double * a, __m256d b){
      _mm256_stream_pd(a,b);
    }


  };



  struct Vset{
    // Complex float 
    inline __m256 operator()(Grid::ComplexF *a){
      return _mm256_set_ps(a[3].imag(),a[3].real(),a[2].imag(),a[2].real(),a[1].imag(),a[1].real(),a[0].imag(),a[0].real());
    }
    // Complex double 
    inline __m256d operator()(Grid::ComplexD *a){
      return _mm256_set_pd(a[1].imag(),a[1].real(),a[0].imag(),a[0].real());
    }
    // Real float 
    inline __m256 operator()(float *a){
      return _mm256_set_ps(a[7],a[6],a[5],a[4],a[3],a[2],a[1],a[0]);
    }
    // Real double
    inline __m256d operator()(double *a){
      return _mm256_set_pd(a[3],a[2],a[1],a[0]);
    }
    // Integer
    inline __m256i operator()(Integer *a){
      return _mm256_set_epi32(a[7],a[6],a[5],a[4],a[3],a[2],a[1],a[0]);
    }


  };

  template <typename Out_type, typename In_type>
  struct Reduce{
    //Need templated class to overload output type
    //General form must generate error if compiled
    inline Out_type operator()(In_type in){
      printf("Error, using wrong Reduce function\n");
      exit(1);
      return 0;
    }
  };


 

  /////////////////////////////////////////////////////
  // Arithmetic operations
  /////////////////////////////////////////////////////
  struct Sum{
    //Complex/Real float
    inline __m256 operator()(__m256 a, __m256 b){
      return _mm256_add_ps(a,b);
    }
    //Complex/Real double
    inline __m256d operator()(__m256d a, __m256d b){
      return _mm256_add_pd(a,b);
    }
    //Integer
    inline __m256i operator()(__m256i a, __m256i b){
#if defined (AVX1) 
          __m128i a0,a1;
          __m128i b0,b1;
          a0 = _mm256_extractf128_si256(a,0);
          b0 = _mm256_extractf128_si256(b,0);
          a1 = _mm256_extractf128_si256(a,1);
          b1 = _mm256_extractf128_si256(b,1);
          a0 = _mm_add_epi32(a0,b0);
          a1 = _mm_add_epi32(a1,b1);
          return _mm256_set_m128i(a1,a0);
#endif
#if defined (AVX2)
            return _mm256_add_epi32(a,b);
#endif

    }
  };

  struct Sub{
    //Complex/Real float
    inline __m256 operator()(__m256 a, __m256 b){
      return _mm256_sub_ps(a,b);
    }
    //Complex/Real double
    inline __m256d operator()(__m256d a, __m256d b){
      return _mm256_sub_pd(a,b);
    }
    //Integer
    inline __m256i operator()(__m256i a, __m256i b){
#if defined (AVX1) 
          __m128i a0,a1;
          __m128i b0,b1;
          a0 = _mm256_extractf128_si256(a,0);
          b0 = _mm256_extractf128_si256(b,0);
          a1 = _mm256_extractf128_si256(a,1);
          b1 = _mm256_extractf128_si256(b,1);
          a0 = _mm_sub_epi32(a0,b0);
          a1 = _mm_sub_epi32(a1,b1);
          return _mm256_set_m128i(a1,a0);
#endif
#if defined (AVX2)
            return _mm256_sub_epi32(a,b);
#endif

    }
  };


  struct MultComplex{
    // Complex float
    inline __m256 operator()(__m256 a, __m256 b){
      __m256 ymm0,ymm1,ymm2;
      ymm0 = _mm256_shuffle_ps(a,a,_MM_SHUFFLE(2,2,0,0)); // ymm0 <- ar ar,
      ymm0 = _mm256_mul_ps(ymm0,b);                       // ymm0 <- ar bi, ar br
      // FIXME AVX2 could MAC
      ymm1 = _mm256_shuffle_ps(b,b,_MM_SHUFFLE(2,3,0,1)); // ymm1 <- br,bi
      ymm2 = _mm256_shuffle_ps(a,a,_MM_SHUFFLE(3,3,1,1)); // ymm2 <- ai,ai
      ymm1 = _mm256_mul_ps(ymm1,ymm2);                    // ymm1 <- br ai, ai bi
      return _mm256_addsub_ps(ymm0,ymm1);  
    }
    // Complex double
    inline __m256d operator()(__m256d a, __m256d b){
      //Multiplication of (ak+ibk)*(ck+idk)
      // a + i b can be stored as a data structure
      //From intel optimisation reference guide
      /*
	movsldup xmm0, Src1; load real parts into the destination,
	; a1, a1, a0, a0
	movaps xmm1, src2; load the 2nd pair of complex values, ; i.e. d1, c1, d0, c0
	mulps xmm0, xmm1; temporary results, a1d1, a1c1, a0d0, ; a0c0
	shufps xmm1, xmm1, b1; reorder the real and imaginary ; parts, c1, d1, c0, d0
	movshdup xmm2, Src1; load the imaginary parts into the ; destination, b1, b1, b0, b0
	mulps xmm2, xmm1; temporary results, b1c1, b1d1, b0c0, ; b0d0
	addsubps xmm0, xmm2; b1c1+a1d1, a1c1 -b1d1, b0c0+a0d
	VSHUFPD (VEX.256 encoded version)
	IF IMM0[0] = 0
	THEN DEST[63:0]=SRC1[63:0] ELSE DEST[63:0]=SRC1[127:64] FI;
	IF IMM0[1] = 0
	THEN DEST[127:64]=SRC2[63:0] ELSE DEST[127:64]=SRC2[127:64] FI;
	IF IMM0[2] = 0
	THEN DEST[191:128]=SRC1[191:128] ELSE DEST[191:128]=SRC1[255:192] FI;
	IF IMM0[3] = 0
	THEN DEST[255:192]=SRC2[191:128] ELSE DEST[255:192]=SRC2[255:192] FI; // Ox5 r<->i   ; 0xC unchanged
      */
      
      __m256d ymm0,ymm1,ymm2;
      ymm0 = _mm256_shuffle_pd(a,a,0x0); // ymm0 <- ar ar, ar,ar b'00,00
      ymm0 = _mm256_mul_pd(ymm0,b);      // ymm0 <- ar bi, ar br
      ymm1 = _mm256_shuffle_pd(b,b,0x5); // ymm1 <- br,bi  b'01,01 
      ymm2 = _mm256_shuffle_pd(a,a,0xF); // ymm2 <- ai,ai  b'11,11
      ymm1 = _mm256_mul_pd(ymm1,ymm2);   // ymm1 <- br ai, ai bi
      return _mm256_addsub_pd(ymm0,ymm1);
    }
  };

  struct Mult{
    // Real float
    inline __m256 operator()(__m256 a, __m256 b){
      return _mm256_mul_ps(a,b);
    }
    // Real double
    inline __m256d operator()(__m256d a, __m256d b){
      return _mm256_mul_pd(a,b);
    }
    // Integer
    inline __m256i operator()(__m256i a, __m256i b){
#if defined (AVX1) 
      __m128i a0,a1;
      __m128i b0,b1;
      a0 = _mm256_extractf128_si256(a,0);
      b0 = _mm256_extractf128_si256(b,0);
      a1 = _mm256_extractf128_si256(a,1);
      b1 = _mm256_extractf128_si256(b,1);
      a0 = _mm_mul_epi32(a0,b0);
      a1 = _mm_mul_epi32(a1,b1);
      return _mm256_set_m128i(a1,a0);
#endif
#if defined (AVX2)
      return _mm256_mul_epi32(a,b);
#endif

    }
  };


  struct Conj{
    // Complex single
    inline __m256 operator()(__m256 in){
      return _mm256_xor_ps(_mm256_addsub_ps(_mm256_setzero_ps(),in), _mm256_set1_ps(-0.f));
    }
    // Complex double
    inline __m256d operator()(__m256d in){
      return _mm256_xor_pd(_mm256_addsub_pd(_mm256_setzero_pd(),in), _mm256_set1_pd(-0.f));//untested
      /*
	// original 
	//      addsubps 0, inv=>0+in.v[3] 0-in.v[2], 0+in.v[1], 0-in.v[0], ...
	__m256d tmp = _mm256_addsub_pd(_mm256_setzero_pd(),_mm256_shuffle_pd(in,in,0x5));
	return _mm256_shuffle_pd(tmp,tmp,0x5);
      */
    }
    // do not define for integer input
  };

  struct TimesMinusI{
    //Complex single
    inline __m256 operator()(__m256 in, __m256 ret){
      __m256 tmp =_mm256_addsub_ps(_mm256_setzero_ps(),in);   // r,-i
      return _mm256_shuffle_ps(tmp,tmp,_MM_SHUFFLE(2,3,0,1)); //-i,r
    }
    //Complex double
    inline __m256d operator()(__m256d in, __m256d ret){
      __m256d tmp = _mm256_addsub_pd(_mm256_setzero_pd(),in); // r,-i
      return _mm256_shuffle_pd(tmp,tmp,0x5);
    }
  };

  struct TimesI{
    //Complex single
    inline __m256 operator()(__m256 in, __m256 ret){
      __m256 tmp =_mm256_shuffle_ps(in,in,_MM_SHUFFLE(2,3,0,1)); // i,r
      return _mm256_addsub_ps(_mm256_setzero_ps(),tmp);          // i,-r
    }
    //Complex double
    inline __m256d operator()(__m256d in, __m256d ret){
      __m256d tmp = _mm256_shuffle_pd(in,in,0x5);
      return _mm256_addsub_pd(_mm256_setzero_pd(),tmp); // i,-r
    }
  };


  


  //////////////////////////////////////////////
  // Some Template specialization
  template < typename vtype > 
    void permute(vtype a, vtype b, int perm) {
    union { 
      __m256 f;
      vtype v;
    } conv;
    conv.v = b;
    switch (perm){
      // 8x32 bits=>3 permutes
    case 2: conv.f = _mm256_shuffle_ps(conv.f,conv.f,_MM_SHUFFLE(2,3,0,1)); break;
    case 1: conv.f = _mm256_shuffle_ps(conv.f,conv.f,_MM_SHUFFLE(1,0,3,2)); break;
    case 0: conv.f = _mm256_permute2f128_ps(conv.f,conv.f,0x01); break;
    default: assert(0); break;
    }
    a = conv.v;
    
  }
  
  //Complex float Reduce
  template<>
    inline Grid::ComplexF Reduce<Grid::ComplexF, __m256>::operator()(__m256 in){
    __m256 v1,v2;
    Optimization::permute(v1,in,0); // sse 128; paired complex single
    v1 = _mm256_add_ps(v1,in);
    Optimization::permute(v2,v1,1); // avx 256; quad complex single
    v1 = _mm256_add_ps(v1,v2);
    return Grid::ComplexF(v1[0],v1[1]);
  }
  //Real float Reduce
  template<>
  inline Grid::RealF Reduce<Grid::RealF, __m256>::operator()(__m256 in){
    __m256 v1,v2;
    Optimization::permute(v1,in,0); // avx 256; octo-double
    v1 = _mm256_add_ps(v1,in);
    Optimization::permute(v2,v1,1); 
    v1 = _mm256_add_ps(v1,v2);
    Optimization::permute(v2,v1,2); 
    v1 = _mm256_add_ps(v1,v2);
    return v1[0];
  }
  
  
  //Complex double Reduce
  template<>
  inline Grid::ComplexD Reduce<Grid::ComplexD, __m256d>::operator()(__m256d in){
    __m256d v1;
    Optimization::permute(v1,in,0); // sse 128; paired complex single
    v1 = _mm256_add_pd(v1,in);
    return Grid::ComplexD(v1[0],v1[1]);
  }
  
  //Real double Reduce
  template<>
  inline Grid::RealD Reduce<Grid::RealD, __m256d>::operator()(__m256d in){
    __m256d v1,v2;
    Optimization::permute(v1,in,0); // avx 256; quad double
    v1 = _mm256_add_pd(v1,in);
    Optimization::permute(v2,v1,1); 
    v1 = _mm256_add_pd(v1,v2);
    return v1[0];
  }

  //Integer Reduce
  template<>
  inline Integer Reduce<Integer, __m256i>::operator()(__m256i in){
    // FIXME unimplemented
    printf("Reduce : Missing integer implementation -> FIX\n");
    assert(0);
  }
  
  
}

//////////////////////////////////////////////////////////////////////////////////////
// Here assign types 
namespace Grid {
  typedef __m256  SIMD_Ftype; // Single precision type
  typedef __m256d SIMD_Dtype; // Double precision type
  typedef __m256i SIMD_Itype; // Integer type


  // Function name aliases
  typedef Optimization::Vsplat   VsplatSIMD;
  typedef Optimization::Vstore   VstoreSIMD;
  typedef Optimization::Vset     VsetSIMD;
  typedef Optimization::Vstream  VstreamSIMD;
  template <typename S, typename T> using ReduceSIMD = Optimization::Reduce<S,T>;


  // Arithmetic operations
  typedef Optimization::Sum         SumSIMD;
  typedef Optimization::Sub         SubSIMD;
  typedef Optimization::Mult        MultSIMD;
  typedef Optimization::MultComplex MultComplexSIMD;
  typedef Optimization::Conj        ConjSIMD;
  typedef Optimization::TimesMinusI TimesMinusISIMD;
  typedef Optimization::TimesI      TimesISIMD;

}
