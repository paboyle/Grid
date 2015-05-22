//----------------------------------------------------------------------
/*! @file Grid_sse4.h
  @brief Optimization libraries for SSE4 instructions set

  Using intrinsics
*/
// Time-stamp: <2015-05-21 18:06:30 neo>
//----------------------------------------------------------------------

#include <pmmintrin.h>

namespace Optimization {
  
  struct Vsplat{
    //Complex float
    inline __m128 operator()(float a, float b){
      return _mm_set_ps(b,a,b,a);
    }
    // Real float
    inline __m128 operator()(float a){
      return _mm_set_ps(a,a,a,a);
    }
    //Complex double
    inline __m128d operator()(double a, double b){
      return _mm_set_pd(b,a);
    }
    //Real double
    inline __m128d operator()(double a){
      return _mm_set_pd(a,a);
    }
    //Integer
    inline __m128i operator()(Integer a){
      return _mm_set1_epi32(a);
    }
  };

  struct Vstore{
    //Float 
    inline void operator()(__m128 a, float* F){
      _mm_store_ps(F,a);
    }
    //Double
    inline void operator()(__m128d a, double* D){
      _mm_store_pd(D,a);
    }
    //Integer
    inline void operator()(__m128i a, Integer* I){
      _mm_store_si128((__m128i *)I,a);
    }

  };


  struct Vstream{
    //Float
    inline void operator()(float * a, __m128 b){
      _mm_stream_ps(a,b);
    }
    //Double
    inline void operator()(double * a, __m128d b){
      _mm_stream_pd(a,b);
    }


  };



  struct Vset{
    // Complex float 
    inline __m128 operator()(Grid::ComplexF *a){
      return _mm_set_ps(a[1].imag(), a[1].real(),a[0].imag(),a[0].real());
    }
    // Complex double 
    inline __m128d operator()(Grid::ComplexD *a){
      return _mm_set_pd(a[0].imag(),a[0].real());
    }
    // Real float 
    inline __m128 operator()(float *a){
      return _mm_set_ps(a[3],a[2],a[1],a[0]);
    }
    // Real double
    inline __m128d operator()(double *a){
      return _mm_set_pd(a[1],a[0]);
    }
    // Integer
    inline __m128i operator()(Integer *a){
      return _mm_set_epi32(a[0],a[1],a[2],a[3]);
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
    inline __m128 operator()(__m128 a, __m128 b){
      return _mm_add_ps(a,b);
    }
    //Complex/Real double
    inline __m128d operator()(__m128d a, __m128d b){
      return _mm_add_pd(a,b);
    }
    //Integer
    inline __m128i operator()(__m128i a, __m128i b){
      return _mm_add_epi32(a,b);
    }
  };

  struct Sub{
    //Complex/Real float
    inline __m128 operator()(__m128 a, __m128 b){
      return _mm_sub_ps(a,b);
    }
    //Complex/Real double
    inline __m128d operator()(__m128d a, __m128d b){
      return _mm_sub_pd(a,b);
    }
    //Integer
    inline __m128i operator()(__m128i a, __m128i b){
      return _mm_sub_epi32(a,b);
    }
  };


  struct MultComplex{
    // Complex float
    inline __m128 operator()(__m128 a, __m128 b){
      __m128 ymm0,ymm1,ymm2;
      ymm0 = _mm_shuffle_ps(a,a,_MM_SHUFFLE(2,2,0,0)); // ymm0 <- ar ar,
      ymm0 = _mm_mul_ps(ymm0,b);                       // ymm0 <- ar bi, ar br
      ymm1 = _mm_shuffle_ps(b,b,_MM_SHUFFLE(2,3,0,1)); // ymm1 <- br,bi
      ymm2 = _mm_shuffle_ps(a,a,_MM_SHUFFLE(3,3,1,1)); // ymm2 <- ai,ai
      ymm1 = _mm_mul_ps(ymm1,ymm2);                    // ymm1 <- br ai, ai bi
      return _mm_addsub_ps(ymm0,ymm1);    
    }
    // Complex double
    inline __m128d operator()(__m128d a, __m128d b){
      __m128d ymm0,ymm1,ymm2;
      ymm0 = _mm_shuffle_pd(a,a,0x0);   // ymm0 <- ar ar,
      ymm0 = _mm_mul_pd(ymm0,b);        // ymm0 <- ar bi, ar br
      ymm1 = _mm_shuffle_pd(b,b,0x1);   // ymm1 <- br,bi   b01
      ymm2 = _mm_shuffle_pd(a,a,0x3);   // ymm2 <- ai,ai   b11
      ymm1 = _mm_mul_pd(ymm1,ymm2);     // ymm1 <- br ai, ai bi
      return _mm_addsub_pd(ymm0,ymm1);  
    }
  };

  struct Mult{
    // Real float
    inline __m128 operator()(__m128 a, __m128 b){
      return _mm_mul_ps(a,b);
    }
    // Real double
    inline __m128d operator()(__m128d a, __m128d b){
      return _mm_mul_pd(a,b);
    }
    // Integer
    inline __m128i operator()(__m128i a, __m128i b){
      return _mm_mul_epi32(a,b);
    }
  };


  struct Conj{
    // Complex single
    inline __m128 operator()(__m128 in){
      return _mm_xor_ps(_mm_addsub_ps(_mm_setzero_ps(),in), _mm_set1_ps(-0.f));
    }
    // Complex double
    inline __m128d operator()(__m128d in){
      return _mm_xor_pd(_mm_addsub_pd(_mm_setzero_pd(),in), _mm_set1_pd(-0.f));//untested
    }
    // do not define for integer input
  };

  struct TimesMinusI{
    //Complex single
    inline __m128 operator()(__m128 in, __m128 ret){
      __m128 tmp =_mm_addsub_ps(_mm_setzero_ps(),in); // r,-i
      return _mm_shuffle_ps(tmp,tmp,_MM_SHUFFLE(2,3,0,1));
    }
    //Complex double
    inline __m128d operator()(__m128d in, __m128d ret){
      __m128d tmp =_mm_addsub_pd(_mm_setzero_pd(),in); // r,-i
      return _mm_shuffle_pd(tmp,tmp,0x1);
    }


  };

  struct TimesI{
    //Complex single
    inline __m128 operator()(__m128 in, __m128 ret){
      __m128 tmp =_mm_shuffle_ps(in,in,_MM_SHUFFLE(2,3,0,1));
      return _mm_addsub_ps(_mm_setzero_ps(),tmp); // r,-i
    }
    //Complex double
    inline __m128d operator()(__m128d in, __m128d ret){
      __m128d tmp = _mm_shuffle_pd(in,in,0x1);
      return _mm_addsub_pd(_mm_setzero_pd(),tmp); // r,-i
    }


  };


  


  //////////////////////////////////////////////
  // Some Template specialization
  
  //Complex float Reduce
  template<>
  inline Grid::ComplexF Reduce<Grid::ComplexF, __m128>::operator()(__m128 in){
    union {
      __m128 v1;  
      float f[4]; 
    } u128;
    u128.v1 = _mm_add_ps(in, _mm_shuffle_ps(in,in, 0b01001110)); // FIXME Prefer to use _MM_SHUFFLE macros
    return Grid::ComplexF(u128.f[0], u128.f[1]);   
  }
  //Real float Reduce
  template<>
  inline Grid::RealF Reduce<Grid::RealF, __m128>::operator()(__m128 in){
    // FIXME Hack
    const Grid::RealF * ptr = (const Grid::RealF *) &in;
    Grid::RealF ret = 0; 
    for(int i=0;i< 4 ;i++){ // 4 number of simd lanes for float
      ret = ret+ptr[i];
    }
    return ret;
  }
  
  
  //Complex double Reduce
  template<>
  inline Grid::ComplexD Reduce<Grid::ComplexD, __m128d>::operator()(__m128d in){
    printf("Reduce : Missing good complex double implementation -> FIX\n");
    return Grid::ComplexD(in[0], in[1]); // inefficient
  }
  
  //Real double Reduce
  template<>
  inline Grid::RealD Reduce<Grid::RealD, __m128d>::operator()(__m128d in){
    // FIXME Hack
    const Grid::RealD * ptr =(const Grid::RealD *)  &in;
    Grid::RealD ret = 0; 
    for(int i=0;i< 2 ;i++){// 2 number of simd lanes for float
      ret = ret+ptr[i];
    }
    return ret;
  }

  //Integer Reduce
  template<>
  inline Integer Reduce<Integer, __m128i>::operator()(__m128i in){
    // FIXME unimplemented
   printf("Reduce : Missing integer implementation -> FIX\n");
    assert(0);
  }
  
  
}

//////////////////////////////////////////////////////////////////////////////////////
// Here assign types 
namespace Grid {
  typedef __m128 SIMD_Ftype;  // Single precision type
  typedef __m128d SIMD_Dtype; // Double precision type
  typedef __m128i SIMD_Itype; // Integer type


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
