//----------------------------------------------------------------------
/*! @file Grid_sse4.h
  @brief Optimization libraries
*/
// Time-stamp: <2015-05-19 17:06:51 neo>
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

  struct Reduce{
    //Complex float
    inline Grid::ComplexF operator()(__m128 in){
      union {
	__m128 v1;  
	float f[4]; 
      } u128;
      u128.v1 = _mm_add_ps(in, _mm_shuffle_ps(in,in, 0b01001110)); // FIXME Prefer to use _MM_SHUFFLE macros
      return Grid::ComplexF(u128.f[0], u128.f[1]);   
    }
    //Complex double
    inline Grid::ComplexD operator()(__m128d in){
      printf("Missing complex double implementation -> FIX\n");
      return Grid::ComplexD(0,0); // FIXME wrong
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
      ymm0 = _mm_mul_ps(ymm0,b);        // ymm0 <- ar bi, ar br
      ymm1 = _mm_shuffle_ps(b,b,_MM_SHUFFLE(2,3,0,1)); // ymm1 <- br,bi
      ymm2 = _mm_shuffle_ps(a,a,_MM_SHUFFLE(3,3,1,1)); // ymm2 <- ai,ai
      ymm1 = _mm_mul_ps(ymm1,ymm2);       // ymm1 <- br ai, ai bi
      return _mm_addsub_ps(ymm0,ymm1);    
    }
    // Complex double
    inline __m128d operator()(__m128d a, __m128d b){
      __m128d ymm0,ymm1,ymm2;
      ymm0 = _mm_shuffle_pd(a,a,0x0); // ymm0 <- ar ar,
      ymm0 = _mm_mul_pd(ymm0,b);        // ymm0 <- ar bi, ar br
      ymm1 = _mm_shuffle_pd(b,b,0x1); // ymm1 <- br,bi   b01
      ymm2 = _mm_shuffle_pd(a,a,0x3); // ymm2 <- ai,ai   b11
      ymm1 = _mm_mul_pd(ymm1,ymm2);       // ymm1 <- br ai, ai bi
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




}

// Here assign types 
namespace Grid {
  typedef __m128 SIMD_Ftype;  // Single precision type
  typedef __m128d SIMD_Dtype; // Double precision type
  typedef __m128i SIMD_Itype; // Integer type


  // Function names
  typedef Optimization::Vsplat VsplatSIMD;
  typedef Optimization::Vstore VstoreSIMD;

  // Arithmetic operations
  typedef Optimization::Sum         SumSIMD;
  typedef Optimization::Sub         SubSIMD;
  typedef Optimization::Mult        MultSIMD;
  typedef Optimization::MultComplex MultComplexSIMD;
  typedef Optimization::Vset        VsetSIMD;

}
