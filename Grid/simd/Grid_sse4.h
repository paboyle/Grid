    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/simd/Grid_sse4.h

    Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: neo <cossu@post.kek.jp>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

    See the full license in the file "LICENSE" in the top level distribution directory
    *************************************************************************************/
    /*  END LEGAL */
//----------------------------------------------------------------------
/*! @file Grid_sse4.h
  @brief Optimization libraries for SSE4 instructions set

  Using intrinsics
*/
// Time-stamp: <2015-06-16 23:27:54 neo>
//----------------------------------------------------------------------

#include <pmmintrin.h>

namespace Grid {
namespace Optimization {

  template<class vtype>
  union uconv {
    __m128 f;
    vtype v;
  };

  union u128f {
    __m128 v;
    float f[4];
  };
  union u128d {
    __m128d v;
    double f[2];
  };
  
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
      return _mm_set_epi32(a[3],a[2],a[1],a[0]);
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

  struct MultRealPart{
    inline __m128 operator()(__m128 a, __m128 b){
      __m128 ymm0;
      ymm0  = _mm_shuffle_ps(a,a,_MM_SELECT_FOUR_FOUR(2,2,0,0)); // ymm0 <- ar ar,
      return  _mm_mul_ps(ymm0,b);                       // ymm0 <- ar bi, ar br
    }
    inline __m128d operator()(__m128d a, __m128d b){
      __m128d ymm0;
      ymm0 = _mm_shuffle_pd(a,a,0x0); // ymm0 <- ar ar, ar,ar b'00,00
      return _mm_mul_pd(ymm0,b);      // ymm0 <- ar bi, ar br
    }
  };
  struct MaddRealPart{
    inline __m128 operator()(__m128 a, __m128 b, __m128 c){
      __m128 ymm0 =  _mm_shuffle_ps(a,a,_MM_SELECT_FOUR_FOUR(2,2,0,0)); // ymm0 <- ar ar,
      return _mm_add_ps(_mm_mul_ps( ymm0, b),c);                         
    }
    inline __m128d operator()(__m128d a, __m128d b, __m128d c){
      __m128d ymm0 = _mm_shuffle_pd( a, a, 0x0 );
      return _mm_add_pd(_mm_mul_pd( ymm0, b),c);                         
    }
  };

  struct MultComplex{
    // Complex float
    inline __m128 operator()(__m128 a, __m128 b){
      __m128 ymm0,ymm1,ymm2;
      ymm0 = _mm_shuffle_ps(a,a,_MM_SELECT_FOUR_FOUR(2,2,0,0)); // ymm0 <- ar ar,
      ymm0 = _mm_mul_ps(ymm0,b);                       // ymm0 <- ar bi, ar br
      ymm1 = _mm_shuffle_ps(b,b,_MM_SELECT_FOUR_FOUR(2,3,0,1)); // ymm1 <- br,bi
      ymm2 = _mm_shuffle_ps(a,a,_MM_SELECT_FOUR_FOUR(3,3,1,1)); // ymm2 <- ai,ai
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

    inline void mac(__m128 &a, __m128 b, __m128 c){
      a= _mm_add_ps(_mm_mul_ps(b,c),a);
    }

    inline void mac(__m128d &a, __m128d b, __m128d c){
      a= _mm_add_pd(_mm_mul_pd(b,c),a);
    }

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
      return _mm_mullo_epi32(a,b);
    }
  };

  struct Div{
    // Real float
    inline __m128 operator()(__m128 a, __m128 b){
      return _mm_div_ps(a,b);
    }
    // Real double
    inline __m128d operator()(__m128d a, __m128d b){
      return _mm_div_pd(a,b);
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
      return _mm_shuffle_ps(tmp,tmp,_MM_SELECT_FOUR_FOUR(2,3,0,1));
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
      __m128 tmp =_mm_shuffle_ps(in,in,_MM_SELECT_FOUR_FOUR(2,3,0,1));
      return _mm_addsub_ps(_mm_setzero_ps(),tmp); // r,-i
    }
    //Complex double
    inline __m128d operator()(__m128d in, __m128d ret){
      __m128d tmp = _mm_shuffle_pd(in,in,0x1);
      return _mm_addsub_pd(_mm_setzero_pd(),tmp); // r,-i
    }
  };

  struct Permute{

    static inline __m128 Permute0(__m128 in){
      return _mm_shuffle_ps(in,in,_MM_SELECT_FOUR_FOUR(1,0,3,2)); //AB CD -> CD AB
    };
    static inline __m128 Permute1(__m128 in){
      return _mm_shuffle_ps(in,in,_MM_SELECT_FOUR_FOUR(2,3,0,1)); //AB CD -> BA DC
    };
    static inline __m128 Permute2(__m128 in){
      return in;
    };
    static inline __m128 Permute3(__m128 in){
      return in;
    };

    static inline __m128d Permute0(__m128d in){ //AB -> BA
      return _mm_shuffle_pd(in,in,0x1);
    };
    static inline __m128d Permute1(__m128d in){
      return in;
    };
    static inline __m128d Permute2(__m128d in){
      return in;
    };
    static inline __m128d Permute3(__m128d in){
      return in;
    };
  };

  
#define _my_alignr_epi32(a,b,n) _mm_alignr_epi8(a,b,(n*4)%16)
#define _my_alignr_epi64(a,b,n) _mm_alignr_epi8(a,b,(n*8)%16)

#ifdef SFW_FP16

  struct Grid_half {
    Grid_half(){}
    Grid_half(uint16_t raw) : x(raw) {}
    uint16_t x;
  };
  union FP32 {
    unsigned int u;
    float f;
  };

  // PAB - Lifted and adapted from Eigen, which is GPL V2
  inline float sfw_half_to_float(Grid_half h) {
    const FP32 magic = { 113 << 23 };
    const unsigned int shifted_exp = 0x7c00 << 13; // exponent mask after shift
    FP32 o;
    o.u = (h.x & 0x7fff) << 13;             // exponent/mantissa bits
    unsigned int exp = shifted_exp & o.u;   // just the exponent
    o.u += (127 - 15) << 23;                // exponent adjust
    // handle exponent special cases
    if (exp == shifted_exp) {     // Inf/NaN?
      o.u += (128 - 16) << 23;    // extra exp adjust
    } else if (exp == 0) {        // Zero/Denormal?
      o.u += 1 << 23;             // extra exp adjust
      o.f -= magic.f;             // renormalize
    }
    o.u |= (h.x & 0x8000) << 16;    // sign bit
    return o.f;
  }
  inline Grid_half sfw_float_to_half(float ff) {
    FP32 f; f.f = ff;
    const FP32 f32infty = { 255 << 23 };
    const FP32 f16max = { (127 + 16) << 23 };
    const FP32 denorm_magic = { ((127 - 15) + (23 - 10) + 1) << 23 };
    unsigned int sign_mask = 0x80000000u;
    Grid_half o;
    
    o.x = static_cast<unsigned short>(0x0u);
    unsigned int sign = f.u & sign_mask;
    f.u ^= sign;
    // NOTE all the integer compares in this function can be safely
    // compiled into signed compares since all operands are below
    // 0x80000000. Important if you want fast straight SSE2 code
    // (since there's no unsigned PCMPGTD).
    if (f.u >= f16max.u) {  // result is Inf or NaN (all exponent bits set)
      o.x = (f.u > f32infty.u) ? 0x7e00 : 0x7c00; // NaN->qNaN and Inf->Inf
    } else {  // (De)normalized number or zero
      if (f.u < (113 << 23)) {  // resulting FP16 is subnormal or zero
	// use a magic value to align our 10 mantissa bits at the bottom of
	// the float. as long as FP addition is round-to-nearest-even this
	// just works.
	f.f += denorm_magic.f;
	// and one integer subtract of the bias later, we have our final float!
	o.x = static_cast<unsigned short>(f.u - denorm_magic.u);
      } else {
	unsigned int mant_odd = (f.u >> 13) & 1; // resulting mantissa is odd
	
	// update exponent, rounding bias part 1
	f.u += ((unsigned int)(15 - 127) << 23) + 0xfff;
	// rounding bias part 2
	f.u += mant_odd;
	// take the bits!
	o.x = static_cast<unsigned short>(f.u >> 13);
      }
    } 
    o.x |= static_cast<unsigned short>(sign >> 16);
    return o;
  }
  static inline __m128i Grid_mm_cvtps_ph(__m128 f,int discard) {
    __m128i ret=(__m128i)_mm_setzero_ps();
    float *fp = (float *)&f;
    Grid_half *hp = (Grid_half *)&ret;
    hp[0] = sfw_float_to_half(fp[0]);
    hp[1] = sfw_float_to_half(fp[1]);
    hp[2] = sfw_float_to_half(fp[2]);
    hp[3] = sfw_float_to_half(fp[3]);
    return ret;
  }
  static inline __m128 Grid_mm_cvtph_ps(__m128i h,int discard) {
    __m128 ret=_mm_setzero_ps();
    float *fp = (float *)&ret;
    Grid_half  *hp = (Grid_half *)&h;
    fp[0] = sfw_half_to_float(hp[0]);
    fp[1] = sfw_half_to_float(hp[1]);
    fp[2] = sfw_half_to_float(hp[2]);
    fp[3] = sfw_half_to_float(hp[3]);
    return ret;
  }
#else 
#define Grid_mm_cvtps_ph _mm_cvtps_ph
#define Grid_mm_cvtph_ps _mm_cvtph_ps
#endif
  struct PrecisionChange {
    static inline __m128i StoH (__m128 a,__m128 b) {
      __m128i ha = Grid_mm_cvtps_ph(a,0);
      __m128i hb = Grid_mm_cvtps_ph(b,0);
      __m128i h =(__m128i) _mm_shuffle_ps((__m128)ha,(__m128)hb,_MM_SELECT_FOUR_FOUR(1,0,1,0));
      return h;
    }
    static inline void  HtoS (__m128i h,__m128 &sa,__m128 &sb) {
      sa = Grid_mm_cvtph_ps(h,0); 
      h =  (__m128i)_my_alignr_epi32((__m128i)h,(__m128i)h,2);
      sb = Grid_mm_cvtph_ps(h,0);
    }
    static inline __m128 DtoS (__m128d a,__m128d b) {
      __m128 sa = _mm_cvtpd_ps(a);
      __m128 sb = _mm_cvtpd_ps(b);
      __m128 s = _mm_shuffle_ps(sa,sb,_MM_SELECT_FOUR_FOUR(1,0,1,0));
      return s;
    }
    static inline void StoD (__m128 s,__m128d &a,__m128d &b) {
      a = _mm_cvtps_pd(s);
      s = (__m128)_my_alignr_epi32((__m128i)s,(__m128i)s,2);
      b = _mm_cvtps_pd(s);
    }
    static inline __m128i DtoH (__m128d a,__m128d b,__m128d c,__m128d d) {
      __m128 sa,sb;
      sa = DtoS(a,b);
      sb = DtoS(c,d);
      return StoH(sa,sb);
    }
    static inline void HtoD (__m128i h,__m128d &a,__m128d &b,__m128d &c,__m128d &d) {
      __m128 sa,sb;
      HtoS(h,sa,sb);
      StoD(sa,a,b);
      StoD(sb,c,d);
    }
  };

  struct Exchange{
    // 3210 ordering
    static inline void Exchange0(__m128 &out1,__m128 &out2,__m128 in1,__m128 in2){
      out1= _mm_shuffle_ps(in1,in2,_MM_SELECT_FOUR_FOUR(1,0,1,0));
      out2= _mm_shuffle_ps(in1,in2,_MM_SELECT_FOUR_FOUR(3,2,3,2));
    };
    static inline void Exchange1(__m128 &out1,__m128 &out2,__m128 in1,__m128 in2){
      out1= _mm_shuffle_ps(in1,in2,_MM_SELECT_FOUR_FOUR(2,0,2,0)); /*ACEG*/
      out2= _mm_shuffle_ps(in1,in2,_MM_SELECT_FOUR_FOUR(3,1,3,1)); /*BDFH*/
      out1= _mm_shuffle_ps(out1,out1,_MM_SELECT_FOUR_FOUR(3,1,2,0)); /*AECG*/
      out2= _mm_shuffle_ps(out2,out2,_MM_SELECT_FOUR_FOUR(3,1,2,0)); /*AECG*/
    };
    static inline void Exchange2(__m128 &out1,__m128 &out2,__m128 in1,__m128 in2){
      assert(0);
      return;
    };
    static inline void Exchange3(__m128 &out1,__m128 &out2,__m128 in1,__m128 in2){
      assert(0);
      return;
    };

    static inline void Exchange0(__m128d &out1,__m128d &out2,__m128d in1,__m128d in2){
      out1= _mm_shuffle_pd(in1,in2,0x0);
      out2= _mm_shuffle_pd(in1,in2,0x3);
    };
    static inline void Exchange1(__m128d &out1,__m128d &out2,__m128d in1,__m128d in2){
      assert(0);
      return;
    };
    static inline void Exchange2(__m128d &out1,__m128d &out2,__m128d in1,__m128d in2){
      assert(0);
      return;
    };
    static inline void Exchange3(__m128d &out1,__m128d &out2,__m128d in1,__m128d in2){
      assert(0);
      return;
    };
  };

  struct Rotate{

    static inline __m128 rotate(__m128 in,int n){ 
      switch(n){
      case 0: return tRotate<0>(in);break;
      case 1: return tRotate<1>(in);break;
      case 2: return tRotate<2>(in);break;
      case 3: return tRotate<3>(in);break;
      default: assert(0);
      }
    }
    static inline __m128d rotate(__m128d in,int n){ 
      switch(n){
      case 0: return tRotate<0>(in);break;
      case 1: return tRotate<1>(in);break;
      default: assert(0);
      }
    }

    template<int n> static inline __m128  tRotate(__m128  in){ return (__m128)_my_alignr_epi32((__m128i)in,(__m128i)in,n); };
    template<int n> static inline __m128d tRotate(__m128d in){ return (__m128d)_my_alignr_epi64((__m128i)in,(__m128i)in,n); };

  };
  //////////////////////////////////////////////
  // Some Template specialization


  //Complex float Reduce
  template<>
  inline Grid::ComplexF Reduce<Grid::ComplexF, __m128>::operator()(__m128 in){
    __m128 v1; // two complex
    v1= Optimization::Permute::Permute0(in); 
    v1= _mm_add_ps(v1,in);
    u128f conv;    conv.v=v1;
    return Grid::ComplexF(conv.f[0],conv.f[1]);
  }
  //Real float Reduce
  template<>
  inline Grid::RealF Reduce<Grid::RealF, __m128>::operator()(__m128 in){
    __m128 v1,v2; // quad single
    v1= Optimization::Permute::Permute0(in); 
    v1= _mm_add_ps(v1,in);
    v2= Optimization::Permute::Permute1(v1); 
    v1 = _mm_add_ps(v1,v2);
    u128f conv; conv.v=v1;
    return conv.f[0];
  }
  
  
  //Complex double Reduce
  template<>
  inline Grid::ComplexD Reduce<Grid::ComplexD, __m128d>::operator()(__m128d in){
    u128d conv; conv.v = in;
    return Grid::ComplexD(conv.f[0],conv.f[1]);
  }
  
  //Real double Reduce
  template<>
  inline Grid::RealD Reduce<Grid::RealD, __m128d>::operator()(__m128d in){
    __m128d v1;
    v1 = Optimization::Permute::Permute0(in); 
    v1 = _mm_add_pd(v1,in);
    u128d conv; conv.v = v1;
    return conv.f[0];
  }

  //Integer Reduce
  template<>
  inline Integer Reduce<Integer, __m128i>::operator()(__m128i in){
    __m128i v1 = _mm_hadd_epi32(in, in);
    __m128i v2 = _mm_hadd_epi32(v1, v1);
    return _mm_cvtsi128_si32(v2);
  }
}



//////////////////////////////////////////////////////////////////////////////////////
// Here assign types 

  typedef __m128i SIMD_Htype;  // Single precision type
  typedef __m128  SIMD_Ftype;  // Single precision type
  typedef __m128d SIMD_Dtype; // Double precision type
  typedef __m128i SIMD_Itype; // Integer type

  // prefetch utilities
  inline void v_prefetch0(int size, const char *ptr){};
  inline void prefetch_HINT_T0(const char *ptr){
    _mm_prefetch(ptr,_MM_HINT_T0);
  }

  // Function name aliases
  typedef Optimization::Vsplat   VsplatSIMD;
  typedef Optimization::Vstore   VstoreSIMD;
  typedef Optimization::Vset     VsetSIMD;
  typedef Optimization::Vstream  VstreamSIMD;
  template <typename S, typename T> using ReduceSIMD = Optimization::Reduce<S,T>;

 


  // Arithmetic operations
  typedef Optimization::Sum         SumSIMD;
  typedef Optimization::Sub         SubSIMD;
  typedef Optimization::Div         DivSIMD;
  typedef Optimization::Mult        MultSIMD;
  typedef Optimization::MultComplex MultComplexSIMD;
  typedef Optimization::MultRealPart MultRealPartSIMD;
  typedef Optimization::MaddRealPart MaddRealPartSIMD;
  typedef Optimization::Conj        ConjSIMD;
  typedef Optimization::TimesMinusI TimesMinusISIMD;
  typedef Optimization::TimesI      TimesISIMD;

}
