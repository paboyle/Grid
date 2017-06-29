/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./lib/simd/Grid_neon.h

    Copyright (C) 2015

Author: Nils Meyer <nils.meyer@ur.de>
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
/*

  ARMv8 NEON intrinsics layer by

  Nils Meyer <nils.meyer@ur.de>,
  University of Regensburg, Germany
  SFB/TRR55

*/
//----------------------------------------------------------------------
//#ifndef ARM_NEON
//#define ARM_NEON

#ifndef GEN_SIMD_WIDTH
#define GEN_SIMD_WIDTH 16u
#endif

#include "Grid_generic_types.h"
#include <arm_neon.h>

namespace Grid {
namespace Optimization {

  template<class vtype>
  union uconv {
    float32x4_t f;
    vtype v;
  };
  union u128f {
    float32x4_t v;
    float f[4];
  };
  union u128d {
    float64x2_t v;
    double f[2];
  };
  // half precision
  union u128h {
    float16x8_t v;
    uint16_t f[8];
  };

  struct Vsplat{
    //Complex float
    inline float32x4_t operator()(float a, float b){
      float tmp[4]={a,b,a,b};
      return vld1q_f32(tmp);
    }
    // Real float
    inline float32x4_t operator()(float a){
      return vdupq_n_f32(a);
    }
    //Complex double
    inline float64x2_t operator()(double a, double b){
      double tmp[2]={a,b};
      return vld1q_f64(tmp);
    }
    //Real double // N:tbc
    inline float64x2_t operator()(double a){
      return vdupq_n_f64(a);
    }
    //Integer // N:tbc
    inline uint32x4_t operator()(Integer a){
      return vdupq_n_u32(a);
    }
  };

  struct Vstore{
    //Float
    inline void operator()(float32x4_t a, float* F){
      vst1q_f32(F, a);
    }
    //Double
    inline void operator()(float64x2_t a, double* D){
      vst1q_f64(D, a);
    }
    //Integer
    inline void operator()(uint32x4_t a, Integer* I){
      vst1q_u32(I, a);
    }

  };

  struct Vstream{ // N:equivalents to _mm_stream_p* in NEON?
    //Float // N:generic
    inline void operator()(float * a, float32x4_t b){
      memcpy(a,&b,4*sizeof(float));
    }
    //Double // N:generic
    inline void operator()(double * a, float64x2_t b){
      memcpy(a,&b,2*sizeof(double));
    }


  };

  // Nils: Vset untested; not used currently in Grid at all;
  // git commit 4a8c4ccfba1d05159348d21a9698028ea847e77b
  struct Vset{
    // Complex float // N:ok
    inline float32x4_t operator()(Grid::ComplexF *a){
      float tmp[4]={a[1].imag(),a[1].real(),a[0].imag(),a[0].real()};
      return vld1q_f32(tmp);
    }
    // Complex double // N:ok
    inline float64x2_t operator()(Grid::ComplexD *a){
      double tmp[2]={a[0].imag(),a[0].real()};
      return vld1q_f64(tmp);
    }
    // Real float // N:ok
    inline float32x4_t operator()(float *a){
      float tmp[4]={a[3],a[2],a[1],a[0]};
      return vld1q_f32(tmp);
    }
    // Real double // N:ok
    inline float64x2_t operator()(double *a){
      double tmp[2]={a[1],a[0]};
      return vld1q_f64(tmp);
    }
    // Integer // N:ok
    inline uint32x4_t operator()(Integer *a){
      return vld1q_dup_u32(a);
    }
  };

  // N:leaving as is
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
    inline float32x4_t operator()(float32x4_t a, float32x4_t b){
      return vaddq_f32(a,b);
    }
    //Complex/Real double
    inline float64x2_t operator()(float64x2_t a, float64x2_t b){
      return vaddq_f64(a,b);
    }
    //Integer
    inline uint32x4_t operator()(uint32x4_t a, uint32x4_t b){
      return vaddq_u32(a,b);
    }
  };

  struct Sub{
    //Complex/Real float
    inline float32x4_t operator()(float32x4_t a, float32x4_t b){
      return vsubq_f32(a,b);
    }
    //Complex/Real double
    inline float64x2_t operator()(float64x2_t a, float64x2_t b){
      return vsubq_f64(a,b);
    }
    //Integer
    inline uint32x4_t operator()(uint32x4_t a, uint32x4_t b){
      return vsubq_u32(a,b);
    }
  };

  struct MultRealPart{
    inline float32x4_t operator()(float32x4_t a, float32x4_t b){
      float32x4_t re = vtrn1q_f32(a, a);
      return vmulq_f32(re, b);
    }
    inline float64x2_t operator()(float64x2_t a, float64x2_t b){
      float64x2_t re = vzip1q_f64(a, a);
      return vmulq_f64(re, b);
    }
  };

  struct MaddRealPart{
    inline float32x4_t operator()(float32x4_t a, float32x4_t b, float32x4_t c){
      float32x4_t re = vtrn1q_f32(a, a);
      return vfmaq_f32(c, re, b);
    }
    inline float64x2_t operator()(float64x2_t a, float64x2_t b, float64x2_t c){
      float64x2_t re = vzip1q_f64(a, a);
      return vfmaq_f64(c, re, b);
    }
  };

  struct Div{
    // Real float
    inline float32x4_t operator()(float32x4_t a, float32x4_t b){
      return vdivq_f32(a, b);
    }
    // Real double
    inline float64x2_t operator()(float64x2_t a, float64x2_t b){
      return vdivq_f64(a, b);
    }
  };

  struct MultComplex{
    // Complex float
    inline float32x4_t operator()(float32x4_t a, float32x4_t b){

      float32x4_t r0, r1, r2, r3, r4;

      // a = ar ai Ar Ai
      // b = br bi Br Bi
      // collect real/imag part, negate bi and Bi
      r0 = vtrn1q_f32(b, b);       //  br  br  Br  Br
      r1 = vnegq_f32(b);           // -br -bi -Br -Bi
      r2 = vtrn2q_f32(b, r1);      //  bi -bi  Bi -Bi

      // the fun part
      r3 = vmulq_f32(r2, a);       //  bi*ar -bi*ai ...
      r4 = vrev64q_f32(r3);        // -bi*ai  bi*ar ...

      // fma(a,b,c) = a+b*c
      return vfmaq_f32(r4, r0, a); //  ar*br-ai*bi ai*br+ar*bi ...

      // no fma, use mul and add
      //float32x4_t r5;
      //r5 = vmulq_f32(r0, a);
      //return vaddq_f32(r4, r5);
    }
    // Complex double
    inline float64x2_t operator()(float64x2_t a, float64x2_t b){

      float64x2_t r0, r1, r2, r3, r4;

      // b = br bi
      // collect real/imag part, negate bi
      r0 = vtrn1q_f64(b, b);       //  br  br
      r1 = vnegq_f64(b);           // -br -bi
      r2 = vtrn2q_f64(b, r1);      //  bi -bi

      // the fun part
      r3 = vmulq_f64(r2, a);       //  bi*ar -bi*ai
      r4 = vextq_f64(r3,r3,1);     // -bi*ai  bi*ar

      // fma(a,b,c) = a+b*c
      return vfmaq_f64(r4, r0, a); //  ar*br-ai*bi ai*br+ar*bi

      // no fma, use mul and add
      //float64x2_t r5;
      //r5 = vmulq_f64(r0, a);
      //return vaddq_f64(r4, r5);
    }
  };

  struct Mult{
    // Real float
    inline float32x4_t mac(float32x4_t a, float32x4_t b, float32x4_t c){
      //return vaddq_f32(vmulq_f32(b,c),a);
      return vfmaq_f32(a, b, c);
    }
    inline float64x2_t mac(float64x2_t a, float64x2_t b, float64x2_t c){
      //return vaddq_f64(vmulq_f64(b,c),a);
      return vfmaq_f64(a, b, c);
    }
    inline float32x4_t operator()(float32x4_t a, float32x4_t b){
      return vmulq_f32(a,b);
    }
    // Real double
    inline float64x2_t operator()(float64x2_t a, float64x2_t b){
      return vmulq_f64(a,b);
    }
    // Integer
    inline uint32x4_t operator()(uint32x4_t a, uint32x4_t b){
      return vmulq_u32(a,b);
    }
  };

  struct Conj{
    // Complex single
    inline float32x4_t operator()(float32x4_t in){
      // ar ai br bi -> ar -ai br -bi
      float32x4_t r0, r1;
      r0 = vnegq_f32(in);        // -ar -ai -br -bi
      r1 = vrev64q_f32(r0);      // -ai -ar -bi -br
      return vtrn1q_f32(in, r1); //  ar -ai  br -bi
    }
    // Complex double
    inline float64x2_t operator()(float64x2_t in){

      float64x2_t r0, r1;
      r0 = vextq_f64(in, in, 1);    //  ai  ar
      r1 = vnegq_f64(r0);           // -ai -ar
      return vextq_f64(r0, r1, 1);  //  ar -ai
    }
    // do not define for integer input
  };

  struct TimesMinusI{
    //Complex single
    inline float32x4_t operator()(float32x4_t in, float32x4_t ret){
      // ar ai br bi -> ai -ar ai -br
      float32x4_t r0, r1;
      r0 = vnegq_f32(in);        // -ar -ai -br -bi
      r1 = vrev64q_f32(in);      //  ai  ar  bi  br
      return vtrn1q_f32(r1, r0); //  ar -ai  br -bi
    }
    //Complex double
    inline float64x2_t operator()(float64x2_t in, float64x2_t ret){
      // a ib -> b -ia
      float64x2_t tmp;
      tmp = vnegq_f64(in);
      return vextq_f64(in, tmp, 1);
    }
  };

  struct TimesI{
    //Complex single
    inline float32x4_t operator()(float32x4_t in, float32x4_t ret){
      // ar ai br bi -> -ai ar -bi br
      float32x4_t r0, r1;
      r0 = vnegq_f32(in);        // -ar -ai -br -bi
      r1 = vrev64q_f32(r0);      // -ai -ar -bi -br
      return vtrn1q_f32(r1, in); // -ai  ar -bi  br
    }
    //Complex double
    inline float64x2_t operator()(float64x2_t in, float64x2_t ret){
      // a ib -> -b ia
      float64x2_t tmp;
      tmp = vnegq_f64(in);
      return vextq_f64(tmp, in, 1);
    }
  };

  struct Permute{

    static inline float32x4_t Permute0(float32x4_t in){ // N:ok
      // AB CD -> CD AB
      return vextq_f32(in, in, 2);
    };
    static inline float32x4_t Permute1(float32x4_t in){ // N:ok
      // AB CD -> BA DC
      return vrev64q_f32(in);
    };
    static inline float32x4_t Permute2(float32x4_t in){ // N:not used by Boyle
      return in;
    };
    static inline float32x4_t Permute3(float32x4_t in){ // N:not used by Boyle
      return in;
    };

    static inline float64x2_t Permute0(float64x2_t in){ // N:ok
      // AB -> BA
      return vextq_f64(in, in, 1);
    };
    static inline float64x2_t Permute1(float64x2_t in){ // N:not used by Boyle
      return in;
    };
    static inline float64x2_t Permute2(float64x2_t in){ // N:not used by Boyle
      return in;
    };
    static inline float64x2_t Permute3(float64x2_t in){ // N:not used by Boyle
      return in;
    };

  };

  struct Rotate{

    static inline float32x4_t rotate(float32x4_t in,int n){ // N:ok
      switch(n){
      case 0: // AB CD -> AB CD
        return tRotate<0>(in);
        break;
      case 1: // AB CD -> BC DA
        return tRotate<1>(in);
        break;
      case 2: // AB CD -> CD AB
        return tRotate<2>(in);
        break;
      case 3: // AB CD -> DA BC
        return tRotate<3>(in);
        break;
      default: assert(0);
      }
    }
    static inline float64x2_t rotate(float64x2_t in,int n){ // N:ok
      switch(n){
      case 0: // AB -> AB
        return tRotate<0>(in);
        break;
      case 1: // AB -> BA
        return tRotate<1>(in);
        break;
      default: assert(0);
      }
    }

// working, but no restriction on n
//    template<int n> static inline float32x4_t tRotate(float32x4_t in){ return vextq_f32(in,in,n); };
//    template<int n> static inline float64x2_t tRotate(float64x2_t in){ return vextq_f64(in,in,n); };

// restriction on n
    template<int n> static inline float32x4_t tRotate(float32x4_t in){ return vextq_f32(in,in,n%4); };
    template<int n> static inline float64x2_t tRotate(float64x2_t in){ return vextq_f64(in,in,n%2); };

  };

  struct PrecisionChange {

    static inline float16x8_t StoH (const float32x4_t &a,const float32x4_t &b) {
      float16x4_t h = vcvt_f16_f32(a);
      return vcvt_high_f16_f32(h, b);
    }
    static inline void  HtoS (float16x8_t h,float32x4_t &sa,float32x4_t &sb) {
      sb = vcvt_high_f32_f16(h);
      // there is no direct conversion from lower float32x4_t to float64x2_t
      // vextq_f16 not supported by clang 3.8 / 4.0 / arm clang
      //float16x8_t h1 = vextq_f16(h, h, 4); // correct, but not supported by clang
      // workaround for clang
      uint32x4_t h1u = reinterpret_cast<uint32x4_t>(h);
      float16x8_t h1 = reinterpret_cast<float16x8_t>(vextq_u32(h1u, h1u, 2));
      sa = vcvt_high_f32_f16(h1);
    }
    static inline float32x4_t DtoS (float64x2_t a,float64x2_t b) {
      float32x2_t s = vcvt_f32_f64(a);
      return vcvt_high_f32_f64(s, b);

    }
    static inline void StoD (float32x4_t s,float64x2_t &a,float64x2_t &b) {
      b = vcvt_high_f64_f32(s);
      // there is no direct conversion from lower float32x4_t to float64x2_t
      float32x4_t s1 = vextq_f32(s, s, 2);
      a = vcvt_high_f64_f32(s1);

    }
    static inline float16x8_t DtoH (float64x2_t a,float64x2_t b,float64x2_t c,float64x2_t d) {
      float32x4_t s1 = DtoS(a, b);
      float32x4_t s2 = DtoS(c, d);
      return StoH(s1, s2);
    }
    static inline void HtoD (float16x8_t h,float64x2_t &a,float64x2_t &b,float64x2_t &c,float64x2_t &d) {
      float32x4_t s1, s2;
      HtoS(h, s1, s2);
      StoD(s1, a, b);
      StoD(s2, c, d);
    }
  };

  //////////////////////////////////////////////
  // Exchange support

  struct Exchange{
    static inline void Exchange0(float32x4_t &out1,float32x4_t &out2,float32x4_t in1,float32x4_t in2){
      // in1: ABCD -> out1: ABEF
      // in2: EFGH -> out2: CDGH

      // z: CDAB
      float32x4_t z = vextq_f32(in1, in1, 2);
      // out1: ABEF
      out1 = vextq_f32(z, in2, 2);

      // z: GHEF
      z = vextq_f32(in2, in2, 2);
      // out2: CDGH
      out2 = vextq_f32(in1, z, 2);
    };

    static inline void Exchange1(float32x4_t &out1,float32x4_t &out2,float32x4_t in1,float32x4_t in2){
      // in1: ABCD -> out1: AECG
      // in2: EFGH -> out2: BFDH
      out1 = vtrn1q_f32(in1, in2);
      out2 = vtrn2q_f32(in1, in2);
    };
    static inline void Exchange2(float32x4_t &out1,float32x4_t &out2,float32x4_t in1,float32x4_t in2){
      assert(0);
      return;
    };
    static inline void Exchange3(float32x4_t &out1,float32x4_t &out2,float32x4_t in1,float32x4_t in2){
      assert(0);
      return;
    };
    // double precision
    static inline void Exchange0(float64x2_t &out1,float64x2_t &out2,float64x2_t in1,float64x2_t in2){
      // in1: AB -> out1: AC
      // in2: CD -> out2: BD
      out1 = vzip1q_f64(in1, in2);
      out2 = vzip2q_f64(in1, in2);
    };
    static inline void Exchange1(float64x2_t &out1,float64x2_t &out2,float64x2_t in1,float64x2_t in2){
      assert(0);
      return;
    };
    static inline void Exchange2(float64x2_t &out1,float64x2_t &out2,float64x2_t in1,float64x2_t in2){
      assert(0);
      return;
    };
    static inline void Exchange3(float64x2_t &out1,float64x2_t &out2,float64x2_t in1,float64x2_t in2){
      assert(0);
      return;
    };
  };

  //////////////////////////////////////////////
  // Some Template specialization


  //Complex float Reduce
  template<>
  inline Grid::ComplexF Reduce<Grid::ComplexF, float32x4_t>::operator()(float32x4_t in){
    float32x4_t v1; // two complex
    v1 = Optimization::Permute::Permute0(in);
    v1 = vaddq_f32(v1,in);
    u128f conv;    conv.v=v1;
    return Grid::ComplexF(conv.f[0],conv.f[1]);
  }
  //Real float Reduce
  template<>
  inline Grid::RealF Reduce<Grid::RealF, float32x4_t>::operator()(float32x4_t in){
    return vaddvq_f32(in);
  }


  //Complex double Reduce
  template<> // N:by Boyle
  inline Grid::ComplexD Reduce<Grid::ComplexD, float64x2_t>::operator()(float64x2_t in){
    u128d conv; conv.v = in;
    return Grid::ComplexD(conv.f[0],conv.f[1]);
  }

  //Real double Reduce
  template<>
  inline Grid::RealD Reduce<Grid::RealD, float64x2_t>::operator()(float64x2_t in){
    return vaddvq_f64(in);
  }

  //Integer Reduce
  template<>
  inline Integer Reduce<Integer, uint32x4_t>::operator()(uint32x4_t in){
    // FIXME unimplemented
    printf("Reduce : Missing integer implementation -> FIX\n");
    assert(0);
  }
}

//////////////////////////////////////////////////////////////////////////////////////
// Here assign types

// typedef Optimization::vech SIMD_Htype; // Reduced precision type
  typedef float16x8_t  SIMD_Htype; // Half precision type
  typedef float32x4_t  SIMD_Ftype; // Single precision type
  typedef float64x2_t  SIMD_Dtype; // Double precision type
  typedef uint32x4_t   SIMD_Itype; // Integer type

  inline void v_prefetch0(int size, const char *ptr){};  // prefetch utilities
  inline void prefetch_HINT_T0(const char *ptr){};


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

//#endif // ARM_NEON
