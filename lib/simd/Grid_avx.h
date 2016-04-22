    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/simd/Grid_avx.h

    Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
Author: Guido Cossu <cossu@iroiro-pc.kek.jp>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: neo <cossu@post.kek.jp>
Author: paboyle <paboyle@ph.ed.ac.uk>

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
/*! @file Grid_avx.h
  @brief Optimization libraries for AVX1/2 instructions set

  Using intrinsics
*/
// Time-stamp: <2015-06-16 23:30:41 neo>
//----------------------------------------------------------------------

#include <immintrin.h>
#ifdef AVXFMA4
#include <x86intrin.h>
#endif
// _mm256_set_m128i(hi,lo); // not defined in all versions of immintrin.h
#ifndef _mm256_set_m128i
#define _mm256_set_m128i(hi,lo) _mm256_insertf128_si256(_mm256_castsi128_si256(lo),(hi),1)
#endif

namespace Grid {
namespace Optimization {

  template<class vtype>
  union uconv {
    __m256 f;
    vtype v;
  };

  union u256f {
    __m256 v;
    float f[8];
  };

  union u256d {
    __m256d v;
    double f[4];
  };

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
#if defined (AVX1) || defined (AVXFMA4)
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
#if defined (AVX1) || defined (AVXFMA4)
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
#if defined (AVX1) 
      __m256 ymm0,ymm1,ymm2;
      ymm0 = _mm256_shuffle_ps(a,a,_MM_SELECT_FOUR_FOUR(2,2,0,0)); // ymm0 <- ar ar,
      ymm0 = _mm256_mul_ps(ymm0,b);                       // ymm0 <- ar bi, ar br
      // FIXME AVX2 could MAC
      ymm1 = _mm256_shuffle_ps(b,b,_MM_SELECT_FOUR_FOUR(2,3,0,1)); // ymm1 <- br,bi
      ymm2 = _mm256_shuffle_ps(a,a,_MM_SELECT_FOUR_FOUR(3,3,1,1)); // ymm2 <- ai,ai
      ymm1 = _mm256_mul_ps(ymm1,ymm2);                    // ymm1 <- br ai, ai bi
      return _mm256_addsub_ps(ymm0,ymm1);  
#endif
#if defined (AVXFMA4)
      __m256 a_real = _mm256_shuffle_ps(a,a,_MM_SELECT_FOUR_FOUR(2,2,0,0)); // ar ar,
      __m256 a_imag = _mm256_shuffle_ps(a,a,_MM_SELECT_FOUR_FOUR(3,3,1,1)); // ai ai
      __m256 tmp = _mm256_shuffle_ps( b,b, _MM_SELECT_FOUR_FOUR(2,3,0,1));
      a_imag = _mm256_mul_ps( a_imag,tmp  );  // (Ai, Ai) * (Bi, Br) = Ai Bi, Ai Br
      return _mm256_maddsub_ps( a_real, b, a_imag ); // Ar Br , Ar Bi   +- Ai Bi             = ArBr-AiBi , ArBi+AiBr
#endif
#if defined (AVX2)
      __m256 a_real = _mm256_moveldup_ps( a ); // Ar Ar
      __m256 a_imag = _mm256_movehdup_ps( a ); // Ai Ai
      a_imag = _mm256_mul_ps( a_imag, _mm256_shuffle_ps( b,b, _MM_SELECT_FOUR_FOUR(2,3,0,1) ));  // (Ai, Ai) * (Bi, Br) = Ai Bi, Ai Br
      return _mm256_fmaddsub_ps( a_real, b, a_imag ); // Ar Br , Ar Bi   +- Ai Bi             = ArBr-AiBi , ArBi+AiBr
#endif
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
#if defined (AVX1) 
      __m256d ymm0,ymm1,ymm2;
      ymm0 = _mm256_shuffle_pd(a,a,0x0); // ymm0 <- ar ar, ar,ar b'00,00
      ymm0 = _mm256_mul_pd(ymm0,b);      // ymm0 <- ar bi, ar br
      ymm1 = _mm256_shuffle_pd(b,b,0x5); // ymm1 <- br,bi  b'01,01 
      ymm2 = _mm256_shuffle_pd(a,a,0xF); // ymm2 <- ai,ai  b'11,11
      ymm1 = _mm256_mul_pd(ymm1,ymm2);   // ymm1 <- br ai, ai bi
      return _mm256_addsub_pd(ymm0,ymm1);
#endif
#if defined (AVXFMA4)
      __m256d a_real = _mm256_shuffle_pd(a,a,0x0);//arar
      __m256d a_imag = _mm256_shuffle_pd(a,a,0xF);//aiai
      a_imag = _mm256_mul_pd( a_imag, _mm256_permute_pd( b, 0x5 ) );  // (Ai, Ai) * (Bi, Br) = Ai Bi, Ai Br
      return _mm256_maddsub_pd( a_real, b, a_imag ); // Ar Br , Ar Bi   +- Ai Bi             = ArBr-AiBi , ArBi+AiBr
#endif
#if defined (AVX2)
      __m256d a_real = _mm256_movedup_pd( a ); // Ar Ar
      __m256d a_imag = _mm256_shuffle_pd(a,a,0xF);//aiai
      a_imag = _mm256_mul_pd( a_imag, _mm256_permute_pd( b, 0x5 ) );  // (Ai, Ai) * (Bi, Br) = Ai Bi, Ai Br
      return _mm256_fmaddsub_pd( a_real, b, a_imag ); // Ar Br , Ar Bi   +- Ai Bi             = ArBr-AiBi , ArBi+AiBr
#endif
    }


  };

#if 0
  struct ComplexDot {

    inline void Prep(__m256 ari,__m256 &air) {
      cdotRIperm(ari,air);
    }
    inline void Mul(__m256 ari,__m256 air,__m256 b,__m256 &riir,__m256 &iirr) {
      riir=air*b;
      iirr=arr*b;
    };
    inline void Madd(__m256 ari,__m256 air,__m256 b,__m256 &riir,__m256 &iirr) {
      mac(riir,air,b);
      mac(iirr,ari,b);
    }
    inline void End(__m256 ari,__m256 &air) {
      //      cdotRI
    }

  };
#endif

  struct Mult{

    inline void mac(__m256 &a, __m256 b, __m256 c){
#if defined (AVX1)
      a= _mm256_add_ps(_mm256_mul_ps(b,c),a);
#endif
#if defined (AVXFMA4)
      a= _mm256_macc_ps(b,c,a);
#endif
#if defined (AVX2)
      a= _mm256_fmadd_ps( b, c, a);
#endif
    }

    inline void mac(__m256d &a, __m256d b, __m256d c){
#if defined (AVX1)
      a= _mm256_add_pd(_mm256_mul_pd(b,c),a);
#endif
#if defined (AVXFMA4)
      a= _mm256_macc_pd(b,c,a);
#endif
#if defined (AVX2)
      a= _mm256_fmadd_pd( b, c, a);
#endif
    }

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
      return _mm256_mullo_epi32(a,b);
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
      return _mm256_xor_pd(_mm256_addsub_pd(_mm256_setzero_pd(),in), _mm256_set1_pd(-0.f));
    }
    // do not define for integer input
  };

  struct TimesMinusI{
    //Complex single
    inline __m256 operator()(__m256 in, __m256 ret){
      __m256 tmp =_mm256_addsub_ps(_mm256_setzero_ps(),in);   // r,-i
      return _mm256_shuffle_ps(tmp,tmp,_MM_SELECT_FOUR_FOUR(2,3,0,1)); //-i,r
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
      __m256 tmp =_mm256_shuffle_ps(in,in,_MM_SELECT_FOUR_FOUR(2,3,0,1)); // i,r
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
  //////////////////////////////////////////////

  struct Permute{

    static inline __m256 Permute0(__m256 in){
      return _mm256_permute2f128_ps(in,in,0x01); //ABCD EFGH -> EFGH ABCD
    };
    static inline __m256 Permute1(__m256 in){
      return _mm256_shuffle_ps(in,in,_MM_SELECT_FOUR_FOUR(1,0,3,2)); //ABCD EFGH -> CDAB GHEF
    };
    static inline __m256 Permute2(__m256 in){
      return _mm256_shuffle_ps(in,in,_MM_SELECT_FOUR_FOUR(2,3,0,1)); //ABCD EFGH -> BADC FEHG
    };
    static inline __m256 Permute3(__m256 in){
      return in;
    };

    static inline __m256d Permute0(__m256d in){
      return _mm256_permute2f128_pd(in,in,0x01); //AB CD -> CD AB
    };
    static inline __m256d Permute1(__m256d in){ //AB CD -> BA DC
      return _mm256_shuffle_pd(in,in,0x5);
    };
    static inline __m256d Permute2(__m256d in){
      return in;
    };
    static inline __m256d Permute3(__m256d in){
      return in;
    };

  };

#if defined (AVX2) || defined (AVXFMA4) 
#define _mm256_alignr_epi32(ret,a,b,n) ret=(__m256) _mm256_alignr_epi8((__m256i)a,(__m256i)b,(n*4)%16)
#define _mm256_alignr_epi64(ret,a,b,n) ret=(__m256d) _mm256_alignr_epi8((__m256i)a,(__m256i)b,(n*8)%16)
#endif

#if defined (AVX1) 

#define _mm256_alignr_epi32(ret,a,b,n) {	\
    __m128 aa, bb;				\
						\
    aa  = _mm256_extractf128_ps(a,1);		\
    bb  = _mm256_extractf128_ps(b,1);		\
    aa  = (__m128)_mm_alignr_epi8((__m128i)aa,(__m128i)bb,(n*4)%16);	\
    ret = _mm256_insertf128_ps(ret,aa,1);	\
						\
    aa  = _mm256_extractf128_ps(a,0);		\
    bb  = _mm256_extractf128_ps(b,0);		\
    aa  = (__m128)_mm_alignr_epi8((__m128i)aa,(__m128i)bb,(n*4)%16);	\
    ret = _mm256_insertf128_ps(ret,aa,0);	\
  }

#define _mm256_alignr_epi64(ret,a,b,n) {	\
    __m128d aa, bb;				\
						\
    aa  = _mm256_extractf128_pd(a,1);		\
    bb  = _mm256_extractf128_pd(b,1);		\
    aa  = (__m128d)_mm_alignr_epi8((__m128i)aa,(__m128i)bb,(n*8)%16);	\
    ret = _mm256_insertf128_pd(ret,aa,1);	\
						\
    aa  = _mm256_extractf128_pd(a,0);		\
    bb  = _mm256_extractf128_pd(b,0);		\
    aa  = (__m128d)_mm_alignr_epi8((__m128i)aa,(__m128i)bb,(n*8)%16);	\
    ret = _mm256_insertf128_pd(ret,aa,0);	\
  }

#endif

    inline std::ostream & operator << (std::ostream& stream, const __m256 a)
    {
      const float *p=(const float *)&a;
      stream<< "{"<<p[0]<<","<<p[1]<<","<<p[2]<<","<<p[3]<<","<<p[4]<<","<<p[5]<<","<<p[6]<<","<<p[7]<<"}";
      return stream;
    };
    inline std::ostream & operator<< (std::ostream& stream, const __m256d a)
    {
      const double *p=(const double *)&a;
      stream<< "{"<<p[0]<<","<<p[1]<<","<<p[2]<<","<<p[3]<<"}";
      return stream;
    };

  struct Rotate{

    static inline __m256 rotate(__m256 in,int n){ 
      switch(n){
      case 0: return tRotate<0>(in);break;
      case 1: return tRotate<1>(in);break;
      case 2: return tRotate<2>(in);break;
      case 3: return tRotate<3>(in);break;
      case 4: return tRotate<4>(in);break;
      case 5: return tRotate<5>(in);break;
      case 6: return tRotate<6>(in);break;
      case 7: return tRotate<7>(in);break;
      default: assert(0);
      }
    }
    static inline __m256d rotate(__m256d in,int n){ 
      switch(n){
      case 0: return tRotate<0>(in);break;
      case 1: return tRotate<1>(in);break;
      case 2: return tRotate<2>(in);break;
      case 3: return tRotate<3>(in);break;
      default: assert(0);
      }
    }
  
    
    template<int n>
    static inline __m256 tRotate(__m256 in){ 
      __m256 tmp = Permute::Permute0(in);
      __m256 ret;
      if ( n > 3 ) { 
	_mm256_alignr_epi32(ret,in,tmp,n);  
      } else {
        _mm256_alignr_epi32(ret,tmp,in,n);          
      }
      //      std::cout << " align epi32 n=" <<n<<" in "<<tmp<<in<<" -> "<< ret <<std::endl;
      return ret;
    };

    template<int n>
    static inline __m256d tRotate(__m256d in){ 
      __m256d tmp = Permute::Permute0(in);
      __m256d ret;
      if ( n > 1 ) {
	_mm256_alignr_epi64(ret,in,tmp,n);          
      } else {
        _mm256_alignr_epi64(ret,tmp,in,n);          
      }
      //      std::cout << " align epi64 n=" <<n<<" in "<<tmp<<in<<" -> "<< ret <<std::endl;
      return ret;
    };

  };



  //Complex float Reduce
  template<>
    inline Grid::ComplexF Reduce<Grid::ComplexF, __m256>::operator()(__m256 in){
    __m256 v1,v2;
    v1=Optimization::Permute::Permute0(in); // avx 256; quad complex single
    v1= _mm256_add_ps(v1,in);
    v2=Optimization::Permute::Permute1(v1); 
    v1 = _mm256_add_ps(v1,v2);
    u256f conv; conv.v = v1;
    return Grid::ComplexF(conv.f[0],conv.f[1]);
  }

  //Real float Reduce
  template<>
  inline Grid::RealF Reduce<Grid::RealF, __m256>::operator()(__m256 in){
    __m256 v1,v2;
    v1 = Optimization::Permute::Permute0(in); // avx 256; octo-double
    v1 = _mm256_add_ps(v1,in);
    v2 = Optimization::Permute::Permute1(v1); 
    v1 = _mm256_add_ps(v1,v2);
    v2 = Optimization::Permute::Permute2(v1); 
    v1 = _mm256_add_ps(v1,v2);
    u256f conv; conv.v=v1;
    return conv.f[0];
  }
  
  
  //Complex double Reduce
  template<>
  inline Grid::ComplexD Reduce<Grid::ComplexD, __m256d>::operator()(__m256d in){
    __m256d v1;
    v1 = Optimization::Permute::Permute0(in); // sse 128; paired complex single
    v1 = _mm256_add_pd(v1,in);
    u256d conv; conv.v = v1;
    return Grid::ComplexD(conv.f[0],conv.f[1]);
  }
  
  //Real double Reduce
  template<>
  inline Grid::RealD Reduce<Grid::RealD, __m256d>::operator()(__m256d in){
    __m256d v1,v2;
    v1 = Optimization::Permute::Permute0(in); // avx 256; quad double
    v1 = _mm256_add_pd(v1,in);
    v2 = Optimization::Permute::Permute1(v1); 
    v1 = _mm256_add_pd(v1,v2);
    u256d conv; conv.v = v1;
    return conv.f[0];
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

  typedef __m256  SIMD_Ftype; // Single precision type
  typedef __m256d SIMD_Dtype; // Double precision type
  typedef __m256i SIMD_Itype; // Integer type

  // prefecthing 
  inline void v_prefetch0(int size, const char *ptr){
    for(int i=0;i<size;i+=64){ //  Define L1 linesize above
      _mm_prefetch(ptr+i+4096,_MM_HINT_T1);
      _mm_prefetch(ptr+i+512,_MM_HINT_T0);
    }
  }
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
  typedef Optimization::Mult        MultSIMD;
  typedef Optimization::MultComplex MultComplexSIMD;
  typedef Optimization::Conj        ConjSIMD;
  typedef Optimization::TimesMinusI TimesMinusISIMD;
  typedef Optimization::TimesI      TimesISIMD;

}
