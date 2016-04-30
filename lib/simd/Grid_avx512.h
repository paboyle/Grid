    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/simd/Grid_avx512.h

    Copyright (C) 2015

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
/*! @file Grid_knc.h
  @brief Optimization libraries for AVX512 instructions set for KNC

  Using intrinsics
*/
// Time-stamp: <2015-06-09 14:27:28 neo>
//----------------------------------------------------------------------

#include <immintrin.h>


namespace Grid{
namespace Optimization {
  
  struct Vsplat{
    //Complex float
    inline __m512 operator()(float a, float b){
      return _mm512_set_ps(b,a,b,a,b,a,b,a,b,a,b,a,b,a,b,a);
    }
    // Real float
    inline __m512 operator()(float a){
      return _mm512_set1_ps(a);
    }
    //Complex double
    inline __m512d operator()(double a, double b){
      return _mm512_set_pd(b,a,b,a,b,a,b,a);
    }
    //Real double
    inline __m512d operator()(double a){
      return _mm512_set1_pd(a);
    }
    //Integer
    inline __m512i operator()(Integer a){
      return _mm512_set1_epi32(a);
    }
  };

  struct Vstore{
    //Float 
    inline void operator()(__m512 a, float* F){
      _mm512_store_ps(F,a);
    }
    //Double
    inline void operator()(__m512d a, double* D){
      _mm512_store_pd(D,a);
    }
    //Integer
    inline void operator()(__m512i a, Integer* I){
      _mm512_store_si512((__m512i *)I,a);
    }

  };


  struct Vstream{
    //Float
    inline void operator()(float * a, __m512 b){
      //_mm512_stream_ps(a,b);
      _mm512_store_ps(a,b);
    }
    //Double
    inline void operator()(double * a, __m512d b){
      //_mm512_stream_pd(a,b);
      _mm512_store_pd(a,b);
    }

  };



  struct Vset{
    // Complex float 
    inline __m512 operator()(Grid::ComplexF *a){
      return _mm512_set_ps(a[7].imag(),a[7].real(),a[6].imag(),a[6].real(),
			   a[5].imag(),a[5].real(),a[4].imag(),a[4].real(),
			   a[3].imag(),a[3].real(),a[2].imag(),a[2].real(),
			   a[1].imag(),a[1].real(),a[0].imag(),a[0].real());
    }
    // Complex double 
    inline __m512d operator()(Grid::ComplexD *a){
      return _mm512_set_pd(a[3].imag(),a[3].real(),a[2].imag(),a[2].real(),
			   a[1].imag(),a[1].real(),a[0].imag(),a[0].real());
    }
    // Real float 
    inline __m512 operator()(float *a){
      return _mm512_set_ps( a[15],a[14],a[13],a[12],a[11],a[10],a[9],a[8],
			    a[7],a[6],a[5],a[4],a[3],a[2],a[1],a[0]);
    }
    // Real double
    inline __m512d operator()(double *a){
      return _mm512_set_pd(a[7],a[6],a[5],a[4],a[3],a[2],a[1],a[0]);
    }
    // Integer
    inline __m512i operator()(Integer *a){
      return _mm512_set_epi32( a[15],a[14],a[13],a[12],a[11],a[10],a[9],a[8],
			       a[7],a[6],a[5],a[4],a[3],a[2],a[1],a[0]);
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
    inline __m512 operator()(__m512 a, __m512 b){
      return _mm512_add_ps(a,b);
    }
    //Complex/Real double
    inline __m512d operator()(__m512d a, __m512d b){
      return _mm512_add_pd(a,b);
    }
    //Integer
    inline __m512i operator()(__m512i a, __m512i b){
      return _mm512_add_epi32(a,b);
    }
  };

  struct Sub{
    //Complex/Real float
    inline __m512 operator()(__m512 a, __m512 b){
      return _mm512_sub_ps(a,b);
    }
    //Complex/Real double
    inline __m512d operator()(__m512d a, __m512d b){
      return _mm512_sub_pd(a,b);
    }
    //Integer
    inline __m512i operator()(__m512i a, __m512i b){
      return _mm512_sub_epi32(a,b);
    }
  };

  // Note, we can beat the shuf overhead in chain with two temporaries
  // Ar Ai , Br Bi,  Ai Ar  // one shuf
  //tmpr Ar Br,  Ai Bi    // Mul/Mac/Mac
  //tmpi Br Ai,  Bi Ar    // Mul/Mac/Mac
  // add tmpi,shuf(tmpi)
  // sub tmpr,shuf(tmpi)
  // shuf(tmpr,tmpi).    // Could drop/trade for write mask

  // Gives
  //  2mul,4 mac +add+sub = 8 flop type insns
  //  3shuf + 2 (+shuf)   = 5/6 simd perm and 1/2 the load.

  struct MultComplex{
    // Complex float
    inline __m512 operator()(__m512 a, __m512 b){
      // dup, dup, perm, mul, madd
      __m512 a_real = _mm512_moveldup_ps( a ); // Ar Ar
      __m512 a_imag = _mm512_movehdup_ps( a ); // Ai Ai
      a_imag = _mm512_mul_ps( a_imag, _mm512_permute_ps( b, 0xB1 ) );  // (Ai, Ai) * (Bi, Br) = Ai Bi, Ai Br
      return _mm512_fmaddsub_ps( a_real, b, a_imag ); // Ar Br , Ar Bi   +- Ai Bi             = ArBr-AiBi , ArBi+AiBr
    }
    // Complex double
    inline __m512d operator()(__m512d a, __m512d b){
      __m512d a_real = _mm512_shuffle_pd( a, a, 0x00 );
      __m512d a_imag = _mm512_shuffle_pd( a, a, 0xFF );
      a_imag = _mm512_mul_pd( a_imag, _mm512_permute_pd( b, 0x55 ) ); 
      return _mm512_fmaddsub_pd( a_real, b, a_imag );
    }
  };
  
  struct Mult{

    inline void mac(__m512 &a, __m512 b, __m512 c){         
       a= _mm512_fmadd_ps( b, c, a);                         
    }

    inline void mac(__m512d &a, __m512d b, __m512d c){
      a= _mm512_fmadd_pd( b, c, a);                   
    }                                             

    // Real float
    inline __m512 operator()(__m512 a, __m512 b){
      return _mm512_mul_ps(a,b);
    }
    // Real double
    inline __m512d operator()(__m512d a, __m512d b){
      return _mm512_mul_pd(a,b);
    }
    // Integer
    inline __m512i operator()(__m512i a, __m512i b){
      return _mm512_mullo_epi32(a,b);
    }
  };


  struct Conj{
    // Complex single
    inline __m512 operator()(__m512 in){
      return _mm512_mask_sub_ps(in,0xaaaa,_mm512_setzero_ps(),in); // Zero out 0+real 0-imag  
    }
    // Complex double
    inline __m512d operator()(__m512d in){
      return _mm512_mask_sub_pd(in, 0xaa,_mm512_setzero_pd(), in);
    }
    // do not define for integer input
  };

  struct TimesMinusI{
    //Complex single
    inline __m512 operator()(__m512 in, __m512 ret){
      //__m512 tmp = _mm512_mask_sub_ps(in,0xaaaa,_mm512_setzero_ps(),in); // real -imag 
      //return _mm512_shuffle_ps(tmp,tmp,_MM_SELECT_FOUR_FOUR(2,3,1,0));   // 0x4E??
      __m512 tmp = _mm512_shuffle_ps(in,in,_MM_SELECT_FOUR_FOUR(2,3,0,1));
      return _mm512_mask_sub_ps(tmp,0xaaaa,_mm512_setzero_ps(),tmp);
    }
    //Complex double
    inline __m512d operator()(__m512d in, __m512d ret){
      //__m512d tmp = _mm512_mask_sub_pd(in,0xaa,_mm512_setzero_pd(),in); // real -imag 
      //return _mm512_shuffle_pd(tmp,tmp,0x55);
      __m512d tmp = _mm512_shuffle_pd(in,in,0x55);
      return _mm512_mask_sub_pd(tmp,0xaa,_mm512_setzero_pd(),tmp);
    } 
  };

  struct TimesI{
    //Complex single
    inline __m512 operator()(__m512 in, __m512 ret){
      __m512 tmp = _mm512_shuffle_ps(in,in,_MM_SELECT_FOUR_FOUR(2,3,0,1));
      return _mm512_mask_sub_ps(tmp,0x5555,_mm512_setzero_ps(),tmp); 
    }
    //Complex double
    inline __m512d operator()(__m512d in, __m512d ret){
      __m512d tmp = _mm512_shuffle_pd(in,in,0x55);
      return _mm512_mask_sub_pd(tmp,0x55,_mm512_setzero_pd(),tmp); 
    }


  };


  
  // Gpermute utilities consider coalescing into 1 Gpermute
  struct Permute{
    
    static inline __m512 Permute0(__m512 in){
      return _mm512_shuffle_f32x4(in,in,_MM_SELECT_FOUR_FOUR(1,0,3,2));
    };
    static inline __m512 Permute1(__m512 in){
      return _mm512_shuffle_f32x4(in,in,_MM_SELECT_FOUR_FOUR(2,3,0,1));
    };
    static inline __m512 Permute2(__m512 in){
      return _mm512_shuffle_ps(in,in,_MM_SELECT_FOUR_FOUR(1,0,3,2));
    };
    static inline __m512 Permute3(__m512 in){
      return _mm512_shuffle_ps(in,in,_MM_SELECT_FOUR_FOUR(2,3,0,1));
    };

    static inline __m512d Permute0(__m512d in){
      return _mm512_shuffle_f64x2(in,in,_MM_SELECT_FOUR_FOUR(1,0,3,2));
    };
    static inline __m512d Permute1(__m512d in){
      return _mm512_shuffle_f64x2(in,in,_MM_SELECT_FOUR_FOUR(2,3,0,1));
    };
    static inline __m512d Permute2(__m512d in){
      return _mm512_shuffle_pd(in,in,0x55);
    };
    static inline __m512d Permute3(__m512d in){
      return in;
    };

  };


  struct Rotate{

    static inline __m512 rotate(__m512 in,int n){ 
      switch(n){
      case 0: return tRotate<0>(in);break;
      case 1: return tRotate<1>(in);break;
      case 2: return tRotate<2>(in);break;
      case 3: return tRotate<3>(in);break;
      case 4: return tRotate<4>(in);break;
      case 5: return tRotate<5>(in);break;
      case 6: return tRotate<6>(in);break;
      case 7: return tRotate<7>(in);break;

      case 8 : return tRotate<8>(in);break;
      case 9 : return tRotate<9>(in);break;
      case 10: return tRotate<10>(in);break;
      case 11: return tRotate<11>(in);break;
      case 12: return tRotate<12>(in);break;
      case 13: return tRotate<13>(in);break;
      case 14: return tRotate<14>(in);break;
      case 15: return tRotate<15>(in);break;
      default: assert(0);
      }
    }
    static inline __m512d rotate(__m512d in,int n){ 
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

    template<int n> static inline __m512 tRotate(__m512 in){ 
      return (__m512)_mm512_alignr_epi32((__m512i)in,(__m512i)in,n);          
    };

    template<int n> static inline __m512d tRotate(__m512d in){ 
      return (__m512d)_mm512_alignr_epi64((__m512i)in,(__m512i)in,n);          
    };

  };

  //////////////////////////////////////////////
  // Some Template specialization
  
  //Complex float Reduce
  template<>
  inline Grid::ComplexF Reduce<Grid::ComplexF, __m512>::operator()(__m512 in){
    return Grid::ComplexF(_mm512_mask_reduce_add_ps(0x5555, in),_mm512_mask_reduce_add_ps(0xAAAA, in));
  }
  //Real float Reduce
  template<>
  inline Grid::RealF Reduce<Grid::RealF, __m512>::operator()(__m512 in){
    return _mm512_reduce_add_ps(in);
  }
  
  
  //Complex double Reduce
  template<>
  inline Grid::ComplexD Reduce<Grid::ComplexD, __m512d>::operator()(__m512d in){
    return Grid::ComplexD(_mm512_mask_reduce_add_pd(0x55, in),_mm512_mask_reduce_add_pd(0xAA, in));
  }
  
  //Real double Reduce
  template<>
  inline Grid::RealD Reduce<Grid::RealD, __m512d>::operator()(__m512d in){
    return _mm512_reduce_add_pd(in);
  }

  //Integer Reduce
  template<>
  inline Integer Reduce<Integer, __m512i>::operator()(__m512i in){
    // FIXME unimplemented
    printf("Reduce : Missing integer implementation -> FIX\n");
    assert(0);
  }
  
  
}

//////////////////////////////////////////////////////////////////////////////////////
// Here assign types 

  typedef __m512 SIMD_Ftype;  // Single precision type
  typedef __m512d SIMD_Dtype; // Double precision type
  typedef __m512i SIMD_Itype; // Integer type

  // prefecth
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
