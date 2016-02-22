    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/simd/Grid_neon.h

    Copyright (C) 2015

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
  @brief Optimization libraries for NEON (ARM) instructions set ARMv8

  Experimental - Using intrinsics - DEVELOPING! 
*/
// Time-stamp: <2015-07-10 17:45:09 neo>
//----------------------------------------------------------------------

#include <arm_neon.h>

// ARMv8 supports double precision

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
    double f[4];
  };
  
  struct Vsplat{
    //Complex float
    inline float32x4_t operator()(float a, float b){
      float tmp[4]={a,b,a,b};
      return vld1q_f32(tmp);
    }
    // Real float
    inline float32x4_t operator()(float a){
      return vld1q_dup_f32(&a);
    }
    //Complex double
    inline float32x4_t operator()(double a, double b){
      float tmp[4]={(float)a,(float)b,(float)a,(float)b};
      return vld1q_f32(tmp);
    }
    //Real double
    inline float32x4_t operator()(double a){
      return vld1q_dup_f32(&a);
    }
    //Integer
    inline uint32x4_t operator()(Integer a){
      return vld1q_dup_u32(&a);
    }
  };

  struct Vstore{
    //Float 
    inline void operator()(float32x4_t a, float* F){
      vst1q_f32(F, a);
    }
    //Double
    inline void operator()(float32x4_t a, double* D){
      vst1q_f32((float*)D, a);
    }
    //Integer
    inline void operator()(uint32x4_t a, Integer* I){
      vst1q_u32(I, a);
    }

  };

  struct Vstream{
    //Float
    inline void operator()(float * a, float32x4_t b){
    
    }
    //Double
    inline void operator()(double * a, float32x4_t b){
  
    }


  };

  struct Vset{
    // Complex float 
    inline float32x4_t operator()(Grid::ComplexF *a){
      float32x4_t foo;
      return foo;
    }
    // Complex double 
    inline float32x4_t operator()(Grid::ComplexD *a){
      float32x4_t foo;
      return foo;
    }
    // Real float 
    inline float32x4_t operator()(float *a){
      float32x4_t foo;
      return foo;
    }
    // Real double
    inline float32x4_t operator()(double *a){
      float32x4_t foo;
      return foo;
    }
    // Integer
    inline uint32x4_t operator()(Integer *a){
      uint32x4_t foo;
      return foo;
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

  struct MultComplex{
    // Complex float
    inline float32x4_t operator()(float32x4_t a, float32x4_t b){
      float32x4_t foo;
      return foo;
    }
    // Complex double
    inline float64x2_t operator()(float64x2_t a, float64x2_t b){
      float32x4_t foo;
      return foo;
    }
  };

  struct Mult{
    // Real float
    inline float32x4_t mac(float32x4_t a, float32x4_t b, float32x4_t c){
      return vaddq_f32(vmulq_f32(b,c),a);
    }
    inline float64x2_t mac(float64x2_t a, float64x2_t b, float64x2_t c){
      return vaddq_f64(vmulq_f64(b,c),a);
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
      return in;
    }
    // Complex double
    //inline float32x4_t operator()(float32x4_t in){
    // return 0;
    //}
    // do not define for integer input
  };

  struct TimesMinusI{
    //Complex single
    inline float32x4_t operator()(float32x4_t in, float32x4_t ret){
      return in;
    }
    //Complex double
    //inline float32x4_t operator()(float32x4_t in, float32x4_t ret){
    //  return in;
    //}


  };

  struct TimesI{
    //Complex single
    inline float32x4_t operator()(float32x4_t in, float32x4_t ret){
      //need shuffle
      return in;
    }
    //Complex double
    //inline float32x4_t operator()(float32x4_t in, float32x4_t ret){
    //  return 0;
    //}
  };

  //////////////////////////////////////////////
  // Some Template specialization
  template < typename vtype > 
    void permute(vtype &a, vtype b, int perm) {

  }; 

  //Complex float Reduce
  template<>
  inline Grid::ComplexF Reduce<Grid::ComplexF, float32x4_t>::operator()(float32x4_t in){
    return 0;
  }
  //Real float Reduce
  template<>
  inline Grid::RealF Reduce<Grid::RealF, float32x4_t>::operator()(float32x4_t in){
    float32x2_t high = vget_high_f32(in);
    float32x2_t low = vget_low_f32(in);
    float32x2_t tmp = vadd_f32(low, high);
    float32x2_t sum = vpadd_f32(tmp, tmp);
    return vget_lane_f32(sum,0);
  }
  
  
  //Complex double Reduce
  template<>
  inline Grid::ComplexD Reduce<Grid::ComplexD, float64x2_t>::operator()(float64x2_t in){
    return 0;
  }
  
  //Real double Reduce
  template<>
  inline Grid::RealD Reduce<Grid::RealD, float64x2_t>::operator()(float64x2_t in){
    float64x2_t sum = vpaddq_f64(in, in);
    return vgetq_lane_f64(sum,0);
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
namespace Grid {

  typedef float32x4_t  SIMD_Ftype; // Single precision type
  typedef float64x2_t  SIMD_Dtype; // Double precision type
  typedef uint32x4_t   SIMD_Itype; // Integer type

  inline void v_prefetch0(int size, const char *ptr){};  // prefetch utilities
  inline void prefetch_HINT_T0(const char *ptr){};


  // Gpermute function
  template < typename VectorSIMD > 
    inline void Gpermute(VectorSIMD &y,const VectorSIMD &b, int perm ) {
    Optimization::permute(y.v,b.v,perm);
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
