/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./lib/simd/Grid_vector_types.h

    Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
Author: Guido Cossu <cossu@iroiro-pc.kek.jp>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: neo <cossu@post.kek.jp>
Author: paboyle <paboyle@ph.ed.ac.uk>
Author: Michael Marshall <michael.marshall@ed.ac.au>

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
//---------------------------------------------------------------------------
/*! @file Grid_vector_types.h
  @brief Defines templated class Grid_simd to deal with inner vector types
*/
// Time-stamp: <2015-07-10 17:45:33 neo>
//---------------------------------------------------------------------------
#ifndef GRID_VECTOR_TYPES
#define GRID_VECTOR_TYPES

// PAB - Lifted and adapted from Eigen, which is GPL V2
struct Grid_half {
  accelerator Grid_half(){}
  accelerator Grid_half(uint16_t raw) : x(raw) {}
  uint16_t x;
};
union FP32 {
  unsigned int u;
  float f;
};
accelerator_inline float sfw_half_to_float(Grid_half h) {
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
accelerator_inline Grid_half sfw_float_to_half(float ff) {
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


#ifdef GPU_VEC
#include "Grid_gpu_vec.h"
#endif

#ifdef GPU_RRII
#include "Grid_gpu_rrii.h"
#endif

#ifdef GEN
  #if defined(A64FX) || defined(A64FXFIXEDSIZE) // breakout A64FX SVE ACLE here
    #include <arm_sve.h>
    #if defined(A64FX) // VLA
      #pragma message("building A64FX / SVE ACLE VLA")
      #if defined(ARMCLANGCOMPAT)
        #pragma message("applying data types patch")
      #endif
      #include "Grid_a64fx-2.h"
    #endif
    #if defined(A64FXFIXEDSIZE) // fixed size data types
      #pragma message("building for A64FX / SVE ACLE fixed size")
      #include "Grid_a64fx-fixedsize.h"
    #endif
  #else
    #include "Grid_generic.h"
  #endif
#endif

#ifdef A64FX
  #include <arm_sve.h>
  #ifdef __ARM_FEATURE_SVE_BITS
    //#pragma message("building A64FX SVE VLS")
    #include "Grid_a64fx-fixedsize.h"
  #else
    #pragma message("building A64FX SVE VLA")
    #if defined(ARMCLANGCOMPAT)
      #pragma message("applying data types patch")
    #endif
    #include "Grid_a64fx-2.h"
  #endif
#endif

#ifdef SSE4
#include "Grid_sse4.h"
#endif
#if defined(AVX1) || defined (AVXFMA) || defined(AVX2) || defined(AVXFMA4)
#include "Grid_avx.h"
#endif
#if defined AVX512
#include "Grid_avx512.h"
#endif
#if defined IMCI
#include "Grid_imci.h"
#endif
#ifdef NEONV8
#include "Grid_neon.h"
#endif
#if defined QPX
#include "Grid_qpx.h"
#endif

NAMESPACE_BEGIN(Grid);


//////////////////////////////////////
// To take the floating point type of real/complex type
//////////////////////////////////////
template <typename T>
struct RealPart {
  typedef T type;
};
template <typename T>
struct RealPart<complex<T> > {
  typedef T type;
};

#include <type_traits>

//////////////////////////////////////
// demote a vector to real type
//////////////////////////////////////
// type alias used to simplify the syntax of std::enable_if
template <typename T> using Invoke = typename T::type;
template <typename Condition, typename ReturnType = void> using EnableIf    = Invoke<std::enable_if<Condition::value, ReturnType> >;
template <typename Condition, typename ReturnType = void> using NotEnableIf = Invoke<std::enable_if<!Condition::value, ReturnType> >;

////////////////////////////////////////////////////////
// Check for complexity with type traits
template <typename T> struct is_complex : public std::false_type {};
template <> struct is_complex<ComplexD> : public std::true_type {};
template <> struct is_complex<ComplexF> : public std::true_type {};

template <typename T> struct is_ComplexD : public std::false_type {};
template <> struct is_ComplexD<ComplexD> : public std::true_type {};

template <typename T> struct is_ComplexF : public std::false_type {};
template <> struct is_ComplexF<ComplexF> : public std::true_type {};

template<typename T, typename V=void> struct is_real : public std::false_type {};
template<typename T> struct is_real<T, typename std::enable_if<std::is_floating_point<T>::value,
  void>::type> : public std::true_type {};

template<typename T, typename V=void> struct is_integer : public std::false_type {};
template<typename T> struct is_integer<T, typename std::enable_if<std::is_integral<T>::value,
  void>::type> : public std::true_type {};

template <typename T>              using IfReal    = Invoke<std::enable_if<is_real<T>::value, int> >;
template <typename T>              using IfComplex = Invoke<std::enable_if<is_complex<T>::value, int> >;
template <typename T>              using IfInteger = Invoke<std::enable_if<is_integer<T>::value, int> >;
template <typename T1,typename T2> using IfSame    = Invoke<std::enable_if<std::is_same<T1,T2>::value, int> >;

template <typename T>              using IfNotReal    = Invoke<std::enable_if<!is_real<T>::value, int> >;
template <typename T>              using IfNotComplex = Invoke<std::enable_if<!is_complex<T>::value, int> >;
template <typename T>              using IfNotInteger = Invoke<std::enable_if<!is_integer<T>::value, int> >;
template <typename T1,typename T2> using IfNotSame    = Invoke<std::enable_if<!std::is_same<T1,T2>::value, int> >;

////////////////////////////////////////////////////////
// Define the operation templates functors
// general forms to allow for vsplat syntax
// need explicit declaration of types when used since
// clang cannot automatically determine the output type sometimes
template <class Out, class Input1, class Input2, class Input3, class Operation>
Out accelerator_inline  trinary(Input1 src_1, Input2 src_2, Input3 src_3, Operation op) {
  return op(src_1, src_2, src_3);
}
template <class Out, class Input1, class Input2, class Operation>
Out accelerator_inline binary(Input1 src_1, Input2 src_2, Operation op) {
  return op(src_1, src_2);
}
template <class Out, class Input, class Operation>
Out accelerator_inline  unary(Input src, Operation op) {
  return op(src);
}
///////////////////////////////////////////////

/*
  @brief Grid_simd class for the SIMD vector type operations
*/
template <class Scalar_type, class Vector_type>
class Grid_simd {
public:
  typedef typename RealPart<Scalar_type>::type Real;
  typedef Vector_type vector_type;
  typedef Scalar_type scalar_type;

  /*
  typedef union conv_t_union {
    Vector_type v;
    Scalar_type s[sizeof(Vector_type) / sizeof(Scalar_type)];
    accelerator_inline conv_t_union(){};
  } conv_t;
  */
  
  Vector_type v;

  static accelerator_inline constexpr int Nsimd(void) {
    static_assert( (sizeof(Vector_type) / sizeof(Scalar_type) >= 1), " size mismatch " );
    return sizeof(Vector_type) / sizeof(Scalar_type);
  }

  #ifdef ARMCLANGCOMPAT
    template <class S = Scalar_type>
    accelerator_inline Grid_simd &operator=(const Grid_simd<typename std::enable_if<!is_complex<S>::value, S>::type, Vector_type> &&rhs) {
      //v = rhs.v;
      svst1(svptrue_b8(), (Scalar_type*)this, svld1(svptrue_b8(), (Scalar_type*)&(rhs.v)));
      return *this;
    };

    template <class S = Scalar_type>
    accelerator_inline Grid_simd &operator=(const Grid_simd<typename std::enable_if<!is_complex<S>::value, S>::type, Vector_type> &rhs) {
      //v = rhs.v;
      svst1(svptrue_b8(), (Scalar_type*)this, svld1(svptrue_b8(), (Scalar_type*)&(rhs.v)));
      return *this;
    };

    /*
    template <class S = Scalar_type>
    accelerator_inline Grid_simd &operator=(const Grid_simd<typename std::enable_if<is_complex<S>::value, S>::type, Vector_type> &&rhs) {
      //v = rhs.v;
      svst1(svptrue_b8(), (int8_t*)this, svld1(svptrue_b8(), (int8_t*)&(rhs.v)));
      return *this;
    };

    template <class S = Scalar_type>
    accelerator_inline Grid_simd &operator=(const Grid_simd<typename std::enable_if<is_complex<S>::value, S>::type, Vector_type> &rhs) {
      //v = rhs.v;
      svst1(svptrue_b8(), (int8_t*)this, svld1(svptrue_b8(), (int8_t*)&(rhs.v)));
      return *this;
    };
    */

    // ComplexF
    template <class S = Scalar_type>
    accelerator_inline Grid_simd &operator=(const Grid_simd<typename std::enable_if<is_ComplexF<S>::value, S>::type, Vector_type> &&rhs) {
      //v = rhs.v;
      svst1(svptrue_b32(), (float*)this, svld1(svptrue_b32(), (float*)&(rhs.v)));
      return *this;
    };

    template <class S = Scalar_type>
    accelerator_inline Grid_simd &operator=(const Grid_simd<typename std::enable_if<is_ComplexF<S>::value, S>::type, Vector_type> &rhs) {
      //v = rhs.v;
      svst1(svptrue_b32(), (float*)this, svld1(svptrue_b32(), (float*)&(rhs.v)));
      return *this;
    };

    // ComplexD
    template <class S = Scalar_type>
    accelerator_inline Grid_simd &operator=(const Grid_simd<typename std::enable_if<is_ComplexD<S>::value, S>::type, Vector_type> &&rhs) {
      //v = rhs.v;
      svst1(svptrue_b64(), (double*)this, svld1(svptrue_b64(), (double*)&(rhs.v)));
      return *this;
    };

    template <class S = Scalar_type>
    accelerator_inline Grid_simd &operator=(const Grid_simd<typename std::enable_if<is_ComplexD<S>::value, S>::type, Vector_type> &rhs) {
      //v = rhs.v;
      svst1(svptrue_b64(), (double*)this, svld1(svptrue_b64(), (double*)&(rhs.v)));
      return *this;
    };

  #else

  accelerator_inline Grid_simd &operator=(const Grid_simd &&rhs) {
    v = rhs.v;
    return *this;
  };
  accelerator_inline Grid_simd &operator=(const Grid_simd &rhs) {
    v = rhs.v;
    return *this;
  };  // faster than not declaring it and leaving to the compiler

  #endif

  accelerator Grid_simd() = default;

  #ifdef ARMCLANGCOMPAT
    template <class S = Scalar_type>
    accelerator_inline Grid_simd(const Grid_simd<typename std::enable_if<!is_complex<S>::value, S>::type, Vector_type> &rhs) { this->operator=(rhs); }
    template <class S = Scalar_type>
    accelerator_inline Grid_simd(const Grid_simd<typename std::enable_if<!is_complex<S>::value, S>::type, Vector_type> &&rhs) { this->operator=(rhs); }
    template <class S = Scalar_type>
    accelerator_inline Grid_simd(const Grid_simd<typename std::enable_if<is_complex<S>::value, S>::type, Vector_type> &rhs) { this->operator=(rhs); }
    template <class S = Scalar_type>
    accelerator_inline Grid_simd(const Grid_simd<typename std::enable_if<is_complex<S>::value, S>::type, Vector_type> &&rhs) { this->operator=(rhs); }
  #else
    accelerator_inline Grid_simd(const Grid_simd &rhs) : v(rhs.v){};  // compiles in movaps
    accelerator_inline Grid_simd(const Grid_simd &&rhs) : v(rhs.v){};
  #endif
  accelerator_inline Grid_simd(const Real a) { vsplat(*this, Scalar_type(a)); };
  // Enable if complex type
  template <typename S = Scalar_type> accelerator_inline
  Grid_simd(const typename std::enable_if<is_complex<S>::value, S>::type a) {
      vsplat(*this, a);
  };

  /////////////////////////////
  // Constructors
  /////////////////////////////
  accelerator_inline Grid_simd &  operator=(const Zero &z) {
    vzero(*this);
    return (*this);
  }



  ///////////////////////////////////////////////
  // mac, mult, sub, add, adj
  ///////////////////////////////////////////////

  // FIXME -- alias this to an accelerator_inline MAC struct.

  #if defined(A64FX) || defined(A64FXFIXEDSIZE)
  friend accelerator_inline void mac(Grid_simd *__restrict__ y,
				     const Grid_simd *__restrict__ a,
				     const Grid_simd *__restrict__ x) {
    *y = fxmac((*a), (*x), (*y));
  };
  #else
  friend accelerator_inline void mac(Grid_simd *__restrict__ y,
				     const Grid_simd *__restrict__ a,
				     const Grid_simd *__restrict__ x) {
    *y = (*a) * (*x) + (*y);
  };
  #endif

  friend accelerator_inline void mult(Grid_simd *__restrict__ y,
				      const Grid_simd *__restrict__ l,
				      const Grid_simd *__restrict__ r) {
    *y = (*l) * (*r);
  }

  friend accelerator_inline void sub(Grid_simd *__restrict__ y,
				     const Grid_simd *__restrict__ l,
				     const Grid_simd *__restrict__ r) {
    *y = (*l) - (*r);
  }
  friend accelerator_inline void add(Grid_simd *__restrict__ y,
				     const Grid_simd *__restrict__ l,
				     const Grid_simd *__restrict__ r) {
    *y = (*l) + (*r);
  }
  friend accelerator_inline void mac(Grid_simd *__restrict__ y,
				     const Scalar_type *__restrict__ a,
				     const Grid_simd *__restrict__ x) {
    *y = (*a) * (*x) + (*y);
  };
  friend accelerator_inline void mult(Grid_simd *__restrict__ y,
				      const Scalar_type *__restrict__ l,
				      const Grid_simd *__restrict__ r) {
    *y = (*l) * (*r);
  }
  friend accelerator_inline void sub(Grid_simd *__restrict__ y,
				     const Scalar_type *__restrict__ l,
				     const Grid_simd *__restrict__ r) {
    *y = (*l) - (*r);
  }
  friend accelerator_inline void add(Grid_simd *__restrict__ y,
				     const Scalar_type *__restrict__ l,
				     const Grid_simd *__restrict__ r) {
    *y = (*l) + (*r);
  }

  friend accelerator_inline void mac(Grid_simd *__restrict__ y,
				     const Grid_simd *__restrict__ a,
				     const Scalar_type *__restrict__ x) {
    *y = (*a) * (*x) + (*y);
  };
  friend accelerator_inline void mult(Grid_simd *__restrict__ y,
				      const Grid_simd *__restrict__ l,
				      const Scalar_type *__restrict__ r) {
    *y = (*l) * (*r);
  }
  friend accelerator_inline void sub(Grid_simd *__restrict__ y,
				     const Grid_simd *__restrict__ l,
				     const Scalar_type *__restrict__ r) {
    *y = (*l) - (*r);
  }
  friend accelerator_inline void add(Grid_simd *__restrict__ y,
				     const Grid_simd *__restrict__ l,
				     const Scalar_type *__restrict__ r) {
    *y = (*l) + (*r);
  }

  ////////////////////////////////////////////////////////////////////////
  // FIXME:  gonna remove these load/store, get, set, prefetch
  ////////////////////////////////////////////////////////////////////////
  friend accelerator_inline void vset(Grid_simd &ret, Scalar_type *a) {
    ret.v = unary<Vector_type>(a, VsetSIMD());
  }

  ///////////////////////
  // Vstore
  ///////////////////////
  friend accelerator_inline void vstore(const Grid_simd &ret, Scalar_type *a) {
    binary<void>(ret.v, (Real *) a, VstoreSIMD());
  }

  ///////////////////////
  // Vprefetch
  ///////////////////////
  friend accelerator_inline void vprefetch(const Grid_simd &v) {
    prefetch_HINT_T0((const char *)&v.v);
  }

  ///////////////////////
  // Reduce
  ///////////////////////
  friend accelerator_inline Scalar_type Reduce(const Grid_simd &in) {
    return unary<Scalar_type>(in.v, ReduceSIMD<Scalar_type, Vector_type>());
  }

  ////////////////////////////
  // operator scalar * simd
  ////////////////////////////
  friend accelerator_inline Grid_simd operator*(const Scalar_type &a, Grid_simd b) {
    Grid_simd va;
    vsplat(va, a);
    return va * b;
  }
  friend accelerator_inline Grid_simd operator*(Grid_simd b, const Scalar_type &a) {
    return a * b;
  }

  //////////////////////////////////
  // Divides
  //////////////////////////////////
  friend accelerator_inline Grid_simd operator/(const Scalar_type &a, Grid_simd b) {
    Grid_simd va;
    vsplat(va, a);
    return va / b;
  }
  friend accelerator_inline Grid_simd operator/(Grid_simd b, const Scalar_type &a) {
    Grid_simd va;
    vsplat(va, a);
    return b / a;
  }

  ///////////////////////
  // Unary negation
  ///////////////////////
  friend accelerator_inline Grid_simd operator-(const Grid_simd &r) {
    Grid_simd ret;
    vzero(ret);
    ret = ret - r;
    return ret;
  }
  // *=,+=,-= operators
  accelerator_inline Grid_simd &operator*=(const Grid_simd &r) {
    *this = (*this) * r;
    return *this;
    // return (*this)*r; ?
  }
  accelerator_inline Grid_simd &operator+=(const Grid_simd &r) {
    *this = *this + r;
    return *this;
  }
  accelerator_inline Grid_simd &operator-=(const Grid_simd &r) {
    *this = *this - r;
    return *this;
  }

  ///////////////////////////////////////
  // Not all functions are supported
  // through SIMD and must breakout to
  // scalar type and back again. This
  // provides support
  ///////////////////////////////////////

  template <class functor>
  friend accelerator_inline Grid_simd SimdApply(const functor &func, const Grid_simd &v) {
    Grid_simd ret;
    Grid_simd::scalar_type s;

    for (int i = 0; i < Nsimd(); i++) {
      s = v.getlane(i);
      s = func(s);
      ret.putlane(s,i);
    }
    return ret;
  }
  template <class functor>
  friend accelerator_inline Grid_simd SimdApplyBinop(const functor &func,
                                         const Grid_simd &x,
                                         const Grid_simd &y) {
    Grid_simd ret;
    Grid_simd::scalar_type sx,sy;

    for (int i = 0; i < Nsimd(); i++) {
      sx = x.getlane(i);
      sy = y.getlane(i);
      sx = func(sx,sy);
      ret.putlane(sx,i);
    }
    return ret;
  }
  ///////////////////////
  // Exchange
  // Al Ah , Bl Bh -> Al Bl Ah,Bh
  ///////////////////////
  friend accelerator_inline void exchange(Grid_simd &out1,Grid_simd &out2,Grid_simd in1,Grid_simd in2,int n)
  {
    if       (n==3) {
      Optimization::Exchange::Exchange3(out1.v,out2.v,in1.v,in2.v);
    } else if(n==2) {
      Optimization::Exchange::Exchange2(out1.v,out2.v,in1.v,in2.v);
    } else if(n==1) {
      Optimization::Exchange::Exchange1(out1.v,out2.v,in1.v,in2.v);
    } else if(n==0) {
      Optimization::Exchange::Exchange0(out1.v,out2.v,in1.v,in2.v);
    }
  }
  friend accelerator_inline void exchange0(Grid_simd &out1,Grid_simd &out2,Grid_simd in1,Grid_simd in2){
    Optimization::Exchange::Exchange0(out1.v,out2.v,in1.v,in2.v);
  }
  friend accelerator_inline void exchange1(Grid_simd &out1,Grid_simd &out2,Grid_simd in1,Grid_simd in2){
    Optimization::Exchange::Exchange1(out1.v,out2.v,in1.v,in2.v);
  }
  friend accelerator_inline void exchange2(Grid_simd &out1,Grid_simd &out2,Grid_simd in1,Grid_simd in2){
    Optimization::Exchange::Exchange2(out1.v,out2.v,in1.v,in2.v);
  }
  friend accelerator_inline void exchange3(Grid_simd &out1,Grid_simd &out2,Grid_simd in1,Grid_simd in2){
    Optimization::Exchange::Exchange3(out1.v,out2.v,in1.v,in2.v);
  }
  ////////////////////////////////////////////////////////////////////
  // General permute; assumes vector length is same across
  // all subtypes; may not be a good assumption, but could
  // add the vector width as a template param for BG/Q for example
  ////////////////////////////////////////////////////////////////////
  friend accelerator_inline void permute0(Grid_simd &y, Grid_simd b) {
    y.v = Optimization::Permute::Permute0(b.v);
  }
  friend accelerator_inline void permute1(Grid_simd &y, Grid_simd b) {
    y.v = Optimization::Permute::Permute1(b.v);
  }
  friend accelerator_inline void permute2(Grid_simd &y, Grid_simd b) {
    y.v = Optimization::Permute::Permute2(b.v);
  }
  friend accelerator_inline void permute3(Grid_simd &y, Grid_simd b) {
    y.v = Optimization::Permute::Permute3(b.v);
  }
  friend accelerator_inline void permute(Grid_simd &y, Grid_simd b, int perm) {
    if (perm & RotateBit) {
      int dist = perm & 0xF;
      y = rotate(b, dist);
      return;
    }
    else if(perm==3) permute3(y, b);
    else if(perm==2) permute2(y, b);
    else if(perm==1) permute1(y, b);
    else if(perm==0) permute0(y, b);
  }

  ///////////////////////////////
  // Getting single lanes
  ///////////////////////////////
#ifdef GPU_RRII
  template <class S = Scalar_type,IfComplex<S> = 0>
  accelerator_inline Scalar_type getlane(int lane) const {
    return Scalar_type(v.rrrr[lane],v.iiii[lane]);
  }
  template <class S = Scalar_type,IfComplex<S> = 0>
  accelerator_inline void putlane(const Scalar_type &_S, int lane){
    v.rrrr[lane] = real(_S);
    v.iiii[lane] = imag(_S);
  }
  template <class S = Scalar_type,IfNotComplex<S> = 0>
  accelerator_inline Scalar_type getlane(int lane) const {
    return ((S*)&v)[lane];
  }
  template <class S = Scalar_type,IfNotComplex<S> = 0>
  accelerator_inline void putlane(const S &_S, int lane){
    ((Scalar_type*)&v)[lane] = _S;
  }
#else // Can pun to an array of complex
  accelerator_inline Scalar_type getlane(int lane) const {
    return ((Scalar_type*)&v)[lane];
  }
  accelerator_inline void putlane(const Scalar_type &S, int lane){
    ((Scalar_type*)&v)[lane] = S;
  }
#endif

};  // end of Grid_simd class definition


///////////////////////////////
// Define available types
///////////////////////////////

typedef Grid_simd<float  , SIMD_Ftype> vRealF;
typedef Grid_simd<double , SIMD_Dtype> vRealD;
typedef Grid_simd<Integer, SIMD_Itype> vInteger;
typedef Grid_simd<uint16_t,SIMD_Htype> vRealH;

#if defined(GPU_VEC) || defined(GPU_RRII)
typedef Grid_simd<complex<uint16_t>, SIMD_CHtype> vComplexH;
typedef Grid_simd<complex<float>   , SIMD_CFtype> vComplexF;
typedef Grid_simd<complex<double>  , SIMD_CDtype> vComplexD;
#else
typedef Grid_simd<complex<uint16_t>, SIMD_Htype> vComplexH;
typedef Grid_simd<complex<float>   , SIMD_Ftype> vComplexF;
typedef Grid_simd<complex<double>  , SIMD_Dtype> vComplexD;
#endif

/////////////////////////////////////////
// Pointer type to use on extractLane
/////////////////////////////////////////
template<class _scalar> class ExtractTypeMap                      { public: typedef _scalar extract_type;};
#ifdef GPU_VEC
template<>              class ExtractTypeMap< complex<uint16_t> > { public: typedef   half2 extract_type;};
template<>              class ExtractTypeMap< complex<   float> > { public: typedef  float2 extract_type;};
template<>              class ExtractTypeMap< complex<  double> > { public: typedef double2 extract_type;};
#endif

/////////////////////////////////////////
// Permute
/////////////////////////////////////////

accelerator_inline void permute(ComplexD &y,ComplexD b, int perm) {  y=b; }
accelerator_inline void permute(ComplexF &y,ComplexF b, int perm) {  y=b; }
accelerator_inline void permute(RealD &y,RealD b, int perm) {  y=b; }
accelerator_inline void permute(RealF &y,RealF b, int perm) {  y=b; }

////////////////////////////////////////////////////////////////////
// General rotate
////////////////////////////////////////////////////////////////////
template <class S, class V, IfNotComplex<S> = 0>
accelerator_inline Grid_simd<S, V> rotate(Grid_simd<S, V> b, int nrot) {
  nrot = nrot % Grid_simd<S, V>::Nsimd();
  Grid_simd<S, V> ret;
  ret.v = Optimization::Rotate::rotate(b.v, nrot);
  return ret;
}
template <class S, class V, IfComplex<S> = 0>
accelerator_inline Grid_simd<S, V> rotate(Grid_simd<S, V> b, int nrot) {
  nrot = nrot % Grid_simd<S, V>::Nsimd();
  Grid_simd<S, V> ret;
  ret.v = Optimization::Rotate::rotate(b.v, 2 * nrot);
  return ret;
}
template <class S, class V, IfNotComplex<S> =0>
accelerator_inline void rotate( Grid_simd<S,V> &ret,Grid_simd<S,V> b,int nrot)
{
  nrot = nrot % Grid_simd<S,V>::Nsimd();
  ret.v = Optimization::Rotate::rotate(b.v,nrot);
}
template <class S, class V, IfComplex<S> =0>
accelerator_inline void rotate(Grid_simd<S,V> &ret,Grid_simd<S,V> b,int nrot)
{
  nrot = nrot % Grid_simd<S,V>::Nsimd();
  ret.v = Optimization::Rotate::rotate(b.v,2*nrot);
}

template <class S, class V>
accelerator_inline void vbroadcast(Grid_simd<S,V> &ret,const Grid_simd<S,V> &src,int lane){
  S* typepun =(S*) &src;
  vsplat(ret,typepun[lane]);
}
template <class S, class V, IfComplex<S> =0>
accelerator_inline void rbroadcast(Grid_simd<S,V> &ret,const Grid_simd<S,V> &src,int lane){
  S* typepun =(S*) &src;
  ret.v = unary<V>(real(typepun[lane]), VsplatSIMD());
}



///////////////////////
// Splat
///////////////////////

// this is only for the complex version
template <class S, class V, IfComplex<S> = 0, class ABtype>
accelerator_inline void vsplat(Grid_simd<S, V> &ret, ABtype a, ABtype b) {
  ret.v = binary<V>(a, b, VsplatSIMD());
}

// overload if complex
template <class S, class V>
accelerator_inline void vsplat(Grid_simd<S, V> &ret, EnableIf<is_complex<S>, S> c) {
  vsplat(ret, real(c), imag(c));
}
template <class S, class V>
accelerator_inline void rsplat(Grid_simd<S, V> &ret, EnableIf<is_complex<S>, S> c) {
  vsplat(ret, real(c), real(c));
}

// if real fill with a, if complex fill with a in the real part (first function
// above)
template <class S, class V>
accelerator_inline void vsplat(Grid_simd<S, V> &ret, NotEnableIf<is_complex<S>, S> a) {
  ret.v = unary<V>(a, VsplatSIMD());
}
//////////////////////////


///////////////////////////////////////////////
// Initialise to 1,0,i for the correct types
///////////////////////////////////////////////
// For complex types
template <class S, class V, IfComplex<S> = 0>
accelerator_inline void vone(Grid_simd<S, V> &ret) {
  vsplat(ret, S(1.0, 0.0));
}
template <class S, class V, IfComplex<S> = 0>
accelerator_inline void vzero(Grid_simd<S, V> &ret) {
  vsplat(ret, S(0.0, 0.0));
}  // use xor?
template <class S, class V, IfComplex<S> = 0>
accelerator_inline void vcomplex_i(Grid_simd<S, V> &ret) {
  vsplat(ret, S(0.0, 1.0));
}

template <class S, class V, IfComplex<S> = 0>
accelerator_inline void visign(Grid_simd<S, V> &ret) {
  vsplat(ret, S(1.0, -1.0));
}
template <class S, class V, IfComplex<S> = 0>
accelerator_inline void vrsign(Grid_simd<S, V> &ret) {
  vsplat(ret, S(-1.0, 1.0));
}

// if not complex overload here
template <class S, class V, IfReal<S> = 0>
accelerator_inline void vone(Grid_simd<S, V> &ret) {
  vsplat(ret, S(1.0));
}
template <class S, class V, IfReal<S> = 0>
accelerator_inline void vzero(Grid_simd<S, V> &ret) {
  vsplat(ret, S(0.0));
}

// For integral types
template <class S, class V, IfInteger<S> = 0>
accelerator_inline void vone(Grid_simd<S, V> &ret) {
  vsplat(ret, 1);
}
template <class S, class V, IfInteger<S> = 0>
accelerator_inline void vzero(Grid_simd<S, V> &ret) {
  vsplat(ret, 0);
}
template <class S, class V, IfInteger<S> = 0>
accelerator_inline void vtrue(Grid_simd<S, V> &ret) {
  vsplat(ret, 0xFFFFFFFF);
}
template <class S, class V, IfInteger<S> = 0>
accelerator_inline void vfalse(Grid_simd<S, V> &ret) {
  vsplat(ret, 0);
}
template <class S, class V>
accelerator_inline void zeroit(Grid_simd<S, V> &z) {
  vzero(z);
}

///////////////////////
// Vstream
///////////////////////
template <class S, class V, IfReal<S> = 0>
accelerator_inline void vstream(Grid_simd<S, V> &out, const Grid_simd<S, V> &in) {
  binary<void>((S *)&out.v, in.v, VstreamSIMD());
}
template <class S, class V, IfComplex<S> = 0>
accelerator_inline void vstream(Grid_simd<S, V> &out, const Grid_simd<S, V> &in) {
  typedef typename S::value_type T;
  binary<void>((T *)&out.v, in.v, VstreamSIMD());
}
template <class S, class V, IfInteger<S> = 0>
accelerator_inline void vstream(Grid_simd<S, V> &out, const Grid_simd<S, V> &in) {
  out = in;
}

////////////////////////////////////
// Arithmetic operator overloads +,-,*
////////////////////////////////////
template <class S, class V>
accelerator_inline Grid_simd<S, V> operator+(Grid_simd<S, V> a, Grid_simd<S, V> b) {
  Grid_simd<S, V> ret;
  ret.v = binary<V>(a.v, b.v, SumSIMD());
  return ret;
};

template <class S, class V>
accelerator_inline Grid_simd<S, V> operator-(Grid_simd<S, V> a, Grid_simd<S, V> b) {
  Grid_simd<S, V> ret;
  ret.v = binary<V>(a.v, b.v, SubSIMD());
  return ret;
};

// Distinguish between complex types and others
template <class S, class V, IfComplex<S> = 0>
accelerator_inline Grid_simd<S, V> real_mult(Grid_simd<S, V> a, Grid_simd<S, V> b) {
  Grid_simd<S, V> ret;
  ret.v = binary<V>(a.v, b.v, MultRealPartSIMD());
  return ret;
};
template <class S, class V, IfComplex<S> = 0>
accelerator_inline Grid_simd<S, V> real_madd(Grid_simd<S, V> a, Grid_simd<S, V> b, Grid_simd<S,V> c) {
  Grid_simd<S, V> ret;
  ret.v = trinary<V>(a.v, b.v, c.v, MaddRealPartSIMD());
  return ret;
};


// Distinguish between complex types and others
template <class S, class V, IfComplex<S> = 0>
accelerator_inline Grid_simd<S, V> operator*(Grid_simd<S, V> a, Grid_simd<S, V> b) {
  Grid_simd<S, V> ret;
  ret.v = binary<V>(a.v, b.v, MultComplexSIMD());
  return ret;
};

// Real/Integer types
template <class S, class V, IfNotComplex<S> = 0>
accelerator_inline Grid_simd<S, V> operator*(Grid_simd<S, V> a, Grid_simd<S, V> b) {
  Grid_simd<S, V> ret;
  ret.v = binary<V>(a.v, b.v, MultSIMD());
  return ret;
};

// ---------------- A64FX MAC -------------------
// Distinguish between complex types and others
#if defined(A64FX) || defined(A64FXFIXEDSIZE)
template <class S, class V, IfComplex<S> = 0>
accelerator_inline Grid_simd<S, V> fxmac(Grid_simd<S, V> a, Grid_simd<S, V> b, Grid_simd<S, V> c) {
  Grid_simd<S, V> ret;
  ret.v = trinary<V>(a.v, b.v, c.v, MultAddComplexSIMD());
  return ret;
};

// Real/Integer types
template <class S, class V, IfNotComplex<S> = 0>
accelerator_inline Grid_simd<S, V> fxmac(Grid_simd<S, V> a, Grid_simd<S, V> b, Grid_simd<S, V> c) {
  Grid_simd<S, V> ret;
  ret.v = trinary<V>(a.v, b.v, c.v, MultSIMD());
  return ret;
};
#endif
// ----------------------------------------------


///////////////////////
// Conjugate
///////////////////////
template <class S, class V, IfComplex<S> = 0>
accelerator_inline Grid_simd<S, V> conjugate(const Grid_simd<S, V> &in) {
  Grid_simd<S, V> ret;
  ret.v = unary<V>(in.v, ConjSIMD());
  return ret;
}
template <class S, class V, IfNotComplex<S> = 0>
accelerator_inline Grid_simd<S, V> conjugate(const Grid_simd<S, V> &in) {
  return in;  // for real objects
}
// Suppress adj for integer types... // odd; why conjugate above but not adj??
template <class S, class V, IfNotInteger<S> = 0>
accelerator_inline Grid_simd<S, V> adj(const Grid_simd<S, V> &in) {
  return conjugate(in);
}

///////////////////////
// timesMinusI
///////////////////////
template <class S, class V, IfComplex<S> = 0>
accelerator_inline void timesMinusI(Grid_simd<S, V> &ret, const Grid_simd<S, V> &in) {
  ret.v = unary<V>(in.v, TimesMinusISIMD());
}
template <class S, class V, IfComplex<S> = 0>
accelerator_inline Grid_simd<S, V> timesMinusI(const Grid_simd<S, V> &in) {
  Grid_simd<S, V> ret;
  ret.v=unary<V>(in.v, TimesMinusISIMD());
  return ret;
}
template <class S, class V, IfNotComplex<S> = 0>
accelerator_inline Grid_simd<S, V> timesMinusI(const Grid_simd<S, V> &in) {
  return in;
}
///////////////////////
// timesI
///////////////////////
template <class S, class V, IfComplex<S> = 0>
accelerator_inline void timesI(Grid_simd<S, V> &ret, const Grid_simd<S, V> &in) {
  ret.v = unary<V>(in.v, TimesISIMD());
}
template <class S, class V, IfComplex<S> = 0>
accelerator_inline Grid_simd<S, V> timesI(const Grid_simd<S, V> &in) {
  Grid_simd<S, V> ret;
  ret.v= unary<V>(in.v, TimesISIMD());
  return ret;
}
template <class S, class V, IfNotComplex<S> = 0>
accelerator_inline Grid_simd<S, V> timesI(const Grid_simd<S, V> &in) {
  return in;
}


// Distinguish between complex types and others
template <class S, class V, IfComplex<S> = 0>
accelerator_inline Grid_simd<S, V> operator/(Grid_simd<S, V> a, Grid_simd<S, V> b) {
  typedef Grid_simd<S, V> simd;

  simd ret;
  simd den;

  ret = a * conjugate(b) ;
  den = b * conjugate(b) ;

  // duplicates real part
  auto real_den  = toReal(den);
  simd zden;
  memcpy((void *)&zden.v,(void *)&real_den.v,sizeof(zden));
  ret.v=binary<V>(ret.v, zden.v, DivSIMD());
  return ret;
};

// Real/Integer types
template <class S, class V, IfNotComplex<S> = 0>
accelerator_inline Grid_simd<S, V> operator/(Grid_simd<S, V> a, Grid_simd<S, V> b) {
  Grid_simd<S, V> ret;
  ret.v = binary<V>(a.v, b.v, DivSIMD());
  return ret;
};


/////////////////////
// Inner, outer
/////////////////////
template <class S, class V>
accelerator_inline Grid_simd<S, V> innerProduct(const Grid_simd<S, V> &l,const Grid_simd<S, V> &r) {
  return conjugate(l) * r;
}
template <class S, class V>
accelerator_inline Grid_simd<S, V> outerProduct(const Grid_simd<S, V> &l,const Grid_simd<S, V> &r) {
  return l * conjugate(r);
}

template <class S, class V>
accelerator_inline Grid_simd<S, V> trace(const Grid_simd<S, V> &arg) {
  return arg;
}

////////////////////////////////////////////////////////////
// copy/splat complex real parts into real;
// insert real into complex and zero imag;
////////////////////////////////////////////////////////////


template <class T> struct toRealMapper {};
template<> struct toRealMapper<vComplexF> {  typedef vRealF Realified; };
template<> struct toRealMapper<vComplexD> {  typedef vRealD Realified; };
// real = toReal( complex )
template <class Csimd>  // must be a real arg
accelerator_inline typename toRealMapper<Csimd>::Realified toReal(const Csimd &in) {
  typedef typename toRealMapper<Csimd>::Realified Rsimd;
  Rsimd ret;
  int j=0;
  for (int i = 0; i < Rsimd::Nsimd(); i += 2) {
    auto s = real(in.getlane(j++));
    ret.putlane(s,i);
    ret.putlane(s,i+1);
  }
  return ret;
}

template <class T> struct toComplexMapper {};
template<> struct toComplexMapper<vRealF> {  typedef vComplexF Complexified; };
template<> struct toComplexMapper<vRealD> {  typedef vComplexD Complexified; };

// complex = toComplex( real )
template <class Rsimd>  // must be a real arg
accelerator_inline typename toComplexMapper<Rsimd>::Complexified toComplex(const Rsimd &in) {

  typedef typename toComplexMapper<Rsimd>::Complexified   Csimd;
  typedef typename Csimd::scalar_type scalar_type;
  int j=0;
  Csimd ret;
  for (int i = 0; i < Rsimd::Nsimd(); i += 2) {
    auto rr = in.getlane(i);
    auto ri = in.getlane(i+1);
    assert(rr==ri);
    // trap any cases where real was not duplicated
    // indicating the SIMD grids of real and imag assignment did not correctly
    // match
    scalar_type s(rr,0.0);
    ret.putlane(s,j++);
  }
  return ret;
}


accelerator_inline void precisionChange(vRealF    *out,const vRealD    *in,int nvec)
{
  assert((nvec&0x1)==0);
  for(int m=0;m*2<nvec;m++){
    int n=m*2;
    out[m].v=Optimization::PrecisionChange::DtoS(in[n].v,in[n+1].v);
  }
}
accelerator_inline void precisionChange(vRealH    *out,const vRealD    *in,int nvec)
{
  assert((nvec&0x3)==0);
  for(int m=0;m*4<nvec;m++){
    int n=m*4;
    out[m].v=Optimization::PrecisionChange::DtoH(in[n].v,in[n+1].v,in[n+2].v,in[n+3].v);
  }
}
accelerator_inline void precisionChange(vRealH    *out,const vRealF    *in,int nvec)
{
  assert((nvec&0x1)==0);
  for(int m=0;m*2<nvec;m++){
    int n=m*2;
    out[m].v=Optimization::PrecisionChange::StoH(in[n].v,in[n+1].v);
  }
}
accelerator_inline void precisionChange(vRealD    *out,const vRealF    *in,int nvec)
{
  assert((nvec&0x1)==0);
  for(int m=0;m*2<nvec;m++){
    int n=m*2;
    Optimization::PrecisionChange::StoD(in[m].v,out[n].v,out[n+1].v);
    // Bug in gcc 10.0.1 and gcc 10.1 using fixed-size SVE ACLE data types  CAS-159553-Y1K4C6
    // function call results in compile-time error:
    // In function ‘void Grid::precisionChange(Grid::vRealD*, Grid::vRealF*, int)’:
    // .../Grid_vector_types.h:961:56: error:
    // cannot bind non-const lvalue reference of type ‘vecd&’ {aka ‘svfloat64_t&’}
    // to an rvalue of type ‘vecd’ {aka ‘svfloat64_t’}
    // 961 |     Optimization::PrecisionChange::StoD(in[m].v,out[n].v,out[n+1].v);
    //  |                                                 ~~~~~~~^
  }
}
accelerator_inline void precisionChange(vRealD    *out,const vRealH    *in,int nvec)
{
  assert((nvec&0x3)==0);
  for(int m=0;m*4<nvec;m++){
    int n=m*4;
    Optimization::PrecisionChange::HtoD(in[m].v,out[n].v,out[n+1].v,out[n+2].v,out[n+3].v);
  }
}
accelerator_inline void precisionChange(vRealF    *out,const vRealH    *in,int nvec)
{
  assert((nvec&0x1)==0);
  for(int m=0;m*2<nvec;m++){
    int n=m*2;
    Optimization::PrecisionChange::HtoS(in[m].v,out[n].v,out[n+1].v);
  }
}
accelerator_inline void precisionChange(vComplexF *out,const vComplexD *in,int nvec){ precisionChange((vRealF *)out,(vRealD *)in,nvec);}
accelerator_inline void precisionChange(vComplexH *out,const vComplexD *in,int nvec){ precisionChange((vRealH *)out,(vRealD *)in,nvec);}
accelerator_inline void precisionChange(vComplexH *out,const vComplexF *in,int nvec){ precisionChange((vRealH *)out,(vRealF *)in,nvec);}
accelerator_inline void precisionChange(vComplexD *out,const vComplexF *in,int nvec){ precisionChange((vRealD *)out,(vRealF *)in,nvec);}
accelerator_inline void precisionChange(vComplexD *out,const vComplexH *in,int nvec){ precisionChange((vRealD *)out,(vRealH *)in,nvec);}
accelerator_inline void precisionChange(vComplexF *out,const vComplexH *in,int nvec){ precisionChange((vRealF *)out,(vRealH *)in,nvec);}

// Check our vector types are of an appropriate size.

#if defined QPX
static_assert(2*sizeof(SIMD_Ftype) == sizeof(SIMD_Dtype), "SIMD vector lengths incorrect");
static_assert(2*sizeof(SIMD_Ftype) == sizeof(SIMD_Itype), "SIMD vector lengths incorrect");
#else
#ifndef GENERIC_SCALAR
static_assert(sizeof(SIMD_Ftype) == sizeof(SIMD_Dtype), "SIMD vector lengths incorrect");
static_assert(sizeof(SIMD_Ftype) == sizeof(SIMD_Itype), "SIMD vector lengths incorrect");
#endif
#endif

// Fixme need coalesced read gpermute
template<class vobj> void gpermute(vobj & inout,int perm){
  vobj tmp=inout;
  if (perm & 0x1 ) { permute(inout,tmp,0); tmp=inout;}
  if (perm & 0x2 ) { permute(inout,tmp,1); tmp=inout;}
  if (perm & 0x4 ) { permute(inout,tmp,2); tmp=inout;}
  if (perm & 0x8 ) { permute(inout,tmp,3); tmp=inout;}
}

NAMESPACE_END(Grid);

#ifdef GRID_SYCL
template<> struct sycl::is_device_copyable<Grid::vComplexF> : public std::true_type {};
template<> struct sycl::is_device_copyable<Grid::vComplexD> : public std::true_type {};
template<> struct sycl::is_device_copyable<Grid::vRealF   > : public std::true_type {};
template<> struct sycl::is_device_copyable<Grid::vRealD   > : public std::true_type {};
template<> struct sycl::is_device_copyable<Grid::vInteger > : public std::true_type {};
#endif


#endif
