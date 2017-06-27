/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/simd/Grid_vector_type.h

Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
Author: Guido Cossu <cossu@iroiro-pc.kek.jp>
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

See the full license in the file "LICENSE" in the top level distribution
directory
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

#ifdef GEN
#include "Grid_generic.h"
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
#ifdef NEONv8
#include "Grid_neon.h"
#endif
#if defined QPX
#include "Grid_qpx.h"
#endif

#include "l1p.h"

namespace Grid {

//////////////////////////////////////
// To take the floating point type of real/complex type
//////////////////////////////////////
template <typename T>
struct RealPart {
  typedef T type;
};
template <typename T>
struct RealPart<std::complex<T> > {
  typedef T type;
};

#include <type_traits>

//////////////////////////////////////
// demote a vector to real type
//////////////////////////////////////
// type alias used to simplify the syntax of std::enable_if
template <typename T> using Invoke = typename T::type;
template <typename Condition, typename ReturnType> using EnableIf    = Invoke<std::enable_if<Condition::value, ReturnType> >;
template <typename Condition, typename ReturnType> using NotEnableIf = Invoke<std::enable_if<!Condition::value, ReturnType> >;

////////////////////////////////////////////////////////
// Check for complexity with type traits
template <typename T> struct is_complex : public std::false_type {};
template <> struct is_complex<std::complex<double> > : public std::true_type {};
template <> struct is_complex<std::complex<float> > : public std::true_type {};

template <typename T>              using IfReal    = Invoke<std::enable_if<std::is_floating_point<T>::value, int> >;
template <typename T>              using IfComplex = Invoke<std::enable_if<is_complex<T>::value, int> >;
template <typename T>              using IfInteger = Invoke<std::enable_if<std::is_integral<T>::value, int> >;
template <typename T1,typename T2> using IfSame    = Invoke<std::enable_if<std::is_same<T1,T2>::value, int> >;

template <typename T>              using IfNotReal    = Invoke<std::enable_if<!std::is_floating_point<T>::value, int> >;
template <typename T>              using IfNotComplex = Invoke<std::enable_if<!is_complex<T>::value, int> >;
template <typename T>              using IfNotInteger = Invoke<std::enable_if<!std::is_integral<T>::value, int> >;
template <typename T1,typename T2> using IfNotSame    = Invoke<std::enable_if<!std::is_same<T1,T2>::value, int> >;

////////////////////////////////////////////////////////
// Define the operation templates functors
// general forms to allow for vsplat syntax
// need explicit declaration of types when used since
// clang cannot automatically determine the output type sometimes
template <class Out, class Input1, class Input2, class Input3, class Operation>
Out trinary(Input1 src_1, Input2 src_2, Input3 src_3, Operation op) {
  return op(src_1, src_2, src_3);
}

template <class Out, class Input1, class Input2, class Operation>
Out binary(Input1 src_1, Input2 src_2, Operation op) {
  return op(src_1, src_2);
}

template <class Out, class Input, class Operation>
Out unary(Input src, Operation op) {
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

  typedef union conv_t_union {
    Vector_type v;
    Scalar_type s[sizeof(Vector_type) / sizeof(Scalar_type)];
    conv_t_union(){};
  } conv_t;

  Vector_type v;

  static inline constexpr int Nsimd(void) {
    return sizeof(Vector_type) / sizeof(Scalar_type);
  }

  Grid_simd &operator=(const Grid_simd &&rhs) {
    v = rhs.v;
    return *this;
  };
  Grid_simd &operator=(const Grid_simd &rhs) {
    v = rhs.v;
    return *this;
  };  // faster than not declaring it and leaving to the compiler
  Grid_simd() = default;
  Grid_simd(const Grid_simd &rhs) : v(rhs.v){};  // compiles in movaps
  Grid_simd(const Grid_simd &&rhs) : v(rhs.v){};

  /////////////////////////////
  // Constructors
  /////////////////////////////
  Grid_simd &operator=(Zero &z) {
    vzero(*this);
    return (*this);
  }

  // Enable if complex type
  template <typename S = Scalar_type>
  Grid_simd(const typename std::enable_if<is_complex<S>::value, S>::type a) {
    vsplat(*this, a);
  };

  Grid_simd(const Real a) { vsplat(*this, Scalar_type(a)); };

  ///////////////////////////////////////////////
  // mac, mult, sub, add, adj
  ///////////////////////////////////////////////

  // FIXME -- alias this to an inline MAC struct.
  friend inline void mac(Grid_simd *__restrict__ y,
                         const Grid_simd *__restrict__ a,
                         const Grid_simd *__restrict__ x) {
    *y = (*a) * (*x) + (*y);
  };

  friend inline void mult(Grid_simd *__restrict__ y,
                          const Grid_simd *__restrict__ l,
                          const Grid_simd *__restrict__ r) {
    *y = (*l) * (*r);
  }

  friend inline void sub(Grid_simd *__restrict__ y,
                         const Grid_simd *__restrict__ l,
                         const Grid_simd *__restrict__ r) {
    *y = (*l) - (*r);
  }
  friend inline void add(Grid_simd *__restrict__ y,
                         const Grid_simd *__restrict__ l,
                         const Grid_simd *__restrict__ r) {
    *y = (*l) + (*r);
  }
  friend inline void mac(Grid_simd *__restrict__ y,
                         const Scalar_type *__restrict__ a,
                         const Grid_simd *__restrict__ x) {
    *y = (*a) * (*x) + (*y);
  };
  friend inline void mult(Grid_simd *__restrict__ y,
                          const Scalar_type *__restrict__ l,
                          const Grid_simd *__restrict__ r) {
    *y = (*l) * (*r);
  }
  friend inline void sub(Grid_simd *__restrict__ y,
                         const Scalar_type *__restrict__ l,
                         const Grid_simd *__restrict__ r) {
    *y = (*l) - (*r);
  }
  friend inline void add(Grid_simd *__restrict__ y,
                         const Scalar_type *__restrict__ l,
                         const Grid_simd *__restrict__ r) {
    *y = (*l) + (*r);
  }

  friend inline void mac(Grid_simd *__restrict__ y,
                         const Grid_simd *__restrict__ a,
                         const Scalar_type *__restrict__ x) {
    *y = (*a) * (*x) + (*y);
  };
  friend inline void mult(Grid_simd *__restrict__ y,
                          const Grid_simd *__restrict__ l,
                          const Scalar_type *__restrict__ r) {
    *y = (*l) * (*r);
  }
  friend inline void sub(Grid_simd *__restrict__ y,
                         const Grid_simd *__restrict__ l,
                         const Scalar_type *__restrict__ r) {
    *y = (*l) - (*r);
  }
  friend inline void add(Grid_simd *__restrict__ y,
                         const Grid_simd *__restrict__ l,
                         const Scalar_type *__restrict__ r) {
    *y = (*l) + (*r);
  }

  ////////////////////////////////////////////////////////////////////////
  // FIXME:  gonna remove these load/store, get, set, prefetch
  ////////////////////////////////////////////////////////////////////////
  friend inline void vset(Grid_simd &ret, Scalar_type *a) {
    ret.v = unary<Vector_type>(a, VsetSIMD());
  }

  ///////////////////////
  // Vstore
  ///////////////////////
  friend inline void vstore(const Grid_simd &ret, Scalar_type *a) {
    binary<void>(ret.v, (Real *)a, VstoreSIMD());
  }

  ///////////////////////
  // Vprefetch
  ///////////////////////
  friend inline void vprefetch(const Grid_simd &v) {
    prefetch_HINT_T0((const char *)&v.v);
  }

  ///////////////////////
  // Reduce
  ///////////////////////
  friend inline Scalar_type Reduce(const Grid_simd &in) {
    return unary<Scalar_type>(in.v, ReduceSIMD<Scalar_type, Vector_type>());
  }

  ////////////////////////////
  // operator scalar * simd
  ////////////////////////////
  friend inline Grid_simd operator*(const Scalar_type &a, Grid_simd b) {
    Grid_simd va;
    vsplat(va, a);
    return va * b;
  }
  friend inline Grid_simd operator*(Grid_simd b, const Scalar_type &a) {
    return a * b;
  }

  //////////////////////////////////
  // Divides
  //////////////////////////////////
  friend inline Grid_simd operator/(const Scalar_type &a, Grid_simd b) {
    Grid_simd va;
    vsplat(va, a);
    return va / b;
  }
  friend inline Grid_simd operator/(Grid_simd b, const Scalar_type &a) {
    Grid_simd va;
    vsplat(va, a);
    return b / a;
  }

  ///////////////////////
  // Unary negation
  ///////////////////////
  friend inline Grid_simd operator-(const Grid_simd &r) {
    Grid_simd ret;
    vzero(ret);
    ret = ret - r;
    return ret;
  }
  // *=,+=,-= operators
  inline Grid_simd &operator*=(const Grid_simd &r) {
    *this = (*this) * r;
    return *this;
    // return (*this)*r; ?
  }
  inline Grid_simd &operator+=(const Grid_simd &r) {
    *this = *this + r;
    return *this;
  }
  inline Grid_simd &operator-=(const Grid_simd &r) {
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
  friend inline Grid_simd SimdApply(const functor &func, const Grid_simd &v) {
    Grid_simd ret;
    Grid_simd::conv_t conv;
    Grid_simd::scalar_type s;
    
    conv.v = v.v;
    for (int i = 0; i < Nsimd(); i++) {
      s = conv.s[i];
      conv.s[i] = func(s);
    }
    ret.v = conv.v;
    return ret;
  }
  template <class functor>
  friend inline Grid_simd SimdApplyBinop(const functor &func,
                                         const Grid_simd &x,
                                         const Grid_simd &y) {
    Grid_simd ret;
    Grid_simd::conv_t cx;
    Grid_simd::conv_t cy;
    Grid_simd::scalar_type sx,sy;

    cx.v = x.v;
    cy.v = y.v;
    for (int i = 0; i < Nsimd(); i++) {
      sx = cx.s[i];
      sy = cy.s[i];
      cx.s[i] = func(sx,sy);
    }
    ret.v = cx.v;
    return ret;
  }
  ///////////////////////
  // Exchange 
  // Al Ah , Bl Bh -> Al Bl Ah,Bh
  ///////////////////////
  friend inline void exchange(Grid_simd &out1,Grid_simd &out2,Grid_simd in1,Grid_simd in2,int n)
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

  ////////////////////////////////////////////////////////////////////
  // General permute; assumes vector length is same across
  // all subtypes; may not be a good assumption, but could
  // add the vector width as a template param for BG/Q for example
  ////////////////////////////////////////////////////////////////////
  friend inline void permute0(Grid_simd &y, Grid_simd b) {
    y.v = Optimization::Permute::Permute0(b.v);
  }
  friend inline void permute1(Grid_simd &y, Grid_simd b) {
    y.v = Optimization::Permute::Permute1(b.v);
  }
  friend inline void permute2(Grid_simd &y, Grid_simd b) {
    y.v = Optimization::Permute::Permute2(b.v);
  }
  friend inline void permute3(Grid_simd &y, Grid_simd b) {
    y.v = Optimization::Permute::Permute3(b.v);
  }
  friend inline void permute(Grid_simd &y, Grid_simd b, int perm) {
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
  inline Scalar_type getlane(int lane) {
    return ((Scalar_type*)&v)[lane];
  }

  inline void putlane(const Scalar_type &S, int lane){
    ((Scalar_type*)&v)[lane] = S;
  }


  
};  // end of Grid_simd class definition

inline void permute(ComplexD &y,ComplexD b, int perm) {  y=b; }
inline void permute(ComplexF &y,ComplexF b, int perm) {  y=b; }
inline void permute(RealD &y,RealD b, int perm) {  y=b; }
inline void permute(RealF &y,RealF b, int perm) {  y=b; }

////////////////////////////////////////////////////////////////////
// General rotate
////////////////////////////////////////////////////////////////////
template <class S, class V, IfNotComplex<S> = 0>
inline Grid_simd<S, V> rotate(Grid_simd<S, V> b, int nrot) {
  nrot = nrot % Grid_simd<S, V>::Nsimd();
  Grid_simd<S, V> ret;
  ret.v = Optimization::Rotate::rotate(b.v, nrot);
  return ret;
}
template <class S, class V, IfComplex<S> = 0>
inline Grid_simd<S, V> rotate(Grid_simd<S, V> b, int nrot) {
  nrot = nrot % Grid_simd<S, V>::Nsimd();
  Grid_simd<S, V> ret;
  ret.v = Optimization::Rotate::rotate(b.v, 2 * nrot);
  return ret;
}
template <class S, class V, IfNotComplex<S> =0> 
inline void rotate( Grid_simd<S,V> &ret,Grid_simd<S,V> b,int nrot)
{
  nrot = nrot % Grid_simd<S,V>::Nsimd();
  ret.v = Optimization::Rotate::rotate(b.v,nrot);
}
template <class S, class V, IfComplex<S> =0> 
inline void rotate(Grid_simd<S,V> &ret,Grid_simd<S,V> b,int nrot)
{
  nrot = nrot % Grid_simd<S,V>::Nsimd();
  ret.v = Optimization::Rotate::rotate(b.v,2*nrot);
}

template <class S, class V> 
inline void vbroadcast(Grid_simd<S,V> &ret,const Grid_simd<S,V> &src,int lane){
  S* typepun =(S*) &src;
  vsplat(ret,typepun[lane]);
}    
template <class S, class V, IfComplex<S> =0> 
inline void rbroadcast(Grid_simd<S,V> &ret,const Grid_simd<S,V> &src,int lane){
  S* typepun =(S*) &src;
  ret.v = unary<V>(real(typepun[lane]), VsplatSIMD());
}    



///////////////////////
// Splat
///////////////////////

// this is only for the complex version
template <class S, class V, IfComplex<S> = 0, class ABtype>
inline void vsplat(Grid_simd<S, V> &ret, ABtype a, ABtype b) {
  ret.v = binary<V>(a, b, VsplatSIMD());
}

// overload if complex
template <class S, class V>
inline void vsplat(Grid_simd<S, V> &ret, EnableIf<is_complex<S>, S> c) {
  vsplat(ret, real(c), imag(c));
}
template <class S, class V>
inline void rsplat(Grid_simd<S, V> &ret, EnableIf<is_complex<S>, S> c) {
  vsplat(ret, real(c), real(c));
}

// if real fill with a, if complex fill with a in the real part (first function
// above)
template <class S, class V>
inline void vsplat(Grid_simd<S, V> &ret, NotEnableIf<is_complex<S>, S> a) {
  ret.v = unary<V>(a, VsplatSIMD());
}
//////////////////////////

///////////////////////////////////////////////
// Initialise to 1,0,i for the correct types
///////////////////////////////////////////////
// For complex types
template <class S, class V, IfComplex<S> = 0>
inline void vone(Grid_simd<S, V> &ret) {
  vsplat(ret, S(1.0, 0.0));
}
template <class S, class V, IfComplex<S> = 0>
inline void vzero(Grid_simd<S, V> &ret) {
  vsplat(ret, S(0.0, 0.0));
}  // use xor?
template <class S, class V, IfComplex<S> = 0>
inline void vcomplex_i(Grid_simd<S, V> &ret) {
  vsplat(ret, S(0.0, 1.0));
}

template <class S, class V, IfComplex<S> = 0>
inline void visign(Grid_simd<S, V> &ret) {
  vsplat(ret, S(1.0, -1.0));
}
template <class S, class V, IfComplex<S> = 0>
inline void vrsign(Grid_simd<S, V> &ret) {
  vsplat(ret, S(-1.0, 1.0));
}

// if not complex overload here
template <class S, class V, IfReal<S> = 0>
inline void vone(Grid_simd<S, V> &ret) {
  vsplat(ret, S(1.0));
}
template <class S, class V, IfReal<S> = 0>
inline void vzero(Grid_simd<S, V> &ret) {
  vsplat(ret, S(0.0));
}

// For integral types
template <class S, class V, IfInteger<S> = 0>
inline void vone(Grid_simd<S, V> &ret) {
  vsplat(ret, 1);
}
template <class S, class V, IfInteger<S> = 0>
inline void vzero(Grid_simd<S, V> &ret) {
  vsplat(ret, 0);
}
template <class S, class V, IfInteger<S> = 0>
inline void vtrue(Grid_simd<S, V> &ret) {
  vsplat(ret, 0xFFFFFFFF);
}
template <class S, class V, IfInteger<S> = 0>
inline void vfalse(Grid_simd<S, V> &ret) {
  vsplat(ret, 0);
}
template <class S, class V>
inline void zeroit(Grid_simd<S, V> &z) {
  vzero(z);
}

///////////////////////
// Vstream
///////////////////////
template <class S, class V, IfReal<S> = 0>
inline void vstream(Grid_simd<S, V> &out, const Grid_simd<S, V> &in) {
  binary<void>((S *)&out.v, in.v, VstreamSIMD());
}
template <class S, class V, IfComplex<S> = 0>
inline void vstream(Grid_simd<S, V> &out, const Grid_simd<S, V> &in) {
  typedef typename S::value_type T;
  binary<void>((T *)&out.v, in.v, VstreamSIMD());
}
template <class S, class V, IfInteger<S> = 0>
inline void vstream(Grid_simd<S, V> &out, const Grid_simd<S, V> &in) {
  out = in;
}

////////////////////////////////////
// Arithmetic operator overloads +,-,*
////////////////////////////////////
template <class S, class V>
inline Grid_simd<S, V> operator+(Grid_simd<S, V> a, Grid_simd<S, V> b) {
  Grid_simd<S, V> ret;
  ret.v = binary<V>(a.v, b.v, SumSIMD());
  return ret;
};

template <class S, class V>
inline Grid_simd<S, V> operator-(Grid_simd<S, V> a, Grid_simd<S, V> b) {
  Grid_simd<S, V> ret;
  ret.v = binary<V>(a.v, b.v, SubSIMD());
  return ret;
};

// Distinguish between complex types and others
template <class S, class V, IfComplex<S> = 0>
inline Grid_simd<S, V> real_mult(Grid_simd<S, V> a, Grid_simd<S, V> b) {
  Grid_simd<S, V> ret;
  ret.v = binary<V>(a.v, b.v, MultRealPartSIMD());
  return ret;
};
template <class S, class V, IfComplex<S> = 0>
inline Grid_simd<S, V> real_madd(Grid_simd<S, V> a, Grid_simd<S, V> b, Grid_simd<S,V> c) {
  Grid_simd<S, V> ret;
  ret.v = trinary<V>(a.v, b.v, c.v, MaddRealPartSIMD());
  return ret;
};


// Distinguish between complex types and others
template <class S, class V, IfComplex<S> = 0>
inline Grid_simd<S, V> operator*(Grid_simd<S, V> a, Grid_simd<S, V> b) {
  Grid_simd<S, V> ret;
  ret.v = binary<V>(a.v, b.v, MultComplexSIMD());
  return ret;
};

// Real/Integer types
template <class S, class V, IfNotComplex<S> = 0>
inline Grid_simd<S, V> operator*(Grid_simd<S, V> a, Grid_simd<S, V> b) {
  Grid_simd<S, V> ret;
  ret.v = binary<V>(a.v, b.v, MultSIMD());
  return ret;
};

// Distinguish between complex types and others
template <class S, class V, IfComplex<S> = 0>
inline Grid_simd<S, V> operator/(Grid_simd<S, V> a, Grid_simd<S, V> b) {
  typedef Grid_simd<S, V> simd;

  simd ret;
  simd den;
  typename simd::conv_t conv;

  ret = a * conjugate(b) ;
  den = b * conjugate(b) ;

  
  auto real_den = toReal(den);

  ret.v=binary<V>(ret.v, real_den.v, DivSIMD());

  return ret;
};

// Real/Integer types
template <class S, class V, IfNotComplex<S> = 0>
inline Grid_simd<S, V> operator/(Grid_simd<S, V> a, Grid_simd<S, V> b) {
  Grid_simd<S, V> ret;
  ret.v = binary<V>(a.v, b.v, DivSIMD());
  return ret;
};

///////////////////////
// Conjugate
///////////////////////
template <class S, class V, IfComplex<S> = 0>
inline Grid_simd<S, V> conjugate(const Grid_simd<S, V> &in) {
  Grid_simd<S, V> ret;
  ret.v = unary<V>(in.v, ConjSIMD());
  return ret;
}
template <class S, class V, IfNotComplex<S> = 0>
inline Grid_simd<S, V> conjugate(const Grid_simd<S, V> &in) {
  return in;  // for real objects
}
// Suppress adj for integer types... // odd; why conjugate above but not adj??
template <class S, class V, IfNotInteger<S> = 0>
inline Grid_simd<S, V> adj(const Grid_simd<S, V> &in) {
  return conjugate(in);
}

///////////////////////
// timesMinusI
///////////////////////
template <class S, class V, IfComplex<S> = 0>
inline void timesMinusI(Grid_simd<S, V> &ret, const Grid_simd<S, V> &in) {
  ret.v = binary<V>(in.v, ret.v, TimesMinusISIMD());
}
template <class S, class V, IfComplex<S> = 0>
inline Grid_simd<S, V> timesMinusI(const Grid_simd<S, V> &in) {
  Grid_simd<S, V> ret;
  timesMinusI(ret, in);
  return ret;
}
template <class S, class V, IfNotComplex<S> = 0>
inline Grid_simd<S, V> timesMinusI(const Grid_simd<S, V> &in) {
  return in;
}

///////////////////////
// timesI
///////////////////////
template <class S, class V, IfComplex<S> = 0>
inline void timesI(Grid_simd<S, V> &ret, const Grid_simd<S, V> &in) {
  ret.v = binary<V>(in.v, ret.v, TimesISIMD());
}
template <class S, class V, IfComplex<S> = 0>
inline Grid_simd<S, V> timesI(const Grid_simd<S, V> &in) {
  Grid_simd<S, V> ret;
  timesI(ret, in);
  return ret;
}
template <class S, class V, IfNotComplex<S> = 0>
inline Grid_simd<S, V> timesI(const Grid_simd<S, V> &in) {
  return in;
}

/////////////////////
// Inner, outer
/////////////////////

template <class S, class V>
inline Grid_simd<S, V> innerProduct(const Grid_simd<S, V> &l,
                                    const Grid_simd<S, V> &r) {
  return conjugate(l) * r;
}
template <class S, class V>
inline Grid_simd<S, V> outerProduct(const Grid_simd<S, V> &l,
                                    const Grid_simd<S, V> &r) {
  return l * conjugate(r);
}

template <class S, class V>
inline Grid_simd<S, V> trace(const Grid_simd<S, V> &arg) {
  return arg;
}

////////////////////////////////////////////////////////////
// copy/splat complex real parts into real;
// insert real into complex and zero imag;
////////////////////////////////////////////////////////////

// real = toReal( complex )
template <class S, class V, IfReal<S> = 0>
inline Grid_simd<S, V> toReal(const Grid_simd<std::complex<S>, V> &in) {
  typedef Grid_simd<S, V> simd;
  simd ret;
  typename simd::conv_t conv;
  conv.v = in.v;  // copy the vector content (bytewise)
  for (int i = 0; i < simd::Nsimd(); i += 2) {
    conv.s[i + 1] = conv.s[i];  // duplicate (r,r);(r,r);(r,r); etc...
  }
  ret.v = conv.v;
  return ret;
}

// complex = toComplex( real )
template <class R, class V, IfReal<R> = 0>  // must be a real arg
inline Grid_simd<std::complex<R>, V> toComplex(const Grid_simd<R, V> &in) {
  typedef Grid_simd<R, V> Rsimd;
  typedef Grid_simd<std::complex<R>, V> Csimd;
  typename Rsimd::conv_t conv;  // address as real

  conv.v = in.v;
  for (int i = 0; i < Rsimd::Nsimd(); i += 2) {
    assert(conv.s[i + 1] == conv.s[i]);  
    // trap any cases where real was not duplicated
    // indicating the SIMD grids of real and imag assignment did not correctly
    // match
    conv.s[i + 1] = 0.0;  // zero imaginary parts
  }
  Csimd ret;
  ret.v = conv.v;
  return ret;
}

///////////////////////////////
// Define available types
///////////////////////////////
typedef Grid_simd<float, SIMD_Ftype> vRealF;
typedef Grid_simd<double, SIMD_Dtype> vRealD;
typedef Grid_simd<std::complex<float>, SIMD_Ftype> vComplexF;
typedef Grid_simd<std::complex<double>, SIMD_Dtype> vComplexD;
typedef Grid_simd<Integer, SIMD_Itype> vInteger;

// Half precision; no arithmetic support
typedef Grid_simd<uint16_t, SIMD_Htype>               vRealH;
typedef Grid_simd<std::complex<uint16_t>, SIMD_Htype> vComplexH;

inline void precisionChange(vRealF    *out,vRealD    *in,int nvec)
{
  assert((nvec&0x1)==0);
  for(int m=0;m*2<nvec;m++){
    int n=m*2;
    out[m].v=Optimization::PrecisionChange::DtoS(in[n].v,in[n+1].v);
  }
}
inline void precisionChange(vRealH    *out,vRealD    *in,int nvec)
{
  assert((nvec&0x3)==0);
  for(int m=0;m*4<nvec;m++){
    int n=m*4;
    out[m].v=Optimization::PrecisionChange::DtoH(in[n].v,in[n+1].v,in[n+2].v,in[n+3].v);
  }
}
inline void precisionChange(vRealH    *out,vRealF    *in,int nvec)
{
  assert((nvec&0x1)==0);
  for(int m=0;m*2<nvec;m++){
    int n=m*2;
    out[m].v=Optimization::PrecisionChange::StoH(in[n].v,in[n+1].v);
  }
}
inline void precisionChange(vRealD    *out,vRealF    *in,int nvec)
{
  assert((nvec&0x1)==0);
  for(int m=0;m*2<nvec;m++){
    int n=m*2;
    Optimization::PrecisionChange::StoD(in[m].v,out[n].v,out[n+1].v);
  }
}
inline void precisionChange(vRealD    *out,vRealH    *in,int nvec)
{
  assert((nvec&0x3)==0);
  for(int m=0;m*4<nvec;m++){
    int n=m*4;
    Optimization::PrecisionChange::HtoD(in[m].v,out[n].v,out[n+1].v,out[n+2].v,out[n+3].v);
  }
}
inline void precisionChange(vRealF    *out,vRealH    *in,int nvec)
{
  assert((nvec&0x1)==0);
  for(int m=0;m*2<nvec;m++){
    int n=m*2;
    Optimization::PrecisionChange::HtoS(in[m].v,out[n].v,out[n+1].v);
  }
}
inline void precisionChange(vComplexF *out,vComplexD *in,int nvec){ precisionChange((vRealF *)out,(vRealD *)in,nvec);}
inline void precisionChange(vComplexH *out,vComplexD *in,int nvec){ precisionChange((vRealH *)out,(vRealD *)in,nvec);}
inline void precisionChange(vComplexH *out,vComplexF *in,int nvec){ precisionChange((vRealH *)out,(vRealF *)in,nvec);}
inline void precisionChange(vComplexD *out,vComplexF *in,int nvec){ precisionChange((vRealD *)out,(vRealF *)in,nvec);}
inline void precisionChange(vComplexD *out,vComplexH *in,int nvec){ precisionChange((vRealD *)out,(vRealH *)in,nvec);}
inline void precisionChange(vComplexF *out,vComplexH *in,int nvec){ precisionChange((vRealF *)out,(vRealH *)in,nvec);}

// Check our vector types are of an appropriate size.
#if defined QPX
static_assert(2*sizeof(SIMD_Ftype) == sizeof(SIMD_Dtype), "SIMD vector lengths incorrect");
static_assert(2*sizeof(SIMD_Ftype) == sizeof(SIMD_Itype), "SIMD vector lengths incorrect");
#else
static_assert(sizeof(SIMD_Ftype) == sizeof(SIMD_Dtype), "SIMD vector lengths incorrect");
static_assert(sizeof(SIMD_Ftype) == sizeof(SIMD_Itype), "SIMD vector lengths incorrect");
#endif

/////////////////////////////////////////
// Some traits to recognise the types
/////////////////////////////////////////
template <typename T>
struct is_simd : public std::false_type {};
template <> struct is_simd<vRealF>     : public std::true_type {};
template <> struct is_simd<vRealD>     : public std::true_type {};
template <> struct is_simd<vComplexF>  : public std::true_type {};
template <> struct is_simd<vComplexD>  : public std::true_type {};
template <> struct is_simd<vInteger>   : public std::true_type {};

template <typename T> using IfSimd    = Invoke<std::enable_if<is_simd<T>::value, int> >;
template <typename T> using IfNotSimd = Invoke<std::enable_if<!is_simd<T>::value, unsigned> >;
}

#endif
