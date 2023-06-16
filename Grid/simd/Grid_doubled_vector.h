/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./lib/simd/Grid_vector_types.h

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>

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
#pragma once

NAMESPACE_BEGIN(Grid);

template <class Scalar_type, class Vector_type>
class Grid_simd2 {
public:
  typedef typename RealPart<Scalar_type>::type Real;
  typedef Vector_type vector_type;
  typedef Scalar_type scalar_type;

  typedef union conv_t_union {
    Vector_type v;
    Scalar_type s[sizeof(Vector_type) / sizeof(Scalar_type)];
    accelerator_inline conv_t_union(){};
  } conv_t;

  static constexpr int nvec=2;
  Vector_type v[nvec];

  static accelerator_inline constexpr int Nsimd(void) {
    static_assert( (sizeof(Vector_type) / sizeof(Scalar_type) >= 1), " size mismatch " );
    
    return nvec*sizeof(Vector_type) / sizeof(Scalar_type);
  }
  
  accelerator_inline Grid_simd2 &operator=(const Grid_simd2 &&rhs) {
    for(int n=0;n<nvec;n++) v[n] = rhs.v[n];
    return *this;
  };
  accelerator_inline Grid_simd2 &operator=(const Grid_simd2 &rhs) {
    for(int n=0;n<nvec;n++) v[n] = rhs.v[n];
    return *this;
  };  // faster than not declaring it and leaving to the compiler

  accelerator Grid_simd2() = default;
  accelerator_inline Grid_simd2(const Grid_simd2 &rhs) {    for(int n=0;n<nvec;n++) v[n] = rhs.v[n];  };
  accelerator_inline Grid_simd2(const Grid_simd2 &&rhs){    for(int n=0;n<nvec;n++) v[n] = rhs.v[n];  };
  accelerator_inline Grid_simd2(const Real a) { vsplat(*this, Scalar_type(a)); };
  // Enable if complex type
  template <typename S = Scalar_type> accelerator_inline
  Grid_simd2(const typename std::enable_if<is_complex<S>::value, S>::type a) {
      vsplat(*this, a);
  };

  /////////////////////////////
  // Constructors
  /////////////////////////////
  accelerator_inline Grid_simd2 &  operator=(const Zero &z) {
    vzero(*this);
    return (*this);
  }

  ///////////////////////////////////////////////
  // mac, mult, sub, add, adj
  ///////////////////////////////////////////////

  friend accelerator_inline void mac(Grid_simd2 *__restrict__ y,
				     const Grid_simd2 *__restrict__ a,
				     const Grid_simd2 *__restrict__ x) {
    *y = (*a) * (*x) + (*y);
  };

  friend accelerator_inline void mult(Grid_simd2 *__restrict__ y,
				      const Grid_simd2 *__restrict__ l,
				      const Grid_simd2 *__restrict__ r) {
    *y = (*l) * (*r);
  }

  friend accelerator_inline void sub(Grid_simd2 *__restrict__ y,
				     const Grid_simd2 *__restrict__ l,
				     const Grid_simd2 *__restrict__ r) {
    *y = (*l) - (*r);
  }
  friend accelerator_inline void add(Grid_simd2 *__restrict__ y,
				     const Grid_simd2 *__restrict__ l,
				     const Grid_simd2 *__restrict__ r) {
    *y = (*l) + (*r);
  }
  friend accelerator_inline void mac(Grid_simd2 *__restrict__ y,
				     const Scalar_type *__restrict__ a,
				     const Grid_simd2 *__restrict__ x) {
    *y = (*a) * (*x) + (*y);
  };
  friend accelerator_inline void mult(Grid_simd2 *__restrict__ y,
				      const Scalar_type *__restrict__ l,
				      const Grid_simd2 *__restrict__ r) {
    *y = (*l) * (*r);
  }
  friend accelerator_inline void sub(Grid_simd2 *__restrict__ y,
				     const Scalar_type *__restrict__ l,
				     const Grid_simd2 *__restrict__ r) {
    *y = (*l) - (*r);
  }
  friend accelerator_inline void add(Grid_simd2 *__restrict__ y,
				     const Scalar_type *__restrict__ l,
				     const Grid_simd2 *__restrict__ r) {
    *y = (*l) + (*r);
  }

  friend accelerator_inline void mac(Grid_simd2 *__restrict__ y,
				     const Grid_simd2 *__restrict__ a,
				     const Scalar_type *__restrict__ x) {
    *y = (*a) * (*x) + (*y);
  };
  friend accelerator_inline void mult(Grid_simd2 *__restrict__ y,
				      const Grid_simd2 *__restrict__ l,
				      const Scalar_type *__restrict__ r) {
    *y = (*l) * (*r);
  }
  friend accelerator_inline void sub(Grid_simd2 *__restrict__ y,
				     const Grid_simd2 *__restrict__ l,
				     const Scalar_type *__restrict__ r) {
    *y = (*l) - (*r);
  }
  friend accelerator_inline void add(Grid_simd2 *__restrict__ y,
				     const Grid_simd2 *__restrict__ l,
				     const Scalar_type *__restrict__ r) {
    *y = (*l) + (*r);
  }

  ////////////////////////////////////////////////////////////////////////
  // FIXME:  gonna remove these load/store, get, set, prefetch
  ////////////////////////////////////////////////////////////////////////
  friend accelerator_inline void vset(Grid_simd2 &ret, Scalar_type *a) {
    for(int n=0;n<nvec;n++) vset(ret.v[n],a);
  }

  ///////////////////////
  // Vstore
  ///////////////////////
  friend accelerator_inline void vstore(const Grid_simd2 &ret, Scalar_type *a) {
    for(int n=0;n<nvec;n++) vstore(ret.v[n],a);
  }

  ///////////////////////
  // Vprefetch
  ///////////////////////
  friend accelerator_inline void vprefetch(const Grid_simd2 &v) {
    vprefetch(v.v[0]);
  }

  ///////////////////////
  // Reduce
  ///////////////////////
  friend accelerator_inline Scalar_type Reduce(const Grid_simd2 &in) {
    return Reduce(in.v[0])+ Reduce(in.v[1]);
  }

  ////////////////////////////
  // operator scalar * simd
  ////////////////////////////
  friend accelerator_inline Grid_simd2 operator*(const Scalar_type &a, Grid_simd2 b) {
    Grid_simd2 va;
    vsplat(va, a);
    return va * b;
  }
  friend accelerator_inline Grid_simd2 operator*(Grid_simd2 b, const Scalar_type &a) {
    return a * b;
  }

  //////////////////////////////////
  // Divides
  //////////////////////////////////
  friend accelerator_inline Grid_simd2 operator/(const Scalar_type &a, Grid_simd2 b) {
    Grid_simd2 va;
    vsplat(va, a);
    return va / b;
  }
  friend accelerator_inline Grid_simd2 operator/(Grid_simd2 b, const Scalar_type &a) {
    Grid_simd2 va;
    vsplat(va, a);
    return b / a;
  }

  ///////////////////////
  // Unary negation
  ///////////////////////
  friend accelerator_inline Grid_simd2 operator-(const Grid_simd2 &r) {
    Grid_simd2 ret;
    vzero(ret);
    ret = ret - r;
    return ret;
  }
  // *=,+=,-= operators
  accelerator_inline Grid_simd2 &operator*=(const Grid_simd2 &r) {
    *this = (*this) * r;
    return *this;
  }
  accelerator_inline Grid_simd2 &operator+=(const Grid_simd2 &r) {
    *this = *this + r;
    return *this;
  }
  accelerator_inline Grid_simd2 &operator-=(const Grid_simd2 &r) {
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
  friend accelerator_inline Grid_simd2 SimdApply(const functor &func, const Grid_simd2 &v) {
    Grid_simd2 ret;
    for(int n=0;n<nvec;n++){
      ret.v[n]=SimdApply(func,v.v[n]);
    }
    return ret;
  }
  template <class functor>
  friend accelerator_inline Grid_simd2 SimdApplyBinop(const functor &func,
                                         const Grid_simd2 &x,
                                         const Grid_simd2 &y) {
    Grid_simd2 ret;
    for(int n=0;n<nvec;n++){
      ret.v[n]=SimdApplyBinop(func,x.v[n],y.v[n]);
    }
    return ret;
  }
  ///////////////////////
  // Exchange
  // Al Ah , Bl Bh -> Al Bl Ah,Bh
  ///////////////////////
  friend accelerator_inline void exchange0(Grid_simd2 &out1,Grid_simd2 &out2,Grid_simd2 in1,Grid_simd2 in2){
      out1.v[0] = in1.v[0];
      out1.v[1] = in2.v[0];
      out2.v[0] = in1.v[1];
      out2.v[1] = in2.v[1];
  }
  friend accelerator_inline void exchange1(Grid_simd2 &out1,Grid_simd2 &out2,Grid_simd2 in1,Grid_simd2 in2){
    exchange0(out1.v[0],out2.v[0],in1.v[0],in2.v[0]);
    exchange0(out1.v[1],out2.v[1],in1.v[1],in2.v[1]);
  }
  friend accelerator_inline void exchange2(Grid_simd2 &out1,Grid_simd2 &out2,Grid_simd2 in1,Grid_simd2 in2){
    exchange1(out1.v[0],out2.v[0],in1.v[0],in2.v[0]);
    exchange1(out1.v[1],out2.v[1],in1.v[1],in2.v[1]);
  }
  friend accelerator_inline void exchange3(Grid_simd2 &out1,Grid_simd2 &out2,Grid_simd2 in1,Grid_simd2 in2){
    exchange2(out1.v[0],out2.v[0],in1.v[0],in2.v[0]);
    exchange2(out1.v[1],out2.v[1],in1.v[1],in2.v[1]);
  }
  friend accelerator_inline void exchange4(Grid_simd2 &out1,Grid_simd2 &out2,Grid_simd2 in1,Grid_simd2 in2){
    exchange3(out1.v[0],out2.v[0],in1.v[0],in2.v[0]);
    exchange3(out1.v[1],out2.v[1],in1.v[1],in2.v[1]);
  }
  friend accelerator_inline void exchange(Grid_simd2 &out1,Grid_simd2 &out2,Grid_simd2 in1,Grid_simd2 in2,int n)
  {
    if       (n==3) {
      exchange3(out1,out2,in1,in2);
    } else if(n==2) {
      exchange2(out1,out2,in1,in2);
    } else if(n==1) {
      exchange1(out1,out2,in1,in2);
    } else if(n==0) {
      exchange0(out1,out2,in1,in2);
    }
  }
  ////////////////////////////////////////////////////////////////////
  // General permute; assumes vector length is same across
  // all subtypes; may not be a good assumption, but could
  // add the vector width as a template param for BG/Q for example
  ////////////////////////////////////////////////////////////////////
  friend accelerator_inline void permute0(Grid_simd2 &y, Grid_simd2 b) {
    y.v[0]=b.v[1];
    y.v[1]=b.v[0];
  }
  friend accelerator_inline void permute1(Grid_simd2 &y, Grid_simd2 b) {
    permute0(y.v[0],b.v[0]);
    permute0(y.v[1],b.v[1]);
  }
  friend accelerator_inline void permute2(Grid_simd2 &y, Grid_simd2 b) {
    permute1(y.v[0],b.v[0]);
    permute1(y.v[1],b.v[1]);
  }
  friend accelerator_inline void permute3(Grid_simd2 &y, Grid_simd2 b) {
    permute2(y.v[0],b.v[0]);
    permute2(y.v[1],b.v[1]);
  }
  friend accelerator_inline void permute4(Grid_simd2 &y, Grid_simd2 b) {
    permute3(y.v[0],b.v[0]);
    permute3(y.v[1],b.v[1]);
  }
  friend accelerator_inline void permute(Grid_simd2 &y, Grid_simd2 b, int perm) {
    if(perm==3) permute3(y, b);
    else if(perm==2) permute2(y, b);
    else if(perm==1) permute1(y, b);
    else if(perm==0) permute0(y, b);
  }

  ///////////////////////////////
  // Getting single lanes
  ///////////////////////////////
  accelerator_inline Scalar_type getlane(int lane) const {
    if(lane < vector_type::Nsimd() ) return v[0].getlane(lane);
    else                             return v[1].getlane(lane%vector_type::Nsimd());
  }

  accelerator_inline void putlane(const Scalar_type &S, int lane){
    if(lane < vector_type::Nsimd() ) v[0].putlane(S,lane);
    else                             v[1].putlane(S,lane%vector_type::Nsimd());
  }
};  // end of Grid_simd2 class definition

///////////////////////////////
// Define available types
///////////////////////////////

typedef Grid_simd2<complex<double>  , vComplexD>  vComplexD2;
typedef Grid_simd2<double           , vRealD>     vRealD2;



/////////////////////////////////////////
// Some traits to recognise the types
/////////////////////////////////////////
template <typename T>
struct is_simd : public std::false_type {};
template <> struct is_simd<vRealF>     : public std::true_type {};
template <> struct is_simd<vRealD>     : public std::true_type {};
template <> struct is_simd<vRealH>     : public std::true_type {};
template <> struct is_simd<vComplexF>  : public std::true_type {};
template <> struct is_simd<vComplexD>  : public std::true_type {};
template <> struct is_simd<vComplexH>  : public std::true_type {};
template <> struct is_simd<vInteger>   : public std::true_type {};
template <> struct is_simd<vRealD2>    : public std::true_type {};
template <> struct is_simd<vComplexD2> : public std::true_type {};

template <typename T> using IfSimd    = Invoke<std::enable_if<is_simd<T>::value, int> >;
template <typename T> using IfNotSimd = Invoke<std::enable_if<!is_simd<T>::value, unsigned> >;

///////////////////////////////////////////////
// insert / extract with complex support
///////////////////////////////////////////////
template <class S, class V>
accelerator_inline S getlane(const Grid_simd<S, V> &in,int lane) {
  return in.getlane(lane);
}
template <class S, class V>
accelerator_inline void putlane(Grid_simd<S, V> &vec,const S &_S, int lane){
  vec.putlane(_S,lane);
}
template <class S,IfNotSimd<S> = 0 >
accelerator_inline S getlane(const S &in,int lane) {
  return in;
}
template <class S,IfNotSimd<S> = 0 >
accelerator_inline void putlane(S &vec,const S &_S, int lane){
  vec = _S;
}
template <class S, class V>
accelerator_inline S getlane(const Grid_simd2<S, V> &in,int lane) {
  return in.getlane(lane);
}
template <class S, class V>
accelerator_inline void putlane(Grid_simd2<S, V> &vec,const S &_S, int lane){
  vec.putlane(_S,lane);
}


////////////////////////////////////////////////////////////////////
// General rotate
////////////////////////////////////////////////////////////////////

template <class S, class V>
accelerator_inline void vbroadcast(Grid_simd2<S,V> &ret,const Grid_simd2<S,V> &src,int lane){
  S* typepun =(S*) &src;
  vsplat(ret,typepun[lane]);
}
template <class S, class V, IfComplex<S> =0>
accelerator_inline void rbroadcast(Grid_simd2<S,V> &ret,const Grid_simd2<S,V> &src,int lane){
  typedef typename V::vector_type vector_type;
  S* typepun =(S*) &src;
  ret.v[0].v = unary<vector_type>(real(typepun[lane]), VsplatSIMD());
  ret.v[1].v = unary<vector_type>(real(typepun[lane]), VsplatSIMD());
}


///////////////////////
// Splat
///////////////////////

// this is only for the complex version
template <class S, class V, IfComplex<S> = 0, class ABtype>
accelerator_inline void vsplat(Grid_simd2<S, V> &ret, ABtype a, ABtype b) {
  vsplat(ret.v[0],a,b);
  vsplat(ret.v[1],a,b);
}

// overload if complex
template <class S, class V>
accelerator_inline void vsplat(Grid_simd2<S, V> &ret, EnableIf<is_complex<S>, S> c) {
  vsplat(ret, real(c), imag(c));
}
template <class S, class V>
accelerator_inline void rsplat(Grid_simd2<S, V> &ret, EnableIf<is_complex<S>, S> c) {
  vsplat(ret, real(c), real(c));
}

// if real fill with a, if complex fill with a in the real part (first function
// above)
template <class S, class V>
accelerator_inline void vsplat(Grid_simd2<S, V> &ret, NotEnableIf<is_complex<S>, S> a)
{
  vsplat(ret.v[0],a);
  vsplat(ret.v[1],a);
}
//////////////////////////

///////////////////////////////////////////////
// Initialise to 1,0,i for the correct types
///////////////////////////////////////////////
// For complex types
template <class S, class V, IfComplex<S> = 0>
accelerator_inline void vone(Grid_simd2<S, V> &ret) {
  vsplat(ret, S(1.0, 0.0));
}
template <class S, class V, IfComplex<S> = 0>
accelerator_inline void vzero(Grid_simd2<S, V> &ret) {
  vsplat(ret, S(0.0, 0.0));
}  // use xor?
template <class S, class V, IfComplex<S> = 0>
accelerator_inline void vcomplex_i(Grid_simd2<S, V> &ret) {
  vsplat(ret, S(0.0, 1.0));
}

template <class S, class V, IfComplex<S> = 0>
accelerator_inline void visign(Grid_simd2<S, V> &ret) {
  vsplat(ret, S(1.0, -1.0));
}
template <class S, class V, IfComplex<S> = 0>
accelerator_inline void vrsign(Grid_simd2<S, V> &ret) {
  vsplat(ret, S(-1.0, 1.0));
}

// if not complex overload here
template <class S, class V, IfReal<S> = 0>
accelerator_inline void vone(Grid_simd2<S, V> &ret) {
  vsplat(ret, S(1.0));
}
template <class S, class V, IfReal<S> = 0>
accelerator_inline void vzero(Grid_simd2<S, V> &ret) {
  vsplat(ret, S(0.0));
}

// For integral types
template <class S, class V, IfInteger<S> = 0>
accelerator_inline void vone(Grid_simd2<S, V> &ret) {
  vsplat(ret, 1);
}
template <class S, class V, IfInteger<S> = 0>
accelerator_inline void vzero(Grid_simd2<S, V> &ret) {
  vsplat(ret, 0);
}
template <class S, class V, IfInteger<S> = 0>
accelerator_inline void vtrue(Grid_simd2<S, V> &ret) {
  vsplat(ret, 0xFFFFFFFF);
}
template <class S, class V, IfInteger<S> = 0>
accelerator_inline void vfalse(Grid_simd2<S, V> &ret) {
  vsplat(ret, 0);
}
template <class S, class V>
accelerator_inline void zeroit(Grid_simd2<S, V> &z) {
  vzero(z);
}

///////////////////////
// Vstream
///////////////////////
template <class S, class V, IfReal<S> = 0>
accelerator_inline void vstream(Grid_simd2<S, V> &out, const Grid_simd2<S, V> &in) {
  vstream(out.v[0],in.v[0]);
  vstream(out.v[1],in.v[1]);
}
template <class S, class V, IfComplex<S> = 0>
accelerator_inline void vstream(Grid_simd2<S, V> &out, const Grid_simd2<S, V> &in) {
  vstream(out.v[0],in.v[0]);
  vstream(out.v[1],in.v[1]);
}
template <class S, class V, IfInteger<S> = 0>
accelerator_inline void vstream(Grid_simd2<S, V> &out, const Grid_simd2<S, V> &in) {
  vstream(out.v[0],in.v[0]);
  vstream(out.v[1],in.v[1]);
}

////////////////////////////////////
// Arithmetic operator overloads +,-,*
////////////////////////////////////
template <class S, class V>
accelerator_inline Grid_simd2<S, V> operator+(Grid_simd2<S, V> a, Grid_simd2<S, V> b) {
  Grid_simd2<S, V> ret;
  ret.v[0] = a.v[0]+b.v[0];
  ret.v[1] = a.v[1]+b.v[1];
  return ret;
};

template <class S, class V>
accelerator_inline Grid_simd2<S, V> operator-(Grid_simd2<S, V> a, Grid_simd2<S, V> b) {
  Grid_simd2<S, V> ret;
  ret.v[0] = a.v[0]-b.v[0];
  ret.v[1] = a.v[1]-b.v[1];
  return ret;
};

// Distinguish between complex types and others
template <class S, class V, IfComplex<S> = 0>
accelerator_inline Grid_simd2<S, V> real_mult(Grid_simd2<S, V> a, Grid_simd2<S, V> b) {
  Grid_simd2<S, V> ret;
  ret.v[0] =real_mult(a.v[0],b.v[0]);
  ret.v[1] =real_mult(a.v[1],b.v[1]);
  return ret;
};
template <class S, class V, IfComplex<S> = 0>
accelerator_inline Grid_simd2<S, V> real_madd(Grid_simd2<S, V> a, Grid_simd2<S, V> b, Grid_simd2<S,V> c) {
  Grid_simd2<S, V> ret;
  ret.v[0] =real_madd(a.v[0],b.v[0],c.v[0]);
  ret.v[1] =real_madd(a.v[1],b.v[1],c.v[1]);
  return ret;
};


// Distinguish between complex types and others
template <class S, class V>
accelerator_inline Grid_simd2<S, V> operator*(Grid_simd2<S, V> a, Grid_simd2<S, V> b) {
  Grid_simd2<S, V> ret;
  ret.v[0] = a.v[0]*b.v[0];
  ret.v[1] = a.v[1]*b.v[1];
  return ret;
};

// Distinguish between complex types and others
template <class S, class V>
accelerator_inline Grid_simd2<S, V> operator/(Grid_simd2<S, V> a, Grid_simd2<S, V> b) {
  Grid_simd2<S, V> ret;
  ret.v[0] = a.v[0]/b.v[0];
  ret.v[1] = a.v[1]/b.v[1];
  return ret;
};

///////////////////////
// Conjugate
///////////////////////
template <class S, class V>
accelerator_inline Grid_simd2<S, V> conjugate(const Grid_simd2<S, V> &in) {
  Grid_simd2<S, V> ret;
  ret.v[0] = conjugate(in.v[0]);
  ret.v[1] = conjugate(in.v[1]);
  return ret;
}
template <class S, class V, IfNotInteger<S> = 0>
accelerator_inline Grid_simd2<S, V> adj(const Grid_simd2<S, V> &in) {
  return conjugate(in);
}

///////////////////////
// timesMinusI
///////////////////////
template <class S, class V>
accelerator_inline void timesMinusI(Grid_simd2<S, V> &ret, const Grid_simd2<S, V> &in) {
  timesMinusI(ret.v[0],in.v[0]);
  timesMinusI(ret.v[1],in.v[1]);
}
template <class S, class V>
accelerator_inline Grid_simd2<S, V> timesMinusI(const Grid_simd2<S, V> &in) {
  Grid_simd2<S, V> ret;
  timesMinusI(ret.v[0],in.v[0]);
  timesMinusI(ret.v[1],in.v[1]);
  return ret;
}

///////////////////////
// timesI
///////////////////////
template <class S, class V>
accelerator_inline void timesI(Grid_simd2<S, V> &ret, const Grid_simd2<S, V> &in) {
  timesI(ret.v[0],in.v[0]);
  timesI(ret.v[1],in.v[1]);
}
template <class S, class V>
accelerator_inline Grid_simd2<S, V> timesI(const Grid_simd2<S, V> &in) {
  Grid_simd2<S, V> ret;
  timesI(ret.v[0],in.v[0]);
  timesI(ret.v[1],in.v[1]);
  return ret;
}

/////////////////////
// Inner, outer
/////////////////////
template <class S, class V>
accelerator_inline Grid_simd2<S, V> innerProduct(const Grid_simd2<S, V> &l,const Grid_simd2<S, V> &r) {
  return conjugate(l) * r;
}
template <class S, class V>
accelerator_inline Grid_simd2<S, V> outerProduct(const Grid_simd2<S, V> &l,const Grid_simd2<S, V> &r) {
  return l * conjugate(r);
}

template <class S, class V>
accelerator_inline Grid_simd2<S, V> trace(const Grid_simd2<S, V> &arg) {
  return arg;
}

////////////////////////////////////////////////////////////
// copy/splat complex real parts into real;
// insert real into complex and zero imag;
////////////////////////////////////////////////////////////
accelerator_inline void precisionChange(vComplexD2 &out,const vComplexF  &in){
  Optimization::PrecisionChange::StoD(in.v,out.v[0].v,out.v[1].v);
}
accelerator_inline void precisionChange(vComplexF  &out,const vComplexD2 &in){
  out.v=Optimization::PrecisionChange::DtoS(in.v[0].v,in.v[1].v);
}
accelerator_inline void precisionChange(vComplexD2 *out,const vComplexF  *in,int nvec){
  for(int m=0;m<nvec;m++){ precisionChange(out[m],in[m]); }
}
accelerator_inline void precisionChange(vComplexF  *out,const vComplexD2 *in,int nvec){
  for(int m=0;m<nvec;m++){ precisionChange(out[m],in[m]); }
}

accelerator_inline void precisionChange(vRealD2 &out,const vRealF  &in){
  Optimization::PrecisionChange::StoD(in.v,out.v[0].v,out.v[1].v);
}
accelerator_inline void precisionChange(vRealF  &out,const vRealD2 &in){
  out.v=Optimization::PrecisionChange::DtoS(in.v[0].v,in.v[1].v);
}
accelerator_inline void precisionChange(vRealD2 *out,const vRealF  *in,int nvec){
  for(int m=0;m<nvec;m++){ precisionChange(out[m],in[m]); }
}
accelerator_inline void precisionChange(vRealF  *out,const vRealD2 *in,int nvec){
  for(int m=0;m<nvec;m++){ precisionChange(out[m],in[m]); }
}

NAMESPACE_END(Grid);


