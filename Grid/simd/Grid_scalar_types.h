#pragma once
NAMESPACE_BEGIN(Grid);
template <class Scalar_type>
class Grid_simd1 {
public:
  typedef typename RealPart<Scalar_type>::type Real;
  typedef Scalar_type vector_type;
  typedef Scalar_type scalar_type;

  vector_type v;

  static accelerator_inline constexpr int Nsimd(void) {
    return 1;
  }
  accelerator_inline Grid_simd1 &operator=(const Grid_simd1 &&rhs) {
    v = rhs.v;
    return *this;
  };
  accelerator_inline Grid_simd1 &operator=(const Grid_simd1 &rhs) {
    v = rhs.v;
    return *this;
  };  // faster than not declaring it and leaving to the compiler

  accelerator Grid_simd1() = default;
  accelerator_inline Grid_simd1(const Grid_simd1 &rhs)  : v(rhs.v){};  // compiles in movaps
  accelerator_inline Grid_simd1(const Grid_simd1 &&rhs) : v(rhs.v){};

  accelerator_inline Grid_simd1(const Real a) { v=Scalar_type(a); };

  template <typename S = Scalar_type> accelerator_inline
  Grid_simd1(const typename std::enable_if<is_complex<S>::value, S>::type a) {
    v=Scalar_type(a);
  };

  /////////////////////////////
  // Constructors
  /////////////////////////////
  accelerator_inline Grid_simd1 &  operator=(const Zero &z) {
    v=scalar_type(0);
    return *this;
  }

  ///////////////////////////////////////////////
  // mac, mult, sub, add, adj
  ///////////////////////////////////////////////

  friend accelerator_inline void mac(Grid_simd1 *__restrict__ y,
				     const Grid_simd1 *__restrict__ a,
				     const Grid_simd1 *__restrict__ x) {
    *y = (*a) * (*x) + (*y);
  };

  friend accelerator_inline void mult(Grid_simd1 *__restrict__ y,
				      const Grid_simd1 *__restrict__ l,
				      const Grid_simd1 *__restrict__ r) {
    *y = (*l) * (*r);
  }

  friend accelerator_inline void sub(Grid_simd1 *__restrict__ y,
				     const Grid_simd1 *__restrict__ l,
				     const Grid_simd1 *__restrict__ r) {
    *y = (*l) - (*r);
  }
  friend accelerator_inline void add(Grid_simd1 *__restrict__ y,
				     const Grid_simd1 *__restrict__ l,
				     const Grid_simd1 *__restrict__ r) {
    *y = (*l) + (*r);
  }
  friend accelerator_inline void mac(Grid_simd1 *__restrict__ y,
				     const Scalar_type *__restrict__ a,
				     const Grid_simd1 *__restrict__ x) {
    *y = (*a) * (*x) + (*y);
  };
  friend accelerator_inline void mult(Grid_simd1 *__restrict__ y,
				      const Scalar_type *__restrict__ l,
				      const Grid_simd1 *__restrict__ r) {
    *y = (*l) * (*r);
  }
  friend accelerator_inline void sub(Grid_simd1 *__restrict__ y,
				     const Scalar_type *__restrict__ l,
				     const Grid_simd1 *__restrict__ r) {
    *y = (*l) - (*r);
  }
  friend accelerator_inline void add(Grid_simd1 *__restrict__ y,
				     const Scalar_type *__restrict__ l,
				     const Grid_simd1 *__restrict__ r) {
    *y = (*l) + (*r);
  }

  friend accelerator_inline void mac(Grid_simd1 *__restrict__ y,
				     const Grid_simd1 *__restrict__ a,
				     const Scalar_type *__restrict__ x) {
    *y = (*a) * (*x) + (*y);
  };
  friend accelerator_inline void mult(Grid_simd1 *__restrict__ y,
				      const Grid_simd1 *__restrict__ l,
				      const Scalar_type *__restrict__ r) {
    *y = (*l) * (*r);
  }
  friend accelerator_inline void sub(Grid_simd1 *__restrict__ y,
				     const Grid_simd1 *__restrict__ l,
				     const Scalar_type *__restrict__ r) {
    *y = (*l) - (*r);
  }
  friend accelerator_inline void add(Grid_simd1 *__restrict__ y,
				     const Grid_simd1 *__restrict__ l,
				     const Scalar_type *__restrict__ r) {
    *y = (*l) + (*r);
  }

  ////////////////////////////////////////////////////////////////////////
  // FIXME:  gonna remove these load/store, get, set, prefetch
  ////////////////////////////////////////////////////////////////////////
  friend accelerator_inline void vset(Grid_simd1 &ret, Scalar_type *a) {
    ret.v = *a;
  }

  ///////////////////////
  // Vstore
  ///////////////////////
  friend accelerator_inline void vstore(const Grid_simd1 &ret, Scalar_type *a) {
    *a=ret.v;
  }

  ///////////////////////
  // Vprefetch
  ///////////////////////
  friend accelerator_inline void vprefetch(const Grid_simd1 &v) { }

  ///////////////////////
  // Reduce
  ///////////////////////
  friend accelerator_inline Scalar_type Reduce(const Grid_simd1 &in) {
    return in.v;
  }
  ////////////////////////////
  // operator scalar * simd
  ////////////////////////////
  friend accelerator_inline Grid_simd1 operator*(const Scalar_type &a, Grid_simd1 b) {
    Grid_simd1 va;
    va.v=a;
    return va * b;
  }
  friend accelerator_inline Grid_simd1 operator*(Grid_simd1 b, const Scalar_type &a) {
    return a * b;
  }

  //////////////////////////////////
  // Divides
  //////////////////////////////////
  friend accelerator_inline Grid_simd1 operator/(const Scalar_type &a, Grid_simd1 b) {
    Grid_simd1 va;
    va.v = a;
    return va / b;
  }
  friend accelerator_inline Grid_simd1 operator/(Grid_simd1 b, const Scalar_type &a) {
    Grid_simd1 va;
    va.v=a;
    return b / va;
  }
  ///////////////////////
  // Unary negation
  ///////////////////////
  friend accelerator_inline Grid_simd1 operator-(const Grid_simd1 &r) {
    Grid_simd1 ret;
    ret.v = scalar_type(0);
    ret = ret - r;
    return ret;
  }
  // *=,+=,-= operators
  accelerator_inline Grid_simd1 &operator*=(const Grid_simd1 &r) {
    *this = (*this) * r;
    return *this;
  }
  accelerator_inline Grid_simd1 &operator+=(const Grid_simd1 &r) {
    *this = *this + r;
    return *this;
  }
  accelerator_inline Grid_simd1 &operator-=(const Grid_simd1 &r) {
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
  friend accelerator_inline Grid_simd1 SimdApply(const functor &func, const Grid_simd1 &v) {
    Grid_simd1 ret;
    ret.v = func(v.v);
    return ret;
  }
  template <class functor>
  friend accelerator_inline Grid_simd1 SimdApplyBinop(const functor &func,
                                         const Grid_simd1 &x,
                                         const Grid_simd1 &y) {
    Grid_simd1 ret;
    ret.v = func(x.v,y.v);
    return ret;
  }
  ///////////////////////
  // Exchange
  // Al Ah , Bl Bh -> Al Bl Ah,Bh
  ///////////////////////
  friend accelerator_inline void exchange(Grid_simd1 &out1,Grid_simd1 &out2,Grid_simd1 in1,Grid_simd1 in2,int n)
  {
    assert(0);
  }
  friend accelerator_inline void exchange0(Grid_simd1 &out1,Grid_simd1 &out2,Grid_simd1 in1,Grid_simd1 in2){
    assert(0);
  }
  friend accelerator_inline void exchange1(Grid_simd1 &out1,Grid_simd1 &out2,Grid_simd1 in1,Grid_simd1 in2){
    assert(0);
  }
  friend accelerator_inline void exchange2(Grid_simd1 &out1,Grid_simd1 &out2,Grid_simd1 in1,Grid_simd1 in2){
    assert(0);
  }
  friend accelerator_inline void exchange3(Grid_simd1 &out1,Grid_simd1 &out2,Grid_simd1 in1,Grid_simd1 in2){
    assert(0);
  }
  ////////////////////////////////////////////////////////////////////
  // Permute: unreachable at Nsimd=1
  ////////////////////////////////////////////////////////////////////
  friend accelerator_inline void permute0(Grid_simd1 &y, Grid_simd1 b) {
    assert(0);
  }
  friend accelerator_inline void permute1(Grid_simd1 &y, Grid_simd1 b) {
    assert(0);
  }
  friend accelerator_inline void permute2(Grid_simd1 &y, Grid_simd1 b) {
    assert(0);
  }
  friend accelerator_inline void permute3(Grid_simd1 &y, Grid_simd1 b) {
    assert(0);
  }
  friend accelerator_inline void permute(Grid_simd1 &y, Grid_simd1 b, int perm) {
    assert(0);
  }

  ///////////////////////////////
  // Getting single lanes
  ///////////////////////////////
  accelerator_inline Scalar_type getlane(int lane) const {
    return v;
  }
  accelerator_inline void putlane(const Scalar_type &S, int lane){
    v=S;
  }

};
template <class Scalar>
inline std::ostream& operator<< (std::ostream& stream, const Grid_simd1<Scalar> &o){
  stream<<"<"<<o.v<<">";
  return stream;
}

typedef Grid_simd1<RealF> sRealF;
typedef Grid_simd1<RealD> sRealD;
typedef Grid_simd1<ComplexF> sComplexF;
typedef Grid_simd1<ComplexD> sComplexD;
typedef Grid_simd1<Integer>  sInteger;

/////////////////////////////////////////
// Permute
/////////////////////////////////////////

//accelerator_inline void permute(sComplexD &y,sComplexD b, int perm) {  y=b; }
//accelerator_inline void permute(sComplexF &y,sComplexF b, int perm) {  y=b; }
//accelerator_inline void permute(sRealD &y,sRealD b, int perm) {  y=b; }
//accelerator_inline void permute(sRealF &y,sRealF b, int perm) {  y=b; }

////////////////////////////////////////////////////////////////////
// General rotate
////////////////////////////////////////////////////////////////////
template <class S, IfNotComplex<S> = 0>
accelerator_inline Grid_simd1<S> rotate(Grid_simd1<S> b, int nrot) { return b; }
template <class S, IfComplex<S> = 0>
accelerator_inline Grid_simd1<S> rotate(Grid_simd1<S> b, int nrot) {  return b; }
template <class S, IfNotComplex<S> =0>
accelerator_inline void rotate( Grid_simd1<S> &ret,Grid_simd1<S> b,int nrot)
{
  ret = b;
}
template <class S, IfComplex<S> =0>
accelerator_inline void rotate(Grid_simd1<S> &ret,Grid_simd1<S> b,int nrot)
{
  ret = b;
}

template <class S>
accelerator_inline void vbroadcast(Grid_simd1<S> &ret,const Grid_simd1<S> &src,int lane){
  ret = src;
}
template <class S,  IfComplex<S> =0>
accelerator_inline void rbroadcast(Grid_simd1<S> &ret,const Grid_simd1<S> &src,int lane){
  ret.v = real(src.v);
}

///////////////////////
// Splat
///////////////////////

// this is only for the complex version
template <class S,  IfComplex<S> = 0, class ABtype>
accelerator_inline void vsplat(Grid_simd1<S> &ret, ABtype a, ABtype b) {
  ret.v = S(a,b);
}

// overload if complex
template <class S>
accelerator_inline void vsplat(Grid_simd1<S> &ret, EnableIf<is_complex<S>, S> c) {
  vsplat(ret, real(c), imag(c));
}
template <class S>
accelerator_inline void rsplat(Grid_simd1<S> &ret, EnableIf<is_complex<S>, S> c) {
  vsplat(ret, real(c), real(c));
}
// if real fill with a, if complex fill with a in the real part (first function
// above)
template <class S>
accelerator_inline void vsplat(Grid_simd1<S> &ret, NotEnableIf<is_complex<S>, S> a) {
  ret.v = a;
}
//////////////////////////


///////////////////////////////////////////////
// Initialise to 1,0,i for the correct types
///////////////////////////////////////////////
// For complex types
template <class S,  IfComplex<S> = 0>
accelerator_inline void vone(Grid_simd1<S> &ret) {
  vsplat(ret, S(1.0, 0.0));
}
template <class S,  IfComplex<S> = 0>
accelerator_inline void vzero(Grid_simd1<S> &ret) {
  vsplat(ret, S(0.0, 0.0));
}  // use xor?
template <class S,  IfComplex<S> = 0>
accelerator_inline void vcomplex_i(Grid_simd1<S> &ret) {
  vsplat(ret, S(0.0, 1.0));
}

template <class S,  IfComplex<S> = 0>
accelerator_inline void visign(Grid_simd1<S> &ret) {
  vsplat(ret, S(1.0, -1.0));
}
template <class S,  IfComplex<S> = 0>
accelerator_inline void vrsign(Grid_simd1<S> &ret) {
  vsplat(ret, S(-1.0, 1.0));
}

// if not complex overload here
template <class S,  IfReal<S> = 0>
accelerator_inline void vone(Grid_simd1<S> &ret) {
  vsplat(ret, S(1.0));
}
template <class S,  IfReal<S> = 0>
accelerator_inline void vzero(Grid_simd1<S> &ret) {
  vsplat(ret, S(0.0));
}

// For integral types
template <class S,  IfInteger<S> = 0>
accelerator_inline void vone(Grid_simd1<S> &ret) {
  vsplat(ret, 1);
}
template <class S,  IfInteger<S> = 0>
accelerator_inline void vzero(Grid_simd1<S> &ret) {
  vsplat(ret, 0);
}
template <class S,  IfInteger<S> = 0>
accelerator_inline void vtrue(Grid_simd1<S> &ret) {
  vsplat(ret, 0xFFFFFFFF);
}
template <class S,  IfInteger<S> = 0>
accelerator_inline void vfalse(Grid_simd1<S> &ret) {
  vsplat(ret, 0);
}
template <class S>
accelerator_inline void zeroit(Grid_simd1<S> &z) {
  vzero(z);
}

///////////////////////
// Vstream
///////////////////////
template <class S,  IfReal<S> = 0>
accelerator_inline void vstream(Grid_simd1<S> &out, const Grid_simd1<S> &in) {
  out = in;
}
template <class S,  IfComplex<S> = 0>
accelerator_inline void vstream(Grid_simd1<S> &out, const Grid_simd1<S> &in) {
  out = in;
}
template <class S,  IfInteger<S> = 0>
accelerator_inline void vstream(Grid_simd1<S> &out, const Grid_simd1<S> &in) {
  out = in;
}

////////////////////////////////////
// Arithmetic operator overloads +,-,*
////////////////////////////////////
template <class S>
accelerator_inline Grid_simd1<S> operator+(Grid_simd1<S> a, Grid_simd1<S> b) {
  Grid_simd1<S> ret;
  ret.v = a.v+b.v;
  return ret;
};

template <class S>
accelerator_inline Grid_simd1<S> operator-(Grid_simd1<S> a, Grid_simd1<S> b) {
  Grid_simd1<S> ret;
  ret.v = a.v-b.v;
  return ret;
};

// Distinguish between complex types and others
template <class S,  IfComplex<S> = 0>
accelerator_inline Grid_simd1<S> real_mult(Grid_simd1<S> a, Grid_simd1<S> b) {
  Grid_simd1<S> ret;
  ret.v = S(real(a.v)*real(b.v),real(a.v)*imag(b.v));
  return ret;
};
template <class S,  IfComplex<S> = 0>
accelerator_inline Grid_simd1<S> real_madd(Grid_simd1<S> a, Grid_simd1<S> b, Grid_simd1<S> c) {
  Grid_simd1<S> ret;
  ret = real_mult(a,b) + c;
  return ret;
};


// Distinguish between complex types and others
template <class S>
accelerator_inline Grid_simd1<S> operator*(Grid_simd1<S> a, Grid_simd1<S> b) {
  Grid_simd1<S> ret;
#ifndef STRICT_COMPLEX_MUL
  // Direct product, matching the vector types. std::complex adds an inf/NaN
  // recovery branch. Define STRICT_COMPLEX_MUL to restore std::complex.
  if constexpr ( is_complex<S>::value ) {
    ret.v = S(real(a.v)*real(b.v) - imag(a.v)*imag(b.v),
	      real(a.v)*imag(b.v) + imag(a.v)*real(b.v));
  } else {
    ret.v = a.v*b.v;
  }
#else
  ret.v = a.v*b.v;
#endif
  return ret;
};
///////////////////////
// Conjugate
///////////////////////
template <class S,  IfComplex<S> = 0>
accelerator_inline Grid_simd1<S> conjugate(const Grid_simd1<S> &in) {
  Grid_simd1<S> ret;
  ret.v = S(real(in.v),-imag(in.v));
  return ret;
}
template <class S,  IfNotComplex<S> = 0>
accelerator_inline Grid_simd1<S> conjugate(const Grid_simd1<S> &in) {
  return in;  // for real objects
}
// Suppress adj for integer types... // odd; why conjugate above but not adj??
template <class S,  IfNotInteger<S> = 0>
accelerator_inline Grid_simd1<S> adj(const Grid_simd1<S> &in) {
  return conjugate(in);
}

///////////////////////
// timesMinusI
///////////////////////
template <class S,  IfComplex<S> = 0>
accelerator_inline void timesMinusI(Grid_simd1<S> &ret, const Grid_simd1<S> &in) {
  ret.v = S(imag(in.v),-real(in.v));
}
template <class S,  IfComplex<S> = 0>
accelerator_inline Grid_simd1<S> timesMinusI(const Grid_simd1<S> &in) {
  Grid_simd1<S> ret;
  timesMinusI(ret,in);
  return ret;
}
template <class S,  IfNotComplex<S> = 0>
accelerator_inline Grid_simd1<S> timesMinusI(const Grid_simd1<S> &in) {
  return in;
}
///////////////////////
// timesI
///////////////////////
template <class S,  IfComplex<S> = 0>
accelerator_inline void timesI(Grid_simd1<S> &ret, const Grid_simd1<S> &in) {
  ret.v = S(-imag(in.v),real(in.v));
}
template <class S,  IfComplex<S> = 0>
accelerator_inline Grid_simd1<S> timesI(const Grid_simd1<S> &in) {
  Grid_simd1<S> ret;
  timesI(ret,in);
  return ret;
}
template <class S,  IfNotComplex<S> = 0>
accelerator_inline Grid_simd1<S> timesI(const Grid_simd1<S> &in) {
  return in;
}

// Distinguish between complex types and others
template <class S>
accelerator_inline Grid_simd1<S> operator/(Grid_simd1<S> a, Grid_simd1<S> b) {
  Grid_simd1<S> ret;
  ret.v = a.v/b.v;
  return ret;
};


/////////////////////
// Inner, outer
/////////////////////
template <class S>
accelerator_inline Grid_simd1<S> innerProduct(const Grid_simd1<S> &l,const Grid_simd1<S> &r) {
  return conjugate(l) * r;
}
template <class S>
accelerator_inline Grid_simd1<S> outerProduct(const Grid_simd1<S> &l,const Grid_simd1<S> &r) {
  return l * conjugate(r);
}

template <class S>
accelerator_inline Grid_simd1<S> trace(const Grid_simd1<S> &arg) {
  return arg;
}
////////////////////////////////////////////////////////////
// copy/splat complex real parts into real;
// insert real into complex and zero imag;
////////////////////////////////////////////////////////////
accelerator_inline sRealF toReal(const sComplexF &in) {
  sRealF ret;
  ret.v=real(in.v);
  return ret;
}
accelerator_inline sRealD toReal(const sComplexD &in) {
  sRealD ret;
  ret.v=real(in.v);
  return ret;
}
accelerator_inline sComplexF toComplex(sRealF &in)
{
  sComplexF ret;
  ret.v = in.v;
  return ret;
}
accelerator_inline sComplexD toComplex(sRealD &in)
{
  sComplexD ret;
  ret.v = in.v;
  return ret;
}

accelerator_inline void precisionChange(sRealF *out,const sRealD *in,int nvec){ assert(nvec==1); out->v = in->v;}
accelerator_inline void precisionChange(sRealD *out,const sRealF *in,int nvec){ assert(nvec==1); out->v = in->v;}
accelerator_inline void precisionChange(sComplexF *out,const sComplexD *in,int nvec){ assert(nvec==1); out->v = in->v;}
accelerator_inline void precisionChange(sComplexD *out,const sComplexF *in,int nvec){ assert(nvec==1); out->v = in->v;}


NAMESPACE_END(Grid);

