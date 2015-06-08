#ifndef GRID_VECTOR_UNOPS
#define GRID_VECTOR_UNOPS

namespace Grid { 

  template<class scalar> struct SqrtRealFunctor {
    scalar operator()(const scalar &a) const {
      return sqrt(real(a));
    }
  };

  template<class scalar> struct RSqrtRealFunctor {
    scalar operator()(const scalar &a)  const {
      return scalar(1.0/sqrt(real(a)));
    }
  };

  template<class scalar> struct CosRealFunctor {
    scalar operator()(const scalar &a)  const {
      return cos(real(a));
    }
  };

  template<class scalar> struct SinRealFunctor {
    scalar operator()(const scalar &a)  const {
      return sin(real(a));
    }
  };

  template<class scalar> struct PowRealFunctor {
    double y;
  PowRealFunctor(double _y) : y(_y) {};
    scalar operator()(const scalar &a)  const {
      return pow(real(a),y);
    }
  };

  template < class S, class V > 
  inline Grid_simd<S,V> sqrt(const Grid_simd<S,V> &r) {
    return SimdApply(SqrtRealFunctor<S>(),r);
  }
  template < class S, class V > 
  inline Grid_simd<S,V> rsqrt(const Grid_simd<S,V> &r) {
    return SimdApply(RSqrtRealFunctor<S>(),r);
  }
  template < class S, class V > 
  inline Grid_simd<S,V> cos(const Grid_simd<S,V> &r) {
    return SimdApply(CosRealFunctor<S>(),r);
  }
  template < class S, class V > 
  inline Grid_simd<S,V> sin(const Grid_simd<S,V> &r) {
    return SimdApply(CosRealFunctor<S>(),r);
  }
  template < class S, class V > 
  inline Grid_simd<S,V> pow(const Grid_simd<S,V> &r,double y) {
    return SimdApply(PowRealFunctor<S>(y),r);
  }

}
#endif
