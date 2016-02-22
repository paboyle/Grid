    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/simd/Grid_vector_unops.h

    Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
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
#ifndef GRID_VECTOR_UNOPS
#define GRID_VECTOR_UNOPS

#include <cmath>

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

  template<class scalar> struct LogRealFunctor {
    scalar operator()(const scalar &a)  const {
      return log(real(a));
    }
  };

  template<class scalar> struct ExpRealFunctor {
    scalar operator()(const scalar &a)  const {
      return exp(real(a));
    }
  };
  template<class scalar> struct NotFunctor {
    scalar operator()(const scalar &a)  const {
      return (!a);
    }
  };
  template<class scalar> struct AbsRealFunctor {
    scalar operator()(const scalar &a)  const {
      return std::abs(real(a));
    }
  };

  template<class scalar> struct PowRealFunctor {
    double y;
  PowRealFunctor(double _y) : y(_y) {};
    scalar operator()(const scalar &a)  const {
      return pow(real(a),y);
    }
  };

  template<class scalar> struct ModIntFunctor {
    Integer y;
  ModIntFunctor(Integer _y) : y(_y) {};
    scalar operator()(const scalar &a)  const {
      return Integer(a)%y;
    }
  };

  template<class scalar> struct DivIntFunctor {
    Integer y;
  DivIntFunctor(Integer _y) : y(_y) {};
    scalar operator()(const scalar &a)  const {
      return Integer(a)/y;
    }
  };

  template<class scalar> struct RealFunctor {
    scalar operator()(const scalar &a)  const {
      return real(a);
    }
  };
  template<class scalar> struct ImagFunctor {
    scalar operator()(const scalar &a)  const {
      return imag(a);
    }
  };
  template < class S, class V > 
  inline Grid_simd<S,V> real(const Grid_simd<S,V> &r) {
    return SimdApply(RealFunctor<S>(),r);
  }
  template < class S, class V > 
  inline Grid_simd<S,V> imag(const Grid_simd<S,V> &r) {
    return SimdApply(ImagFunctor<S>(),r);
  }

  template < class S, class V > 
  inline Grid_simd<S,V> sqrt(const Grid_simd<S,V> &r) {
    return SimdApply(SqrtRealFunctor<S>(),r);
  }
  template < class S, class V > 
  inline Grid_simd<S,V> rsqrt(const Grid_simd<S,V> &r) {
    return SimdApply(RSqrtRealFunctor<S>(),r);
  }
  template < class Scalar > 
  inline Scalar rsqrt(const Scalar &r) {
    return (RSqrtRealFunctor<Scalar>(),r);
  }

  template < class S, class V > 
  inline Grid_simd<S,V> cos(const Grid_simd<S,V> &r) {
    return SimdApply(CosRealFunctor<S>(),r);
  }
  template < class S, class V > 
  inline Grid_simd<S,V> sin(const Grid_simd<S,V> &r) {
    return SimdApply(SinRealFunctor<S>(),r);
  }
  template < class S, class V > 
  inline Grid_simd<S,V> log(const Grid_simd<S,V> &r) {
    return SimdApply(LogRealFunctor<S>(),r);
  }
  template < class S, class V > 
  inline Grid_simd<S,V> abs(const Grid_simd<S,V> &r) {
    return SimdApply(AbsRealFunctor<S>(),r);
  }
  template < class S, class V > 
  inline Grid_simd<S,V> exp(const Grid_simd<S,V> &r) {
    return SimdApply(ExpRealFunctor<S>(),r);
  }
  template < class S, class V > 
  inline Grid_simd<S,V> Not(const Grid_simd<S,V> &r) {
    return SimdApply(NotFunctor<S>(),r);
  }
  template < class S, class V > 
  inline Grid_simd<S,V> pow(const Grid_simd<S,V> &r,double y) {
    return SimdApply(PowRealFunctor<S>(y),r);
  }
  template < class S, class V > 
  inline Grid_simd<S,V> mod(const Grid_simd<S,V> &r,Integer y) {
    return SimdApply(ModIntFunctor<S>(y),r);
  }
  template < class S, class V > 
  inline Grid_simd<S,V> div(const Grid_simd<S,V> &r,Integer y) {
    return SimdApply(DivIntFunctor<S>(y),r);
  }
  ////////////////////////////////////////////////////////////////////////////
  // Allows us to assign into **conformable** real vectors from complex
  ////////////////////////////////////////////////////////////////////////////
  //  template < class S, class V > 
  //  inline auto ComplexRemove(const Grid_simd<S,V> &c) -> Grid_simd<Grid_simd<S,V>::Real,V> {
  //    Grid_simd<Grid_simd<S,V>::Real,V> ret;
  //    ret.v = c.v;
  //    return ret;
  //  }
  template<class scalar> struct AndFunctor {
    scalar operator()(const scalar &x, const scalar &y)  const {
      return x & y;
    }
  };
  template<class scalar> struct OrFunctor {
    scalar operator()(const scalar &x, const scalar &y)  const {
      return x | y;
    }
  };
  template<class scalar> struct AndAndFunctor {
    scalar operator()(const scalar &x, const scalar &y)  const {
      return x && y;
    }
  };
  template<class scalar> struct OrOrFunctor {
    scalar operator()(const scalar &x, const scalar &y)  const {
      return x || y;
    }
  };

  ////////////////////////////////
  // Calls to simd binop functors
  ////////////////////////////////
  template < class S, class V > 
  inline Grid_simd<S,V> operator &(const Grid_simd<S,V> &x,const Grid_simd<S,V> &y) {
    return SimdApplyBinop(AndFunctor<S>(),x,y);
  }
  template < class S, class V > 
  inline Grid_simd<S,V> operator &&(const Grid_simd<S,V> &x,const Grid_simd<S,V> &y) {
    return SimdApplyBinop(AndAndFunctor<S>(),x,y);
  }
  template < class S, class V > 
  inline Grid_simd<S,V> operator |(const Grid_simd<S,V> &x,const Grid_simd<S,V> &y) {
    return SimdApplyBinop(OrFunctor<S>(),x,y);
  }
  template < class S, class V > 
  inline Grid_simd<S,V> operator ||(const Grid_simd<S,V> &x,const Grid_simd<S,V> &y) {
    return SimdApplyBinop(OrOrFunctor<S>(),x,y);
  }

}
#endif
