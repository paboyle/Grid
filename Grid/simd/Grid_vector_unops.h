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

See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
			   /*  END LEGAL */
#ifndef GRID_VECTOR_UNOPS
#define GRID_VECTOR_UNOPS

#include <cmath>

NAMESPACE_BEGIN(Grid);

template <class scalar>
struct SqrtRealFunctor {
  accelerator scalar operator()(const scalar &a) const { return sqrt(real(a)); }
};

template <class scalar>
struct RSqrtRealFunctor {
  accelerator scalar operator()(const scalar &a) const {
    return scalar(1.0 / sqrt(real(a)));
  }
};

template <class scalar>
struct CosRealFunctor {
  accelerator scalar operator()(const scalar &a) const { return cos(real(a)); }
};

template <class scalar>
struct SinRealFunctor {
  accelerator scalar operator()(const scalar &a) const { return sin(real(a)); }
};

template <class scalar>
struct AcosRealFunctor {
  accelerator scalar operator()(const scalar &a) const { return acos(real(a)); }
};

template <class scalar>
struct AsinRealFunctor {
  accelerator scalar operator()(const scalar &a) const { return asin(real(a)); }
};
template <class scalar>
struct LogRealFunctor {
  accelerator scalar operator()(const scalar &a) const { return log(real(a)); }
};
template <class scalar>
struct ExpFunctor {
  accelerator scalar operator()(const scalar &a) const { return exp(a); }
};
template <class scalar>
struct NotFunctor {
  accelerator scalar operator()(const scalar &a) const { return (!a); }
};
template <class scalar>
struct AbsRealFunctor {
  accelerator scalar operator()(const scalar &a) const { return std::abs(real(a)); }
};
template <class scalar>
struct PowRealFunctor {
  double y;
  accelerator PowRealFunctor(double _y) : y(_y){};
  accelerator scalar operator()(const scalar &a) const { return pow(real(a), y); }
};

template <class scalar>
struct ModIntFunctor {
  Integer y;
  accelerator ModIntFunctor(Integer _y) : y(_y){};
  accelerator scalar operator()(const scalar &a) const { return Integer(a) % y; }
};

template <class scalar>
struct DivIntFunctor {
  Integer y;
  accelerator DivIntFunctor(Integer _y) : y(_y){};
  accelerator scalar operator()(const scalar &a) const { return Integer(a) / y; }
};

template <class scalar>
struct RealFunctor {
  accelerator scalar operator()(const scalar &a) const { return real(a); }
};
template <class scalar>
struct ImagFunctor {
  accelerator scalar operator()(const scalar &a) const { return imag(a); }
};
template <class S, class V>
accelerator_inline Grid_simd<S, V> real(const Grid_simd<S, V> &r) {
  return SimdApply(RealFunctor<S>(), r);
}
template <class S, class V>
accelerator_inline Grid_simd<S, V> imag(const Grid_simd<S, V> &r) {
  return SimdApply(ImagFunctor<S>(), r);
}
template <class S, class V>
accelerator_inline Grid_simd<S, V> sqrt(const Grid_simd<S, V> &r) {
  return SimdApply(SqrtRealFunctor<S>(), r);
}
template <class S, class V>
accelerator_inline Grid_simd<S, V> cos(const Grid_simd<S, V> &r) {
  return SimdApply(CosRealFunctor<S>(), r);
}
template <class S, class V>
accelerator_inline Grid_simd<S, V> sin(const Grid_simd<S, V> &r) {
  return SimdApply(SinRealFunctor<S>(), r);
}
template <class S, class V>
accelerator_inline Grid_simd<S, V> acos(const Grid_simd<S, V> &r) {
  return SimdApply(AcosRealFunctor<S>(), r);
}
template <class S, class V>
accelerator_inline Grid_simd<S, V> asin(const Grid_simd<S, V> &r) {
  return SimdApply(AsinRealFunctor<S>(), r);
}
template <class S, class V>
accelerator_inline Grid_simd<S, V> log(const Grid_simd<S, V> &r) {
  return SimdApply(LogRealFunctor<S>(), r);
}
template <class S, class V>
accelerator_inline Grid_simd<S, V> abs(const Grid_simd<S, V> &r) {
  return SimdApply(AbsRealFunctor<S>(), r);
}
template <class S, class V>
accelerator_inline Grid_simd<S, V> exp(const Grid_simd<S, V> &r) {
  return SimdApply(ExpFunctor<S>(), r);
}
template <class S, class V>
accelerator_inline Grid_simd<S, V> Not(const Grid_simd<S, V> &r) {
  return SimdApply(NotFunctor<S>(), r);
}
template <class S, class V>
accelerator_inline Grid_simd<S, V> pow(const Grid_simd<S, V> &r, double y) {
  return SimdApply(PowRealFunctor<S>(y), r);
}
template <class S, class V>
accelerator_inline Grid_simd<S, V> mod(const Grid_simd<S, V> &r, Integer y) {
  return SimdApply(ModIntFunctor<S>(y), r);
}
template <class S, class V>
accelerator_inline Grid_simd<S, V> div(const Grid_simd<S, V> &r, Integer y) {
  return SimdApply(DivIntFunctor<S>(y), r);
}
////////////////////////////////////////////////////////////////////////////
// Allows us to assign into **conformable** real vectors from complex
////////////////////////////////////////////////////////////////////////////
template <class scalar>
struct AndFunctor {
  accelerator scalar operator()(const scalar &x, const scalar &y) const { return x & y; }
};
template <class scalar>
struct OrFunctor {
  accelerator scalar operator()(const scalar &x, const scalar &y) const { return x | y; }
};
template <class scalar>
struct AndAndFunctor {
  accelerator scalar operator()(const scalar &x, const scalar &y) const { return x && y; }
};
template <class scalar>
struct OrOrFunctor {
  accelerator scalar operator()(const scalar &x, const scalar &y) const { return x || y; }
};

////////////////////////////////
// Calls to simd binop functors
////////////////////////////////
template <class S, class V>
accelerator_inline Grid_simd<S, V> operator&(const Grid_simd<S, V> &x,
				 const Grid_simd<S, V> &y) {
  return SimdApplyBinop(AndFunctor<S>(), x, y);
}
template <class S, class V>
accelerator_inline Grid_simd<S, V> operator&&(const Grid_simd<S, V> &x,
				  const Grid_simd<S, V> &y) {
  return SimdApplyBinop(AndAndFunctor<S>(), x, y);
}
template <class S, class V>
accelerator_inline Grid_simd<S, V> operator|(const Grid_simd<S, V> &x,
				 const Grid_simd<S, V> &y) {
  return SimdApplyBinop(OrFunctor<S>(), x, y);
}
template <class S, class V>
accelerator_inline Grid_simd<S, V> operator||(const Grid_simd<S, V> &x,
				  const Grid_simd<S, V> &y) {
  return SimdApplyBinop(OrOrFunctor<S>(), x, y);
}
NAMESPACE_END(Grid);
#endif
