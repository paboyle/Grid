#pragma once
NAMESPACE_BEGIN(Grid);
/////////////
// Unary operations
/////////////
template <class S>
accelerator_inline Grid_simd1<S> real(const Grid_simd1<S> &r) {
  return SimdApply(RealFunctor<S>(), r);
}
template <class S>
accelerator_inline Grid_simd1<S> imag(const Grid_simd1<S> &r) {
  return SimdApply(ImagFunctor<S>(), r);
}
template <class S>
accelerator_inline Grid_simd1<S> sqrt(const Grid_simd1<S> &r) {
  return SimdApply(SqrtRealFunctor<S>(), r);
}
template <class S>
accelerator_inline Grid_simd1<S> cos(const Grid_simd1<S> &r) {
  return SimdApply(CosRealFunctor<S>(), r);
}
template <class S>
accelerator_inline Grid_simd1<S> sin(const Grid_simd1<S> &r) {
  return SimdApply(SinRealFunctor<S>(), r);
}
template <class S>
accelerator_inline Grid_simd1<S> acos(const Grid_simd1<S> &r) {
  return SimdApply(AcosRealFunctor<S>(), r);
}
template <class S>
accelerator_inline Grid_simd1<S> asin(const Grid_simd1<S> &r) {
  return SimdApply(AsinRealFunctor<S>(), r);
}
template <class S>
accelerator_inline Grid_simd1<S> log(const Grid_simd1<S> &r) {
  return SimdApply(LogRealFunctor<S>(), r);
}
template <class S>
accelerator_inline Grid_simd1<S> abs(const Grid_simd1<S> &r) {
  return SimdApply(AbsRealFunctor<S>(), r);
}
template <class S>
accelerator_inline Grid_simd1<S> exp(const Grid_simd1<S> &r) {
  return SimdApply(ExpFunctor<S>(), r);
}
template <class S>
accelerator_inline Grid_simd1<S> Not(const Grid_simd1<S> &r) {
  return SimdApply(NotFunctor<S>(), r);
}
template <class S>
accelerator_inline Grid_simd1<S> pow(const Grid_simd1<S> &r, double y) {
  return SimdApply(PowRealFunctor<S>(y), r);
}
template <class S>
accelerator_inline Grid_simd1<S> mod(const Grid_simd1<S> &r, Integer y) {
  return SimdApply(ModIntFunctor<S>(y), r);
}
template <class S>
accelerator_inline Grid_simd1<S> div(const Grid_simd1<S> &r, Integer y) {
  return SimdApply(DivIntFunctor<S>(y), r);
}
NAMESPACE_END(Grid);
