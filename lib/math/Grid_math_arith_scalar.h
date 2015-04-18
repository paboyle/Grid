#ifndef GRID_MATH_ARITH_SCALAR_H
#define GRID_MATH_ARITH_SCALAR_H

namespace Grid {


//////////////////////////////////////////////////////////////////////////////////////////
// Must support native C++ types Integer, Complex, Real
//////////////////////////////////////////////////////////////////////////////////////////

// multiplication by fundamental scalar type
template<class l,int N> inline iScalar<l> operator * (const iScalar<l>& lhs,const typename iScalar<l>::scalar_type rhs) 
{
  typename iScalar<l>::tensor_reduced srhs(rhs);
  return lhs*srhs;
}
template<class l,int N> inline iScalar<l> operator * (const typename iScalar<l>::scalar_type lhs,const iScalar<l>& rhs) {  return rhs*lhs; }

template<class l,int N> inline iVector<l,N> operator * (const iVector<l,N>& lhs,const typename iScalar<l>::scalar_type rhs) 
{
  typename iVector<l,N>::tensor_reduced srhs(rhs);
  return lhs*srhs;
}
template<class l,int N> inline iVector<l,N> operator * (const typename iScalar<l>::scalar_type lhs,const iVector<l,N>& rhs) {  return rhs*lhs; }

template<class l,int N> inline iMatrix<l,N> operator * (const iMatrix<l,N>& lhs,const typename iScalar<l>::scalar_type &rhs) 
{
  typename iMatrix<l,N>::tensor_reduced srhs(rhs);
  return lhs*srhs;
}
template<class l,int N> inline iMatrix<l,N> operator * (const typename iScalar<l>::scalar_type & lhs,const iMatrix<l,N>& rhs) {  return rhs*lhs; }

////////////////////////////////////////////////////////////////////
// Double support; cast to "scalar_type" through constructor
////////////////////////////////////////////////////////////////////
template<class l> inline iScalar<l> operator * (const iScalar<l>& lhs,double rhs) 
{
  typename iScalar<l>::scalar_type t(rhs);
  typename iScalar<l>::tensor_reduced srhs(t);
  return lhs*srhs;
}
template<class l> inline iScalar<l> operator * (double lhs,const iScalar<l>& rhs) {  return rhs*lhs; }

template<class l,int N> inline iVector<l,N> operator * (const iVector<l,N>& lhs,double rhs) 
{
  typename iScalar<l>::scalar_type t(rhs);
  typename iScalar<l>::tensor_reduced srhs(t);
  return lhs*srhs;
}
template<class l,int N> inline iVector<l,N> operator * (double lhs,const iVector<l,N>& rhs) {  return rhs*lhs; }

template<class l,int N> inline iMatrix<l,N> operator * (const iMatrix<l,N>& lhs,double rhs) 
{
  typename iScalar<l>::scalar_type t(rhs);
  typename iScalar<l>::tensor_reduced srhs(t);
  return lhs*srhs;
}
template<class l,int N> inline iMatrix<l,N> operator * (double lhs,const iMatrix<l,N>& rhs) {  return rhs*lhs; }

////////////////////////////////////////////////////////////////////
// Complex support; cast to "scalar_type" through constructor
////////////////////////////////////////////////////////////////////
template<class l> inline iScalar<l> operator * (const iScalar<l>& lhs,ComplexD rhs) 
{
  typename iScalar<l>::scalar_type t(rhs);
  typename iScalar<l>::tensor_reduced srhs(t);
  return lhs*srhs;
}
template<class l> inline iScalar<l> operator * (ComplexD lhs,const iScalar<l>& rhs) {  return rhs*lhs; }

template<class l,int N> inline iVector<l,N> operator * (const iVector<l,N>& lhs,ComplexD rhs) 
{
  typename iScalar<l>::scalar_type t(rhs);
  typename iScalar<l>::tensor_reduced srhs(t);
  return lhs*srhs;
}
template<class l,int N> inline iVector<l,N> operator * (ComplexD lhs,const iVector<l,N>& rhs) {  return rhs*lhs; }

template<class l,int N> inline iMatrix<l,N> operator * (const iMatrix<l,N>& lhs,ComplexD rhs) 
{
  typename iScalar<l>::scalar_type t(rhs);
  typename iScalar<l>::tensor_reduced srhs(t);
  return lhs*srhs;
}
template<class l,int N> inline iMatrix<l,N> operator * (ComplexD lhs,const iMatrix<l,N>& rhs) {  return rhs*lhs; }

////////////////////////////////////////////////////////////////////
// Integer support; cast to "scalar_type" through constructor
////////////////////////////////////////////////////////////////////
template<class l> inline iScalar<l> operator * (const iScalar<l>& lhs,Integer rhs) 
{
  typename iScalar<l>::scalar_type t(rhs);
  typename iScalar<l>::tensor_reduced srhs(t);
  return lhs*srhs;
}
template<class l> inline iScalar<l> operator * (Integer lhs,const iScalar<l>& rhs) {  return rhs*lhs; }

template<class l,int N> inline iVector<l,N> operator * (const iVector<l,N>& lhs,Integer rhs) 
{
  typename iScalar<l>::scalar_type t(rhs);
  typename iScalar<l>::tensor_reduced srhs(t);
  return lhs*srhs;
}
template<class l,int N> inline iVector<l,N> operator * (Integer lhs,const iVector<l,N>& rhs) {  return rhs*lhs; }

template<class l,int N> inline iMatrix<l,N> operator * (const iMatrix<l,N>& lhs,Integer rhs) 
{
  typename iScalar<l>::scalar_type t(rhs);
  typename iScalar<l>::tensor_reduced srhs(t);
  return lhs*srhs;
}
template<class l,int N> inline iMatrix<l,N> operator * (Integer lhs,const iMatrix<l,N>& rhs) {  return rhs*lhs; }



///////////////////////////////////////////////////////////////////////////////////////////////
// addition by fundamental scalar type applies to matrix(down diag) and scalar
///////////////////////////////////////////////////////////////////////////////////////////////
template<class l,int N> inline iScalar<l> operator + (const iScalar<l>& lhs,const typename iScalar<l>::scalar_type rhs) 
{
  typename iScalar<l>::tensor_reduced srhs(rhs);
  return lhs+srhs;
}
template<class l,int N> inline iScalar<l> operator + (const typename iScalar<l>::scalar_type lhs,const iScalar<l>& rhs) {  return rhs+lhs; }

template<class l,int N> inline iMatrix<l,N> operator + (const iMatrix<l,N>& lhs,const typename iScalar<l>::scalar_type rhs) 
{
  typename iMatrix<l,N>::tensor_reduced srhs(rhs);
  return lhs+srhs;
}
template<class l,int N> inline iMatrix<l,N> operator + (const typename iScalar<l>::scalar_type lhs,const iMatrix<l,N>& rhs) {  return rhs+lhs; }

////////////////////////////////////////////////////////////////////
// Double support; cast to "scalar_type" through constructor
////////////////////////////////////////////////////////////////////
template<class l> inline iScalar<l> operator + (const iScalar<l>& lhs,double rhs) 
{
  typename iScalar<l>::scalar_type t(rhs);
  typename iScalar<l>::tensor_reduced srhs(t);
  return lhs+srhs;
}
template<class l> inline iScalar<l> operator + (double lhs,const iScalar<l>& rhs) {  return rhs+lhs; }

template<class l,int N> inline iMatrix<l,N> operator + (const iMatrix<l,N>& lhs,double rhs) 
{
  typename iScalar<l>::scalar_type t(rhs);
  typename iScalar<l>::tensor_reduced srhs(t);
  return lhs+srhs;
}
template<class l,int N> inline iMatrix<l,N> operator + (double lhs,const iMatrix<l,N>& rhs) {  return rhs+lhs; }


// Integer support cast to scalar type through constructor


template<class l> inline iScalar<l> operator + (const iScalar<l>& lhs,Integer rhs) 
{
  typename iScalar<l>::scalar_type t(rhs);
  typename iScalar<l>::tensor_reduced srhs(t);
  return lhs+srhs;
}

template<class l> inline iScalar<l> operator + (Integer lhs,const iScalar<l>& rhs) {  return rhs+lhs; }

template<class l,int N> inline iMatrix<l,N> operator + (const iMatrix<l,N>& lhs,Integer rhs) 
{
  typename iScalar<l>::scalar_type t(rhs);
  typename iScalar<l>::tensor_reduced srhs(t);
  return lhs+srhs;
}
template<class l,int N> inline iMatrix<l,N> operator + (Integer lhs,const iMatrix<l,N>& rhs) {  return rhs+lhs; }


///////////////////////////////////////////////////////////////////////////////////////////////
// subtraction of fundamental scalar type applies to matrix(down diag) and scalar
///////////////////////////////////////////////////////////////////////////////////////////////
template<class l,int N> inline iScalar<l> operator - (const iScalar<l>& lhs,const typename iScalar<l>::scalar_type rhs) 
{
  typename iScalar<l>::tensor_reduced srhs(rhs);
  return lhs-srhs;
}
template<class l,int N> inline iScalar<l> operator - (const typename iScalar<l>::scalar_type lhs,const iScalar<l>& rhs) 
{
  typename iScalar<l>::tensor_reduced slhs(lhs);
  return slhs-rhs;
}

template<class l,int N> inline iMatrix<l,N> operator - (const iMatrix<l,N>& lhs,const typename iScalar<l>::scalar_type rhs) 
{
  typename iScalar<l>::tensor_reduced srhs(rhs);
  return lhs-srhs;
}
template<class l,int N> inline iMatrix<l,N> operator - (const typename iScalar<l>::scalar_type lhs,const iMatrix<l,N>& rhs) 
{
  typename iScalar<l>::tensor_reduced slhs(lhs);
  return slhs-rhs;
}

////////////////////////////////////////////////////////////////////
// Double support; cast to "scalar_type" through constructor
////////////////////////////////////////////////////////////////////
template<class l> inline iScalar<l> operator - (const iScalar<l>& lhs,double rhs) 
{
  typename iScalar<l>::scalar_type t(rhs);
  typename iScalar<l>::tensor_reduced srhs(t);
  return lhs-srhs;
}
template<class l> inline iScalar<l> operator - (double lhs,const iScalar<l>& rhs) 
{
  typename iScalar<l>::scalar_type t(lhs);
  typename iScalar<l>::tensor_reduced slhs(t);
  return slhs-rhs;
}

template<class l,int N> inline iMatrix<l,N> operator - (const iMatrix<l,N>& lhs,double rhs) 
{
  typename iScalar<l>::scalar_type t(rhs);
  typename iScalar<l>::tensor_reduced srhs(t);
  return lhs-srhs;
}
template<class l,int N> inline iMatrix<l,N> operator - (double lhs,const iMatrix<l,N>& rhs) 
{
  typename iScalar<l>::scalar_type t(lhs);
  typename iScalar<l>::tensor_reduced slhs(t);
  return slhs-rhs;
}

////////////////////////////////////////////////////////////////////
// Integer support; cast to "scalar_type" through constructor
////////////////////////////////////////////////////////////////////
template<class l> inline iScalar<l> operator - (const iScalar<l>& lhs,Integer rhs) 
{
  typename iScalar<l>::scalar_type t(rhs);
  typename iScalar<l>::tensor_reduced srhs(t);
  return lhs-srhs;
}
template<class l> inline iScalar<l> operator - (Integer lhs,const iScalar<l>& rhs) 
{
  typename iScalar<l>::scalar_type t(lhs);
  typename iScalar<l>::tensor_reduced slhs(t);
  return slhs-rhs;
}
template<class l,int N> inline iMatrix<l,N> operator - (const iMatrix<l,N>& lhs,Integer rhs) 
{
  typename iScalar<l>::scalar_type t(rhs);
  typename iScalar<l>::tensor_reduced srhs(t);
  return lhs-srhs;
}
template<class l,int N> inline iMatrix<l,N> operator - (Integer lhs,const iMatrix<l,N>& rhs) 
{
  typename iScalar<l>::scalar_type t(lhs);
  typename iScalar<l>::tensor_reduced slhs(t);
  return slhs-rhs;
}


}
#endif
