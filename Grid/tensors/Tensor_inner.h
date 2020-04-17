/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/tensors/Tensor_inner.h

    Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: Christoph Lehner <christoph@lhnr.de>

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
#ifndef GRID_MATH_INNER_H
#define GRID_MATH_INNER_H

NAMESPACE_BEGIN(Grid);

///////////////////////////////////////////////////////////////////////////////////////
// innerProduct Scalar x Scalar -> Scalar
// innerProduct Vector x Vector -> Scalar
// innerProduct Matrix x Matrix -> Scalar
///////////////////////////////////////////////////////////////////////////////////////
template<class sobj> accelerator_inline RealD norm2(const sobj &arg){
  auto nrm = innerProductD(arg,arg);
  RealD ret = real(nrm);
  return ret;
}

////////////////////////////////////////////////////////////////
// FIXME horizontal sum: single promote to double and sum 
// Not yet used; make sum_cpu in LatticeReduction.h call these
// GPU sum_gpu is using a double precision promotion.
////////////////////////////////////////////////////////////////
#if 0
accelerator_inline ComplexD ReduceD(const ComplexF &l) { return l; };
accelerator_inline ComplexD ReduceD(const ComplexD &l) { return l; };
accelerator_inline RealD    ReduceD(const RealD    &l) { return l; };
accelerator_inline RealD    ReduceD(const RealF    &l) { return l; };

accelerator_inline ComplexD ReduceD(const vComplexD &l) { return Reduce(l); };
accelerator_inline RealD    ReduceD(const vRealD    &l) { return Reduce(l); };
accelerator_inline ComplexD ReduceD(const vComplexF  &l) 
{ 
  vComplexD la,lb;
  Optimization::PrecisionChange::StoD(l.v,la.v,lb.v);
  return Reduce(la)+Reduce(lb); 
}
accelerator_inline RealD    ReduceD(const vRealF    &l) 
{
  vRealD la,lb;
  Optimization::PrecisionChange::StoD(l.v,la.v,lb.v);
  return Reduce(la)+Reduce(lb); 
};
///////////////////////////////////////////////////////
// Now do it for vector, matrix, scalar
///////////////////////////////////////////////////////
template<class l,int N> accelerator_inline
auto ReduceD (const iVector<l,N>& lhs) -> iVector<decltype(ReduceD(lhs._internal[0])),N>
{
  typedef decltype(ReduceD(lhs._internal[0])) ret_t;
  iVector<ret_t,N> ret;
  for(int c1=0;c1<N;c1++){
    ret._internal[c1] = ReduceD(lhs._internal[c1]);
  }
  return ret;
}
template<class l,int N> accelerator_inline
auto ReduceD (const iMatrix<l,N>& lhs) -> iMatrix<decltype(ReduceD(lhs._internal[0][0])),N>
{
  typedef decltype(ReduceD(lhs._internal[0][0])) ret_t;
  iMatrix<ret_t,N> ret;
  for(int c1=0;c1<N;c1++){
  for(int c2=0;c2<N;c2++){
    ret._internal[c1][c2]=ReduceD(lhs._internal[c1][c2]);
  }}
  return ret;
}
template<class l> accelerator_inline
auto ReduceD (const iScalar<l>& lhs) -> iScalar<decltype(ReduceD(lhs._internal))>
{
  typedef decltype(ReduceD(lhs._internal)) ret_t;
  iScalar<ret_t> ret;
  ret._internal = ReduceD(lhs._internal);
  return ret;
}
#endif
///////////////////////////////////////////////////////
// Now do Reduce for vector, matrix, scalar
///////////////////////////////////////////////////////
template<class l,int N> accelerator_inline
auto Reduce (const iVector<l,N>& lhs) -> iVector<decltype(Reduce(lhs._internal[0])),N>
{
  typedef decltype(Reduce(lhs._internal[0])) ret_t;
  iVector<ret_t,N> ret;
  for(int c1=0;c1<N;c1++){
    ret._internal[c1] = Reduce(lhs._internal[c1]);
  }
  return ret;
}
template<class l,int N> accelerator_inline
auto Reduce (const iMatrix<l,N>& lhs) -> iMatrix<decltype(Reduce(lhs._internal[0][0])),N>
{
  typedef decltype(Reduce(lhs._internal[0][0])) ret_t;
  iMatrix<ret_t,N> ret;
  for(int c1=0;c1<N;c1++){
  for(int c2=0;c2<N;c2++){
    ret._internal[c1][c2]=Reduce(lhs._internal[c1][c2]);
  }}
  return ret;
}
template<class l> accelerator_inline
auto Reduce (const iScalar<l>& lhs) -> iScalar<decltype(Reduce(lhs._internal))>
{
  typedef decltype(Reduce(lhs._internal)) ret_t;
  iScalar<ret_t> ret;
  ret._internal = Reduce(lhs._internal);
  return ret;
}



//////////////////////////////////////
// innerProductD : if single promote to double and evaluate with sum 2x
//////////////////////////////////////
accelerator_inline ComplexD innerProductD(const ComplexF &l,const ComplexF &r){  return innerProduct(l,r); }
accelerator_inline ComplexD innerProductD(const ComplexD &l,const ComplexD &r){  return innerProduct(l,r); }
accelerator_inline RealD    innerProductD(const RealD    &l,const RealD    &r){  return innerProduct(l,r); }
accelerator_inline RealD    innerProductD(const RealF    &l,const RealF    &r){  return innerProduct(l,r); }

accelerator_inline vComplexD innerProductD(const vComplexD &l,const vComplexD &r){  return innerProduct(l,r); }
accelerator_inline vRealD    innerProductD(const vRealD    &l,const vRealD    &r){  return innerProduct(l,r); }
accelerator_inline vComplexD innerProductD(const vComplexF &l,const vComplexF &r)
{  
  vComplexD la,lb;
  vComplexD ra,rb;
  Optimization::PrecisionChange::StoD(l.v,la.v,lb.v);
  Optimization::PrecisionChange::StoD(r.v,ra.v,rb.v);
  return innerProduct(la,ra) + innerProduct(lb,rb); 
}
accelerator_inline vRealD innerProductD(const vRealF &l,const vRealF &r)
{  
  vRealD la,lb;
  vRealD ra,rb;
  Optimization::PrecisionChange::StoD(l.v,la.v,lb.v);
  Optimization::PrecisionChange::StoD(r.v,ra.v,rb.v);
  return innerProduct(la,ra) + innerProduct(lb,rb); 
}

// Now do it for vector, matrix, scalar
template<class l,class r,int N> accelerator_inline
auto innerProductD (const iVector<l,N>& lhs,const iVector<r,N>& rhs) -> iScalar<decltype(innerProductD(lhs._internal[0],rhs._internal[0]))>
{
  typedef decltype(innerProductD(lhs._internal[0],rhs._internal[0])) ret_t;
  iScalar<ret_t> ret;
  zeroit(ret);
  for(int c1=0;c1<N;c1++){
    ret._internal += innerProductD(lhs._internal[c1],rhs._internal[c1]);
  }
  return ret;
}
template<class l,class r,int N> accelerator_inline
auto innerProductD (const iMatrix<l,N>& lhs,const iMatrix<r,N>& rhs) -> iScalar<decltype(innerProductD(lhs._internal[0][0],rhs._internal[0][0]))>
{
  typedef decltype(innerProductD(lhs._internal[0][0],rhs._internal[0][0])) ret_t;
  iScalar<ret_t> ret;
  ret=Zero();
  for(int c1=0;c1<N;c1++){
    for(int c2=0;c2<N;c2++){
      ret._internal+=innerProductD(lhs._internal[c1][c2],rhs._internal[c1][c2]);
  }}
  return ret;
}
template<class l,class r> accelerator_inline
auto innerProductD (const iScalar<l>& lhs,const iScalar<r>& rhs) -> iScalar<decltype(innerProductD(lhs._internal,rhs._internal))>
{
  typedef decltype(innerProductD(lhs._internal,rhs._internal)) ret_t;
  iScalar<ret_t> ret;
  ret._internal = innerProductD(lhs._internal,rhs._internal);
  return ret;
}


//////////////////////////////////////
// innerProductD2: precision promotion without inner sum
//////////////////////////////////////
accelerator_inline vComplexD2 TensorRemove(const vComplexD2 & x) { return x; };
accelerator_inline vRealD2 TensorRemove(const vRealD2 & x) { return x; };

accelerator_inline ComplexD innerProductD2(const ComplexF &l,const ComplexF &r){  return innerProduct(l,r); }
accelerator_inline ComplexD innerProductD2(const ComplexD &l,const ComplexD &r){  return innerProduct(l,r); }
accelerator_inline RealD    innerProductD2(const RealD    &l,const RealD    &r){  return innerProduct(l,r); }
accelerator_inline RealD    innerProductD2(const RealF    &l,const RealF    &r){  return innerProduct(l,r); }

accelerator_inline vComplexD innerProductD2(const vComplexD &l,const vComplexD &r){  return innerProduct(l,r); }
accelerator_inline vRealD    innerProductD2(const vRealD    &l,const vRealD    &r){  return innerProduct(l,r); }

accelerator_inline vComplexD2 innerProductD2(const vComplexF &l,const vComplexF &r)
{  
  vComplexD la,lb;
  vComplexD ra,rb;
  Optimization::PrecisionChange::StoD(l.v,la.v,lb.v);
  Optimization::PrecisionChange::StoD(r.v,ra.v,rb.v);
  vComplexD2 ret;
  ret._internal[0] = innerProduct(la,ra);
  ret._internal[1] = innerProduct(lb,rb);
  return ret;
}
accelerator_inline vRealD2 innerProductD2(const vRealF &l,const vRealF &r)
{  
  vRealD la,lb;
  vRealD ra,rb;
  Optimization::PrecisionChange::StoD(l.v,la.v,lb.v);
  Optimization::PrecisionChange::StoD(r.v,ra.v,rb.v);
  vRealD2 ret;
  ret._internal[0]=innerProduct(la,ra);
  ret._internal[1]=innerProduct(lb,rb); 
  return ret;
}

// Now do it for vector, matrix, scalar
template<class l,class r,int N> accelerator_inline
  auto innerProductD2 (const iVector<l,N>& lhs,const iVector<r,N>& rhs) -> iScalar<decltype(innerProductD2(lhs._internal[0],rhs._internal[0]))>
{
  typedef decltype(innerProductD2(lhs._internal[0],rhs._internal[0])) ret_t;
  iScalar<ret_t> ret;
  zeroit(ret);
  for(int c1=0;c1<N;c1++){
    ret._internal += innerProductD2(lhs._internal[c1],rhs._internal[c1]);
  }
  return ret;
}
template<class l,class r,int N> accelerator_inline
  auto innerProductD2 (const iMatrix<l,N>& lhs,const iMatrix<r,N>& rhs) -> iScalar<decltype(innerProductD2(lhs._internal[0][0],rhs._internal[0][0]))>
{
  typedef decltype(innerProductD2(lhs._internal[0][0],rhs._internal[0][0])) ret_t;
  iScalar<ret_t> ret;
  ret=Zero();
  for(int c1=0;c1<N;c1++){
    for(int c2=0;c2<N;c2++){
      ret._internal+=innerProductD2(lhs._internal[c1][c2],rhs._internal[c1][c2]);
    }}
  return ret;
}
template<class l,class r> accelerator_inline
  auto innerProductD2 (const iScalar<l>& lhs,const iScalar<r>& rhs) -> iScalar<decltype(innerProductD2(lhs._internal,rhs._internal))>
{
  typedef decltype(innerProductD2(lhs._internal,rhs._internal)) ret_t;
  iScalar<ret_t> ret;
  ret._internal = innerProductD2(lhs._internal,rhs._internal);
  return ret;
}

//////////////////////
// Keep same precison
//////////////////////
template<class l,class r,int N> accelerator_inline
auto innerProduct (const iVector<l,N>& lhs,const iVector<r,N>& rhs) -> iScalar<decltype(innerProduct(lhs._internal[0],rhs._internal[0]))>
{
  typedef decltype(innerProduct(lhs._internal[0],rhs._internal[0])) ret_t;
  iScalar<ret_t> ret;
  ret=Zero();
  for(int c1=0;c1<N;c1++){
    ret._internal += innerProduct(lhs._internal[c1],rhs._internal[c1]);
  }
  return ret;
}
template<class l,class r,int N> accelerator_inline
auto innerProduct (const iMatrix<l,N>& lhs,const iMatrix<r,N>& rhs) -> iScalar<decltype(innerProduct(lhs._internal[0][0],rhs._internal[0][0]))>
{
  typedef decltype(innerProduct(lhs._internal[0][0],rhs._internal[0][0])) ret_t;
  iScalar<ret_t> ret;
  ret=Zero();
  for(int c1=0;c1<N;c1++){
    for(int c2=0;c2<N;c2++){
      ret._internal+=innerProduct(lhs._internal[c1][c2],rhs._internal[c1][c2]);
    }}
  return ret;
}
template<class l,class r> accelerator_inline
auto innerProduct (const iScalar<l>& lhs,const iScalar<r>& rhs) -> iScalar<decltype(innerProduct(lhs._internal,rhs._internal))>
{
  typedef decltype(innerProduct(lhs._internal,rhs._internal)) ret_t;
  iScalar<ret_t> ret;
  ret._internal = innerProduct(lhs._internal,rhs._internal);
  return ret;
}

NAMESPACE_END(Grid);

#endif
