/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/lattice/Lattice_comparison.h

    Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
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
#ifndef GRID_LATTICE_COMPARISON_H
#define GRID_LATTICE_COMPARISON_H

NAMESPACE_BEGIN(Grid);

//////////////////////////////////////////////////////////////////////////
// relational operators
// 
// Support <,>,<=,>=,==,!=
//
//Query supporting bitwise &, |, ^, !
//Query supporting logical &&, ||, 
//////////////////////////////////////////////////////////////////////////

typedef iScalar<vInteger> vPredicate ;

/*
template <class iobj, class vobj, class robj> accelerator_inline 
vobj predicatedWhere(const iobj &predicate, const vobj &iftrue, const robj &iffalse) 
{
  typename std::remove_const<vobj>::type ret;

  typedef typename vobj::scalar_object scalar_object;
  typedef typename vobj::scalar_type scalar_type;
  typedef typename vobj::vector_type vector_type;

  const int Nsimd = vobj::vector_type::Nsimd();

  ExtractBuffer<Integer> mask(Nsimd);
  ExtractBuffer<scalar_object> truevals(Nsimd);
  ExtractBuffer<scalar_object> falsevals(Nsimd);

  extract(iftrue, truevals);
  extract(iffalse, falsevals);
  extract<vInteger, Integer>(TensorRemove(predicate), mask);

  for (int s = 0; s < Nsimd; s++) {
    if (mask[s]) falsevals[s] = truevals[s];
  }

  merge(ret, falsevals);
  return ret;
}
*/
//////////////////////////////////////////////////////////////////////////
// compare lattice to lattice
//////////////////////////////////////////////////////////////////////////

template<class vfunctor,class lobj,class robj>  
inline Lattice<vPredicate> LLComparison(vfunctor op,const Lattice<lobj> &lhs,const Lattice<robj> &rhs)
{
  Lattice<vPredicate> ret(rhs.Grid());
  auto lhs_v = lhs.View();
  auto rhs_v = rhs.View();
  auto ret_v = ret.View();
  thread_for( ss, rhs_v.size(), {
      ret_v[ss]=op(lhs_v[ss],rhs_v[ss]);
  });
  return ret;
}
//////////////////////////////////////////////////////////////////////////
// compare lattice to scalar
//////////////////////////////////////////////////////////////////////////
template<class vfunctor,class lobj,class robj> 
inline Lattice<vPredicate> LSComparison(vfunctor op,const Lattice<lobj> &lhs,const robj &rhs)
{
  Lattice<vPredicate> ret(lhs.Grid());
  auto lhs_v = lhs.View();
  auto ret_v = ret.View();
  thread_for( ss, lhs_v.size(), {
    ret_v[ss]=op(lhs_v[ss],rhs);
  });
  return ret;
}
//////////////////////////////////////////////////////////////////////////
// compare scalar to lattice
//////////////////////////////////////////////////////////////////////////
template<class vfunctor,class lobj,class robj> 
inline Lattice<vPredicate> SLComparison(vfunctor op,const lobj &lhs,const Lattice<robj> &rhs)
{
  Lattice<vPredicate> ret(rhs.Grid());
  auto rhs_v = rhs.View();
  auto ret_v = ret.View();
  thread_for( ss, rhs_v.size(), {
    ret_v[ss]=op(lhs,rhs_v[ss]);
  });
  return ret;
}
  
//////////////////////////////////////////////////////////////////////////
// Map to functors
//////////////////////////////////////////////////////////////////////////
// Less than
template<class lobj,class robj>
inline Lattice<vPredicate> operator < (const Lattice<lobj> & lhs, const Lattice<robj> & rhs) {
  return LLComparison(vlt<lobj,robj>(),lhs,rhs);
}
template<class lobj,class robj>
inline Lattice<vPredicate> operator < (const Lattice<lobj> & lhs, const robj & rhs) {
  return LSComparison(vlt<lobj,robj>(),lhs,rhs);
}
template<class lobj,class robj>
inline Lattice<vPredicate> operator < (const lobj & lhs, const Lattice<robj> & rhs) {
  return SLComparison(vlt<lobj,robj>(),lhs,rhs);
}
  
// Less than equal
template<class lobj,class robj>
inline Lattice<vPredicate> operator <= (const Lattice<lobj> & lhs, const Lattice<robj> & rhs) {
  return LLComparison(vle<lobj,robj>(),lhs,rhs);
}
template<class lobj,class robj>
inline Lattice<vPredicate> operator <= (const Lattice<lobj> & lhs, const robj & rhs) {
  return LSComparison(vle<lobj,robj>(),lhs,rhs);
}
template<class lobj,class robj>
inline Lattice<vPredicate> operator <= (const lobj & lhs, const Lattice<robj> & rhs) {
  return SLComparison(vle<lobj,robj>(),lhs,rhs);
}
  
// Greater than 
template<class lobj,class robj>
inline Lattice<vPredicate> operator > (const Lattice<lobj> & lhs, const Lattice<robj> & rhs) {
  return LLComparison(vgt<lobj,robj>(),lhs,rhs);
}
template<class lobj,class robj>
inline Lattice<vPredicate> operator > (const Lattice<lobj> & lhs, const robj & rhs) {
  return LSComparison(vgt<lobj,robj>(),lhs,rhs);
}
template<class lobj,class robj>
inline Lattice<vPredicate> operator > (const lobj & lhs, const Lattice<robj> & rhs) {
  return SLComparison(vgt<lobj,robj>(),lhs,rhs);
}
  
  
// Greater than equal
template<class lobj,class robj>
inline Lattice<vPredicate> operator >= (const Lattice<lobj> & lhs, const Lattice<robj> & rhs) {
  return LLComparison(vge<lobj,robj>(),lhs,rhs);
}
template<class lobj,class robj>
inline Lattice<vPredicate> operator >= (const Lattice<lobj> & lhs, const robj & rhs) {
  return LSComparison(vge<lobj,robj>(),lhs,rhs);
}
template<class lobj,class robj>
inline Lattice<vPredicate> operator >= (const lobj & lhs, const Lattice<robj> & rhs) {
  return SLComparison(vge<lobj,robj>(),lhs,rhs);
}
   
// equal
template<class lobj,class robj>
inline Lattice<vPredicate> operator == (const Lattice<lobj> & lhs, const Lattice<robj> & rhs) {
  return LLComparison(veq<lobj,robj>(),lhs,rhs);
}
template<class lobj,class robj>
inline Lattice<vPredicate> operator == (const Lattice<lobj> & lhs, const robj & rhs) {
  return LSComparison(veq<lobj,robj>(),lhs,rhs);
}
template<class lobj,class robj>
inline Lattice<vPredicate> operator == (const lobj & lhs, const Lattice<robj> & rhs) {
  return SLComparison(veq<lobj,robj>(),lhs,rhs);
}
   
   
// not equal
template<class lobj,class robj>
inline Lattice<vPredicate> operator != (const Lattice<lobj> & lhs, const Lattice<robj> & rhs) {
  return LLComparison(vne<lobj,robj>(),lhs,rhs);
}
template<class lobj,class robj>
inline Lattice<vPredicate> operator != (const Lattice<lobj> & lhs, const robj & rhs) {
  return LSComparison(vne<lobj,robj>(),lhs,rhs);
}
template<class lobj,class robj>
inline Lattice<vPredicate> operator != (const lobj & lhs, const Lattice<robj> & rhs) {
  return SLComparison(vne<lobj,robj>(),lhs,rhs);
}
NAMESPACE_END(Grid);
#endif
