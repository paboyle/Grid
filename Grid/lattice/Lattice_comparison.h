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

template<class vobj> using vPredicate = iScalar<IntegerPredicate<vobj> > ;

//////////////////////////////////////////////////////////////////////////
// compare lattice to lattice
//////////////////////////////////////////////////////////////////////////

template<class vfunctor,class lobj,class robj>  
inline Lattice<vPredicate<lobj> > LLComparison(vfunctor op,const Lattice<lobj> &lhs,const Lattice<robj> &rhs)
{
  Lattice<vPredicate<lobj> > ret(rhs.Grid());
  autoView( lhs_v, lhs, CpuRead);
  autoView( rhs_v, rhs, CpuRead);
  autoView( ret_v, ret, CpuWrite);
  thread_for( ss, rhs_v.size(), {
      ret_v[ss]=op(lhs_v[ss],rhs_v[ss]);
  });
  return ret;
}
//////////////////////////////////////////////////////////////////////////
// compare lattice to scalar
//////////////////////////////////////////////////////////////////////////
template<class vfunctor,class lobj,class robj> 
inline Lattice<vPredicate<lobj> > LSComparison(vfunctor op,const Lattice<lobj> &lhs,const robj &rhs)
{
  Lattice<vPredicate<lobj> > ret(lhs.Grid());
  autoView( lhs_v, lhs, CpuRead);
  autoView( ret_v, ret, CpuWrite);
  thread_for( ss, lhs_v.size(), {
    ret_v[ss]=op(lhs_v[ss],rhs);
  });
  return ret;
}
//////////////////////////////////////////////////////////////////////////
// compare scalar to lattice
//////////////////////////////////////////////////////////////////////////
template<class vfunctor,class lobj,class robj> 
inline Lattice<vPredicate<robj> > SLComparison(vfunctor op,const lobj &lhs,const Lattice<robj> &rhs)
{
  Lattice<vPredicate<robj> > ret(rhs.Grid());
  autoView( rhs_v, rhs, CpuRead);
  autoView( ret_v, ret, CpuWrite);
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
inline Lattice<vPredicate<lobj> > operator < (const Lattice<lobj> & lhs, const Lattice<robj> & rhs) {
  return LLComparison(vlt<lobj,robj>(),lhs,rhs);
}
template<class lobj,class robj>
inline Lattice<vPredicate<lobj> > operator < (const Lattice<lobj> & lhs, const robj & rhs) {
  return LSComparison(vlt<lobj,robj>(),lhs,rhs);
}
template<class lobj,class robj>
inline Lattice<vPredicate<robj> > operator < (const lobj & lhs, const Lattice<robj> & rhs) {
  return SLComparison(vlt<lobj,robj>(),lhs,rhs);
}
  
// Less than equal
template<class lobj,class robj>
inline Lattice<vPredicate<lobj> > operator <= (const Lattice<lobj> & lhs, const Lattice<robj> & rhs) {
  return LLComparison(vle<lobj,robj>(),lhs,rhs);
}
template<class lobj,class robj>
inline Lattice<vPredicate<lobj> > operator <= (const Lattice<lobj> & lhs, const robj & rhs) {
  return LSComparison(vle<lobj,robj>(),lhs,rhs);
}
template<class lobj,class robj>
inline Lattice<vPredicate<robj> > operator <= (const lobj & lhs, const Lattice<robj> & rhs) {
  return SLComparison(vle<lobj,robj>(),lhs,rhs);
}
  
// Greater than 
template<class lobj,class robj>
inline Lattice<vPredicate<lobj> > operator > (const Lattice<lobj> & lhs, const Lattice<robj> & rhs) {
  return LLComparison(vgt<lobj,robj>(),lhs,rhs);
}
template<class lobj,class robj>
inline Lattice<vPredicate<lobj> > operator > (const Lattice<lobj> & lhs, const robj & rhs) {
  return LSComparison(vgt<lobj,robj>(),lhs,rhs);
}
template<class lobj,class robj>
inline Lattice<vPredicate<robj> > operator > (const lobj & lhs, const Lattice<robj> & rhs) {
  return SLComparison(vgt<lobj,robj>(),lhs,rhs);
}
  
  
// Greater than equal
template<class lobj,class robj>
inline Lattice<vPredicate<lobj> > operator >= (const Lattice<lobj> & lhs, const Lattice<robj> & rhs) {
  return LLComparison(vge<lobj,robj>(),lhs,rhs);
}
template<class lobj,class robj>
inline Lattice<vPredicate<lobj> > operator >= (const Lattice<lobj> & lhs, const robj & rhs) {
  return LSComparison(vge<lobj,robj>(),lhs,rhs);
}
template<class lobj,class robj>
inline Lattice<vPredicate<robj> > operator >= (const lobj & lhs, const Lattice<robj> & rhs) {
  return SLComparison(vge<lobj,robj>(),lhs,rhs);
}
   
// equal
template<class lobj,class robj>
inline Lattice<vPredicate<lobj> > operator == (const Lattice<lobj> & lhs, const Lattice<robj> & rhs) {
  return LLComparison(veq<lobj,robj>(),lhs,rhs);
}
template<class lobj,class robj>
inline Lattice<vPredicate<lobj> > operator == (const Lattice<lobj> & lhs, const robj & rhs) {
  return LSComparison(veq<lobj,robj>(),lhs,rhs);
}
template<class lobj,class robj>
inline Lattice<vPredicate<robj> > operator == (const lobj & lhs, const Lattice<robj> & rhs) {
  return SLComparison(veq<lobj,robj>(),lhs,rhs);
}
   
   
// not equal
template<class lobj,class robj>
inline Lattice<vPredicate<lobj> > operator != (const Lattice<lobj> & lhs, const Lattice<robj> & rhs) {
  return LLComparison(vne<lobj,robj>(),lhs,rhs);
}
template<class lobj,class robj>
inline Lattice<vPredicate<lobj> > operator != (const Lattice<lobj> & lhs, const robj & rhs) {
  return LSComparison(vne<lobj,robj>(),lhs,rhs);
}
template<class lobj,class robj>
inline Lattice<vPredicate<robj> > operator != (const lobj & lhs, const Lattice<robj> & rhs) {
  return SLComparison(vne<lobj,robj>(),lhs,rhs);
}

//////////////////////////////////////////////////////////////////////////
// Real-part relational comparison for COMPLEX lattices.
//
// Complex has no ordering, so operator<,>,<=,>= are (deliberately) undefined for
// complex operands -- Comparison() in Lattice_comparison_utils.h is IfNotComplex --
// and where() cannot be driven by a complex lattice directly.  It is however often
// useful to threshold on the REAL PART of a complex (scalar/singlet) field, e.g. a
// momentum-magnitude mask phat^2 > pc^2.  These free functions compare the real part
// to a real threshold and return the matching IntegerPredicate, so the result feeds
// where() exactly like the built-in relationals.
//
// Written with the per-lane getlane/putlane accessors (as in FFT.h, PaddedCell.h),
// not extract/merge buffers.  One wrinkle: the predicate type IntegerPredicate<CComplex>
// is vInteger, whose Nsimd (>= the widest real type's) exceeds the complex operand's
// Nsimd by s = Npred/Nsimd.  where() only reads the ii=0 representative of each group,
// at physical lane lane*s (the "s-fold skip" -- cf. extract()'s getlane(i*s)), so
// filling the rest of the group is not functionally required; we replicate the value
// across all s lanes anyway because some Grid code asserts the s replicas are equal
// (and it mirrors what merge() does).  When vInteger is reworked to carry the operand
// Nsimd (scope later) s becomes 1 and the inner loop drops out.
//////////////////////////////////////////////////////////////////////////
template<class scalar> class sRealLt { public:
  accelerator_inline Integer operator()(const scalar &a, RealD b) const { return a.real() <  b ? 1 : 0; } };
template<class scalar> class sRealLe { public:
  accelerator_inline Integer operator()(const scalar &a, RealD b) const { return a.real() <= b ? 1 : 0; } };
template<class scalar> class sRealGt { public:
  accelerator_inline Integer operator()(const scalar &a, RealD b) const { return a.real() >  b ? 1 : 0; } };
template<class scalar> class sRealGe { public:
  accelerator_inline Integer operator()(const scalar &a, RealD b) const { return a.real() >= b ? 1 : 0; } };

template<class sfunctor,class CComplex>
inline Lattice<vPredicate<CComplex> > RealPartComparison(sfunctor op,const Lattice<CComplex> &lhs, RealD thr)
{
  Lattice<vPredicate<CComplex> > ret(lhs.Grid());
  autoView( lv, lhs, AcceleratorRead);
  autoView( rv, ret, AcceleratorWrite);
  typedef typename CComplex::vector_type vsimd;
  const int Nsimd = vsimd::Nsimd();
  const int s     = IntegerPredicate<CComplex>::Nsimd() / Nsimd;   // lane-replication factor (see note)
  accelerator_for(ss, lhs.Grid()->oSites(), Nsimd, {
    vsimd v = TensorRemove(lv[ss]);                 // strip iScalar nest to the bare complex SIMD word
#ifdef GRID_SIMT
    { int lane = acceleratorSIMTlane(Nsimd);        // GPU: this thread == this operand lane
#else
    for(int lane=0;lane<Nsimd;lane++){              // CPU: walk the packed operand lanes
#endif
      Integer p = op(v.getlane(lane), thr);
      for(int ii=0;ii<s;ii++) rv[ss]._internal.putlane(p, lane*s+ii);   // replicate across the group
#ifdef GRID_SIMT
    }
#else
    }
#endif
  });
  return ret;
}

template<class CComplex> inline Lattice<vPredicate<CComplex> > RealPartLessThan    (const Lattice<CComplex> &a, RealD b){ return RealPartComparison(sRealLt<typename CComplex::vector_type::scalar_type>(),a,b); }
template<class CComplex> inline Lattice<vPredicate<CComplex> > RealPartLessEqual    (const Lattice<CComplex> &a, RealD b){ return RealPartComparison(sRealLe<typename CComplex::vector_type::scalar_type>(),a,b); }
template<class CComplex> inline Lattice<vPredicate<CComplex> > RealPartGreaterThan  (const Lattice<CComplex> &a, RealD b){ return RealPartComparison(sRealGt<typename CComplex::vector_type::scalar_type>(),a,b); }
template<class CComplex> inline Lattice<vPredicate<CComplex> > RealPartGreaterEqual (const Lattice<CComplex> &a, RealD b){ return RealPartComparison(sRealGe<typename CComplex::vector_type::scalar_type>(),a,b); }

NAMESPACE_END(Grid);
#endif
