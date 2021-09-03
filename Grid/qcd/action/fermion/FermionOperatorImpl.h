/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/action/fermion/FermionOperatorImpl.h

Copyright (C) 2015

Author: Peter Boyle <pabobyle@ph.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: Peter Boyle <peterboyle@Peters-MacBook-Pro-2.local>
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

#pragma once

NAMESPACE_BEGIN(Grid);

//////////////////////////////////////////////
// Template parameter class constructs to package
// externally control Fermion implementations
// in orthogonal directions
//
// Ultimately need Impl to always define types where XXX is opaque
//
//    typedef typename XXX               Simd;
//    typedef typename XXX     GaugeLinkField;        
//    typedef typename XXX         GaugeField;
//    typedef typename XXX      GaugeActField;
//    typedef typename XXX       FermionField;
//    typedef typename XXX    PropagatorField;
//    typedef typename XXX  DoubledGaugeField;
//    typedef typename XXX         SiteSpinor;
//    typedef typename XXX     SitePropagator;
//    typedef typename XXX     SiteHalfSpinor;	
//    typedef typename XXX         Compressor;	
//
// and Methods:
//    void ImportGauge(GridBase *GaugeGrid,DoubledGaugeField &Uds,const GaugeField &Umu)
//    void DoubleStore(GridBase *GaugeGrid,DoubledGaugeField &Uds,const GaugeField &Umu)
//    void multLink(SiteHalfSpinor &phi,const SiteDoubledGaugeField &U,const SiteHalfSpinor &chi,int mu,StencilEntry *SE,StencilImpl::View_type &St)
//    void InsertForce4D(GaugeField &mat,const FermionField &Btilde,const FermionField &A,int mu)
//    void InsertForce5D(GaugeField &mat,const FermionField &Btilde,const FermionField &A,int mu)
//
//
// To acquire the typedefs from "Base" (either a base class or template param) use:
//
// INHERIT_GIMPL_TYPES(Base)
// INHERIT_FIMPL_TYPES(Base)
// INHERIT_IMPL_TYPES(Base)
//
// The Fermion operators will do the following:
//
// struct MyOpParams { 
//   RealD mass;
// };
//
//
// template<class Impl>
// class MyOp : public<Impl> { 
// public:
//
//    INHERIT_ALL_IMPL_TYPES(Impl);
//
//    MyOp(MyOpParams Myparm, ImplParams &ImplParam) :  Impl(ImplParam)
//    {
//
//    };
//    
//  }
//////////////////////////////////////////////

template <class T> struct SamePrecisionMapper {
  typedef T HigherPrecVector ;
  typedef T LowerPrecVector ;
};
template <class T> struct LowerPrecisionMapper {  };
template <> struct LowerPrecisionMapper<vRealF> {
  typedef vRealF HigherPrecVector ;
  typedef vRealH LowerPrecVector ;
};
template <> struct LowerPrecisionMapper<vRealD> {
  typedef vRealD HigherPrecVector ;
  typedef vRealF LowerPrecVector ;
};
template <> struct LowerPrecisionMapper<vComplexF> {
  typedef vComplexF HigherPrecVector ;
  typedef vComplexH LowerPrecVector ;
};
template <> struct LowerPrecisionMapper<vComplexD> {
  typedef vComplexD HigherPrecVector ;
  typedef vComplexF LowerPrecVector ;
};

struct CoeffReal {
public:
  typedef RealD _Coeff_t;
  static const int Nhcs = 2;
  template<class Simd> using PrecisionMapper = SamePrecisionMapper<Simd>;
};
struct CoeffRealHalfComms {
public:
  typedef RealD _Coeff_t;
  static const int Nhcs = 1;
  template<class Simd> using PrecisionMapper = LowerPrecisionMapper<Simd>;
};
struct CoeffComplex {
public:
  typedef ComplexD _Coeff_t;
  static const int Nhcs = 2;
  template<class Simd> using PrecisionMapper = SamePrecisionMapper<Simd>;
};
struct CoeffComplexHalfComms {
public:
  typedef ComplexD _Coeff_t;
  static const int Nhcs = 1;
  template<class Simd> using PrecisionMapper = LowerPrecisionMapper<Simd>;
};

////////////////////////////////////////////////////////////////////////
// Implementation dependent fermion types
////////////////////////////////////////////////////////////////////////
  
#define INHERIT_FIMPL_TYPES(Impl)\
  typedef typename Impl::Coeff_t                     Coeff_t;           \
  typedef Impl Impl_t;							\
  typedef typename Impl::FermionField           FermionField;		\
  typedef typename Impl::PropagatorField     PropagatorField;		\
  typedef typename Impl::DoubledGaugeField DoubledGaugeField;		\
  typedef typename Impl::SiteDoubledGaugeField SiteDoubledGaugeField;	\
  typedef typename Impl::SiteSpinor               SiteSpinor;		\
  typedef typename Impl::SitePropagator       SitePropagator;		\
  typedef typename Impl::SiteHalfSpinor       SiteHalfSpinor;		\
  typedef typename Impl::Compressor               Compressor;		\
  typedef typename Impl::StencilImpl             StencilImpl;		\
  typedef typename Impl::ImplParams               ImplParams;	        \
  typedef typename Impl::StencilImpl::View_type  StencilView;		\
  typedef const typename ViewMap<FermionField>::Type      FermionFieldView;	\
  typedef const typename ViewMap<DoubledGaugeField>::Type DoubledGaugeFieldView;

#define INHERIT_IMPL_TYPES(Base)		\
  INHERIT_GIMPL_TYPES(Base)			\
  INHERIT_FIMPL_TYPES(Base)

NAMESPACE_END(Grid);
NAMESPACE_CHECK(ImplBase);  
/////////////////////////////////////////////////////////////////////////////
// Single flavour four spinors with colour index
/////////////////////////////////////////////////////////////////////////////
#include <Grid/qcd/action/fermion/WilsonImpl.h> 
NAMESPACE_CHECK(ImplWilson);  
   
////////////////////////////////////////////////////////////////////////////////////////
// Flavour doubled spinors; is Gparity the only? what about C*?
////////////////////////////////////////////////////////////////////////////////////////
#include <Grid/qcd/action/fermion/GparityWilsonImpl.h> 
NAMESPACE_CHECK(ImplGparityWilson);  

/////////////////////////////////////////////////////////////////////////////
// Single flavour one component spinors with colour index
/////////////////////////////////////////////////////////////////////////////
#include <Grid/qcd/action/fermion/StaggeredImpl.h> 
NAMESPACE_CHECK(ImplStaggered);  

/////////////////////////////////////////////////////////////////////////////
// Single flavour one component spinors with colour index. 5d vec
/////////////////////////////////////////////////////////////////////////////
// Deprecate Vec5d
//#include <Grid/qcd/action/fermion/StaggeredVec5dImpl.h> 
//NAMESPACE_CHECK(ImplStaggered5dVec);  


