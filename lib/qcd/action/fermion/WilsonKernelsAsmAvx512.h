/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 



    Source file: ./lib/qcd/action/fermion/WilsonKernelsAsmAvx512.h

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
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


#if defined(AVX512) 
    ///////////////////////////////////////////////////////////
    // If we are AVX512 specialise the single precision routine
    ///////////////////////////////////////////////////////////
#include <simd/Intel512wilson.h>
#include <simd/Intel512single.h>
    
static Vector<vComplexF> signsF;

  template<typename vtype>    
  int setupSigns(Vector<vtype>& signs ){
    Vector<vtype> bother(2);
    signs = bother;
    vrsign(signs[0]);
    visign(signs[1]);
    return 1;
  }

  static int signInitF = setupSigns(signsF);

#define MAYBEPERM(A,perm) if (perm) { A ; }
#define MULT_2SPIN(ptr,pf) MULT_ADDSUB_2SPIN(ptr,pf)
#define COMPLEX_SIGNS(isigns) vComplexF *isigns = &signsF[0];  
  
/////////////////////////////////////////////////////////////////
// XYZT vectorised, undag Kernel, single
/////////////////////////////////////////////////////////////////
#undef KERNEL_DAG
#define INTERIOR_AND_EXTERIOR
#undef INTERIOR
#undef EXTERIOR
template<> void 
WilsonKernels<WilsonImplF>::AsmDhopSite(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

template<> void 
WilsonKernels<ZWilsonImplF>::AsmDhopSite(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

template<> void 
WilsonKernels<WilsonImplFH>::AsmDhopSite(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

template<> void 
WilsonKernels<ZWilsonImplFH>::AsmDhopSite(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

#undef INTERIOR_AND_EXTERIOR
#define INTERIOR
#undef EXTERIOR
template<> void 
WilsonKernels<WilsonImplF>::AsmDhopSiteInt(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

template<> void 
WilsonKernels<ZWilsonImplF>::AsmDhopSiteInt(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

template<> void 
WilsonKernels<WilsonImplFH>::AsmDhopSiteInt(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

template<> void 
WilsonKernels<ZWilsonImplFH>::AsmDhopSiteInt(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>


#undef INTERIOR_AND_EXTERIOR
#undef INTERIOR
#define EXTERIOR
template<> void 
WilsonKernels<WilsonImplF>::AsmDhopSiteExt(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

template<> void 
WilsonKernels<ZWilsonImplF>::AsmDhopSiteExt(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

template<> void 
WilsonKernels<WilsonImplFH>::AsmDhopSiteExt(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

template<> void 
WilsonKernels<ZWilsonImplFH>::AsmDhopSiteExt(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
      
/////////////////////////////////////////////////////////////////
// XYZT vectorised, dag Kernel, single
/////////////////////////////////////////////////////////////////
#define KERNEL_DAG
#define INTERIOR_AND_EXTERIOR
#undef INTERIOR
#undef EXTERIOR
template<> void 
WilsonKernels<WilsonImplF>::AsmDhopSiteDag(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

template<> void 
WilsonKernels<ZWilsonImplF>::AsmDhopSiteDag(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

template<> void 
WilsonKernels<WilsonImplFH>::AsmDhopSiteDag(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

template<> void 
WilsonKernels<ZWilsonImplFH>::AsmDhopSiteDag(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

#undef INTERIOR_AND_EXTERIOR
#define INTERIOR
#undef EXTERIOR
template<> void 
WilsonKernels<WilsonImplF>::AsmDhopSiteDagInt(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

template<> void 
WilsonKernels<ZWilsonImplF>::AsmDhopSiteDagInt(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

template<> void 
WilsonKernels<WilsonImplFH>::AsmDhopSiteDagInt(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

template<> void 
WilsonKernels<ZWilsonImplFH>::AsmDhopSiteDagInt(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

#undef INTERIOR_AND_EXTERIOR
#undef INTERIOR
#define EXTERIOR
template<> void 
WilsonKernels<WilsonImplF>::AsmDhopSiteDagExt(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
				    
template<> void 
WilsonKernels<ZWilsonImplF>::AsmDhopSiteDagExt(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
				    
template<> void 
WilsonKernels<WilsonImplFH>::AsmDhopSiteDagExt(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
				    
template<> void 
WilsonKernels<ZWilsonImplFH>::AsmDhopSiteDagExt(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
				    
#undef MAYBEPERM
#undef MULT_2SPIN
#define MAYBEPERM(A,B) 
#define MULT_2SPIN(ptr,pf) MULT_ADDSUB_2SPIN_LS(ptr,pf)
				    
/////////////////////////////////////////////////////////////////
// Ls vectorised, undag Kernel, single
/////////////////////////////////////////////////////////////////
#undef KERNEL_DAG
#define INTERIOR_AND_EXTERIOR
#undef INTERIOR
#undef EXTERIOR
template<> void 
WilsonKernels<DomainWallVec5dImplF>::AsmDhopSite(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
							 int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
template<> void 
WilsonKernels<ZDomainWallVec5dImplF>::AsmDhopSite(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
							 int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

template<> void 
WilsonKernels<DomainWallVec5dImplFH>::AsmDhopSite(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
							 int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
template<> void 
WilsonKernels<ZDomainWallVec5dImplFH>::AsmDhopSite(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
							 int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

#undef INTERIOR_AND_EXTERIOR
#define INTERIOR
#undef EXTERIOR
template<> void 
WilsonKernels<DomainWallVec5dImplF>::AsmDhopSiteInt(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
							 int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
template<> void 
WilsonKernels<ZDomainWallVec5dImplF>::AsmDhopSiteInt(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
							 int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

template<> void 
WilsonKernels<DomainWallVec5dImplFH>::AsmDhopSiteInt(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
							 int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
template<> void 
WilsonKernels<ZDomainWallVec5dImplFH>::AsmDhopSiteInt(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
							 int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

#undef INTERIOR_AND_EXTERIOR
#undef INTERIOR
#define EXTERIOR
#undef  MULT_2SPIN
#define MULT_2SPIN(ptr,pf) MULT_ADDSUB_2SPIN_LSNOPF(ptr,pf)
template<> void 
WilsonKernels<DomainWallVec5dImplF>::AsmDhopSiteExt(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
							 int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
				    
template<> void 
WilsonKernels<ZDomainWallVec5dImplF>::AsmDhopSiteExt(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
							 int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
				    
template<> void 
WilsonKernels<DomainWallVec5dImplFH>::AsmDhopSiteExt(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
							 int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
				    
template<> void 
WilsonKernels<ZDomainWallVec5dImplFH>::AsmDhopSiteExt(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
							 int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
				    
/////////////////////////////////////////////////////////////////
// Ls vectorised, dag Kernel, single
/////////////////////////////////////////////////////////////////
#define KERNEL_DAG
#define INTERIOR_AND_EXTERIOR
#undef INTERIOR
#undef EXTERIOR
template<> void 
WilsonKernels<DomainWallVec5dImplF>::AsmDhopSiteDag(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U,SiteHalfSpinor *buf,
							    int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
template<> void 
WilsonKernels<ZDomainWallVec5dImplF>::AsmDhopSiteDag(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U,SiteHalfSpinor *buf,
							    int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

template<> void 
WilsonKernels<DomainWallVec5dImplFH>::AsmDhopSiteDag(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U,SiteHalfSpinor *buf,
							    int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
template<> void 
WilsonKernels<ZDomainWallVec5dImplFH>::AsmDhopSiteDag(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U,SiteHalfSpinor *buf,
							    int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

#undef INTERIOR_AND_EXTERIOR
#define INTERIOR
#undef EXTERIOR
template<> void 
WilsonKernels<DomainWallVec5dImplF>::AsmDhopSiteDagInt(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U,SiteHalfSpinor *buf,
							    int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
template<> void 
WilsonKernels<ZDomainWallVec5dImplF>::AsmDhopSiteDagInt(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U,SiteHalfSpinor *buf,
							    int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

template<> void 
WilsonKernels<DomainWallVec5dImplFH>::AsmDhopSiteDagInt(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U,SiteHalfSpinor *buf,
							    int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
template<> void 
WilsonKernels<ZDomainWallVec5dImplFH>::AsmDhopSiteDagInt(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U,SiteHalfSpinor *buf,
							    int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

#undef INTERIOR_AND_EXTERIOR
#undef INTERIOR
#define EXTERIOR
template<> void 
WilsonKernels<DomainWallVec5dImplF>::AsmDhopSiteDagExt(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U,SiteHalfSpinor *buf,
							    int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
template<> void 
WilsonKernels<ZDomainWallVec5dImplF>::AsmDhopSiteDagExt(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U,SiteHalfSpinor *buf,
							    int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

template<> void 
WilsonKernels<DomainWallVec5dImplFH>::AsmDhopSiteDagExt(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U,SiteHalfSpinor *buf,
							    int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
template<> void 
WilsonKernels<ZDomainWallVec5dImplFH>::AsmDhopSiteDagExt(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U,SiteHalfSpinor *buf,
							    int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

#undef COMPLEX_SIGNS
#undef MAYBEPERM
#undef MULT_2SPIN
	


///////////////////////////////////////////////////////////
// If we are AVX512 specialise the double precision routine
///////////////////////////////////////////////////////////

#include <simd/Intel512double.h>
    
static Vector<vComplexD> signsD;
static int signInitD = setupSigns(signsD);
    
#define MAYBEPERM(A,perm) if (perm) { A ; }
#define MULT_2SPIN(ptr,pf) MULT_ADDSUB_2SPIN(ptr,pf)
#define COMPLEX_SIGNS(isigns) vComplexD *isigns = &signsD[0];  


#define INTERIOR_AND_EXTERIOR    
#undef  INTERIOR
#undef  EXTERIOR
  
/////////////////////////////////////////////////////////////////
// XYZT vectorised, undag Kernel, single
/////////////////////////////////////////////////////////////////
#undef KERNEL_DAG
#define INTERIOR_AND_EXTERIOR
#undef INTERIOR
#undef EXTERIOR
template<> void 
WilsonKernels<WilsonImplD>::AsmDhopSite(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
template<> void 
WilsonKernels<ZWilsonImplD>::AsmDhopSite(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

template<> void 
WilsonKernels<WilsonImplDF>::AsmDhopSite(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
template<> void 
WilsonKernels<ZWilsonImplDF>::AsmDhopSite(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

#undef INTERIOR_AND_EXTERIOR
#define INTERIOR
#undef EXTERIOR
template<> void 
WilsonKernels<WilsonImplD>::AsmDhopSiteInt(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
template<> void 
WilsonKernels<ZWilsonImplD>::AsmDhopSiteInt(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

template<> void 
WilsonKernels<WilsonImplDF>::AsmDhopSiteInt(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
template<> void 
WilsonKernels<ZWilsonImplDF>::AsmDhopSiteInt(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

#undef INTERIOR_AND_EXTERIOR
#undef INTERIOR
#define EXTERIOR
template<> void 
WilsonKernels<WilsonImplD>::AsmDhopSiteExt(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
template<> void 
WilsonKernels<ZWilsonImplD>::AsmDhopSiteExt(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
      
template<> void 
WilsonKernels<WilsonImplDF>::AsmDhopSiteExt(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
template<> void 
WilsonKernels<ZWilsonImplDF>::AsmDhopSiteExt(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
      
/////////////////////////////////////////////////////////////////
// XYZT vectorised, dag Kernel, single
/////////////////////////////////////////////////////////////////
#define KERNEL_DAG
#define INTERIOR_AND_EXTERIOR
#undef INTERIOR
#undef EXTERIOR
template<> void 
WilsonKernels<WilsonImplD>::AsmDhopSiteDag(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
template<> void 
WilsonKernels<ZWilsonImplD>::AsmDhopSiteDag(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

template<> void 
WilsonKernels<WilsonImplDF>::AsmDhopSiteDag(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
template<> void 
WilsonKernels<ZWilsonImplDF>::AsmDhopSiteDag(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

#undef INTERIOR_AND_EXTERIOR
#define INTERIOR
#undef EXTERIOR
template<> void 
WilsonKernels<WilsonImplD>::AsmDhopSiteDagInt(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
template<> void 
WilsonKernels<ZWilsonImplD>::AsmDhopSiteDagInt(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

template<> void 
WilsonKernels<WilsonImplDF>::AsmDhopSiteDagInt(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
template<> void 
WilsonKernels<ZWilsonImplDF>::AsmDhopSiteDagInt(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

#undef INTERIOR_AND_EXTERIOR
#undef INTERIOR
#define EXTERIOR
template<> void 
WilsonKernels<WilsonImplD>::AsmDhopSiteDagExt(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
template<> void 
WilsonKernels<ZWilsonImplD>::AsmDhopSiteDagExt(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
				    
template<> void 
WilsonKernels<WilsonImplDF>::AsmDhopSiteDagExt(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
template<> void 
WilsonKernels<ZWilsonImplDF>::AsmDhopSiteDagExt(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
				    
#undef MAYBEPERM
#undef MULT_2SPIN
#define MAYBEPERM(A,B) 
#define MULT_2SPIN(ptr,pf) MULT_ADDSUB_2SPIN_LS(ptr,pf)
				    
/////////////////////////////////////////////////////////////////
// Ls vectorised, undag Kernel, single
/////////////////////////////////////////////////////////////////
#undef KERNEL_DAG
#define INTERIOR_AND_EXTERIOR
#undef INTERIOR
#undef EXTERIOR
template<> void 
WilsonKernels<DomainWallVec5dImplD>::AsmDhopSite(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
							 int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
template<> void 
WilsonKernels<ZDomainWallVec5dImplD>::AsmDhopSite(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
							 int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

template<> void 
WilsonKernels<DomainWallVec5dImplDF>::AsmDhopSite(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
							 int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
template<> void 
WilsonKernels<ZDomainWallVec5dImplDF>::AsmDhopSite(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
							 int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

#undef INTERIOR_AND_EXTERIOR
#define INTERIOR
#undef EXTERIOR
template<> void 
WilsonKernels<DomainWallVec5dImplD>::AsmDhopSiteInt(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
							 int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
template<> void 
WilsonKernels<ZDomainWallVec5dImplD>::AsmDhopSiteInt(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
							 int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

template<> void 
WilsonKernels<DomainWallVec5dImplDF>::AsmDhopSiteInt(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
							 int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
template<> void 
WilsonKernels<ZDomainWallVec5dImplDF>::AsmDhopSiteInt(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
							 int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

#undef INTERIOR_AND_EXTERIOR
#undef INTERIOR
#define EXTERIOR
#undef  MULT_2SPIN
#define MULT_2SPIN(ptr,pf) MULT_ADDSUB_2SPIN_LSNOPF(ptr,pf)
template<> void 
WilsonKernels<DomainWallVec5dImplD>::AsmDhopSiteExt(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
							 int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
template<> void 
WilsonKernels<ZDomainWallVec5dImplD>::AsmDhopSiteExt(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
							 int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
				    
template<> void 
WilsonKernels<DomainWallVec5dImplDF>::AsmDhopSiteExt(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
							 int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
template<> void 
WilsonKernels<ZDomainWallVec5dImplDF>::AsmDhopSiteExt(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
							 int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
				    
/////////////////////////////////////////////////////////////////
// Ls vectorised, dag Kernel, single
/////////////////////////////////////////////////////////////////
#define KERNEL_DAG
#define INTERIOR_AND_EXTERIOR
#undef INTERIOR
#undef EXTERIOR
template<> void 
WilsonKernels<DomainWallVec5dImplD>::AsmDhopSiteDag(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U,SiteHalfSpinor *buf,
							    int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
template<> void 
WilsonKernels<ZDomainWallVec5dImplD>::AsmDhopSiteDag(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U,SiteHalfSpinor *buf,
							    int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

template<> void 
WilsonKernels<DomainWallVec5dImplDF>::AsmDhopSiteDag(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U,SiteHalfSpinor *buf,
							    int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
template<> void 
WilsonKernels<ZDomainWallVec5dImplDF>::AsmDhopSiteDag(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U,SiteHalfSpinor *buf,
							    int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

#undef INTERIOR_AND_EXTERIOR
#define INTERIOR
#undef EXTERIOR
template<> void 
WilsonKernels<DomainWallVec5dImplD>::AsmDhopSiteDagInt(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U,SiteHalfSpinor *buf,
							    int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
template<> void 
WilsonKernels<ZDomainWallVec5dImplD>::AsmDhopSiteDagInt(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U,SiteHalfSpinor *buf,
							    int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

template<> void 
WilsonKernels<DomainWallVec5dImplDF>::AsmDhopSiteDagInt(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U,SiteHalfSpinor *buf,
							    int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
template<> void 
WilsonKernels<ZDomainWallVec5dImplDF>::AsmDhopSiteDagInt(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U,SiteHalfSpinor *buf,
							    int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

#undef INTERIOR_AND_EXTERIOR
#undef INTERIOR
#define EXTERIOR
template<> void 
WilsonKernels<DomainWallVec5dImplD>::AsmDhopSiteDagExt(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U,SiteHalfSpinor *buf,
							    int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
template<> void 
WilsonKernels<ZDomainWallVec5dImplD>::AsmDhopSiteDagExt(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U,SiteHalfSpinor *buf,
							    int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

template<> void 
WilsonKernels<DomainWallVec5dImplDF>::AsmDhopSiteDagExt(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U,SiteHalfSpinor *buf,
							    int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
template<> void 
WilsonKernels<ZDomainWallVec5dImplDF>::AsmDhopSiteDagExt(typename StencilImpl::View_type &st, DoubledGaugeFieldView &U,SiteHalfSpinor *buf,
							    int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

#undef COMPLEX_SIGNS
#undef MAYBEPERM
#undef MULT_2SPIN

#endif //AVX512
