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
WilsonKernels<WilsonImplF>::AsmDhopSite(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

template<> void 
WilsonKernels<ZWilsonImplF>::AsmDhopSite(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

template<> void 
WilsonKernels<WilsonImplFH>::AsmDhopSite(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

template<> void 
WilsonKernels<ZWilsonImplFH>::AsmDhopSite(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

#undef INTERIOR_AND_EXTERIOR
#define INTERIOR
#undef EXTERIOR
template<> void 
WilsonKernels<WilsonImplF>::AsmDhopSiteInt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

template<> void 
WilsonKernels<ZWilsonImplF>::AsmDhopSiteInt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

template<> void 
WilsonKernels<WilsonImplFH>::AsmDhopSiteInt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

template<> void 
WilsonKernels<ZWilsonImplFH>::AsmDhopSiteInt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>


#undef INTERIOR_AND_EXTERIOR
#undef INTERIOR
#define EXTERIOR
template<> void 
WilsonKernels<WilsonImplF>::AsmDhopSiteExt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

template<> void 
WilsonKernels<ZWilsonImplF>::AsmDhopSiteExt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

template<> void 
WilsonKernels<WilsonImplFH>::AsmDhopSiteExt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

template<> void 
WilsonKernels<ZWilsonImplFH>::AsmDhopSiteExt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
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
WilsonKernels<WilsonImplF>::AsmDhopSiteDag(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

template<> void 
WilsonKernels<ZWilsonImplF>::AsmDhopSiteDag(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

template<> void 
WilsonKernels<WilsonImplFH>::AsmDhopSiteDag(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

template<> void 
WilsonKernels<ZWilsonImplFH>::AsmDhopSiteDag(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

#undef INTERIOR_AND_EXTERIOR
#define INTERIOR
#undef EXTERIOR
template<> void 
WilsonKernels<WilsonImplF>::AsmDhopSiteDagInt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

template<> void 
WilsonKernels<ZWilsonImplF>::AsmDhopSiteDagInt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

template<> void 
WilsonKernels<WilsonImplFH>::AsmDhopSiteDagInt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

template<> void 
WilsonKernels<ZWilsonImplFH>::AsmDhopSiteDagInt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

#undef INTERIOR_AND_EXTERIOR
#undef INTERIOR
#define EXTERIOR
template<> void 
WilsonKernels<WilsonImplF>::AsmDhopSiteDagExt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
				    
template<> void 
WilsonKernels<ZWilsonImplF>::AsmDhopSiteDagExt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
				    
template<> void 
WilsonKernels<WilsonImplFH>::AsmDhopSiteDagExt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
				    
template<> void 
WilsonKernels<ZWilsonImplFH>::AsmDhopSiteDagExt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
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
WilsonKernels<DomainWallVec5dImplF>::AsmDhopSite(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
							 int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
template<> void 
WilsonKernels<ZDomainWallVec5dImplF>::AsmDhopSite(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
							 int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

template<> void 
WilsonKernels<DomainWallVec5dImplFH>::AsmDhopSite(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
							 int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
template<> void 
WilsonKernels<ZDomainWallVec5dImplFH>::AsmDhopSite(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
							 int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

#undef INTERIOR_AND_EXTERIOR
#define INTERIOR
#undef EXTERIOR
template<> void 
WilsonKernels<DomainWallVec5dImplF>::AsmDhopSiteInt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
							 int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
template<> void 
WilsonKernels<ZDomainWallVec5dImplF>::AsmDhopSiteInt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
							 int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

template<> void 
WilsonKernels<DomainWallVec5dImplFH>::AsmDhopSiteInt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
							 int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
template<> void 
WilsonKernels<ZDomainWallVec5dImplFH>::AsmDhopSiteInt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
							 int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

#undef INTERIOR_AND_EXTERIOR
#undef INTERIOR
#define EXTERIOR
#undef  MULT_2SPIN
#define MULT_2SPIN(ptr,pf) MULT_ADDSUB_2SPIN_LSNOPF(ptr,pf)
template<> void 
WilsonKernels<DomainWallVec5dImplF>::AsmDhopSiteExt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
							 int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
				    
template<> void 
WilsonKernels<ZDomainWallVec5dImplF>::AsmDhopSiteExt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
							 int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
				    
template<> void 
WilsonKernels<DomainWallVec5dImplFH>::AsmDhopSiteExt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
							 int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
				    
template<> void 
WilsonKernels<ZDomainWallVec5dImplFH>::AsmDhopSiteExt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
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
WilsonKernels<DomainWallVec5dImplF>::AsmDhopSiteDag(StencilView &st, DoubledGaugeFieldView &U,SiteHalfSpinor *buf,
							    int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
template<> void 
WilsonKernels<ZDomainWallVec5dImplF>::AsmDhopSiteDag(StencilView &st, DoubledGaugeFieldView &U,SiteHalfSpinor *buf,
							    int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

template<> void 
WilsonKernels<DomainWallVec5dImplFH>::AsmDhopSiteDag(StencilView &st, DoubledGaugeFieldView &U,SiteHalfSpinor *buf,
							    int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
template<> void 
WilsonKernels<ZDomainWallVec5dImplFH>::AsmDhopSiteDag(StencilView &st, DoubledGaugeFieldView &U,SiteHalfSpinor *buf,
							    int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

#undef INTERIOR_AND_EXTERIOR
#define INTERIOR
#undef EXTERIOR
template<> void 
WilsonKernels<DomainWallVec5dImplF>::AsmDhopSiteDagInt(StencilView &st, DoubledGaugeFieldView &U,SiteHalfSpinor *buf,
							    int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
template<> void 
WilsonKernels<ZDomainWallVec5dImplF>::AsmDhopSiteDagInt(StencilView &st, DoubledGaugeFieldView &U,SiteHalfSpinor *buf,
							    int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

template<> void 
WilsonKernels<DomainWallVec5dImplFH>::AsmDhopSiteDagInt(StencilView &st, DoubledGaugeFieldView &U,SiteHalfSpinor *buf,
							    int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
template<> void 
WilsonKernels<ZDomainWallVec5dImplFH>::AsmDhopSiteDagInt(StencilView &st, DoubledGaugeFieldView &U,SiteHalfSpinor *buf,
							    int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

#undef INTERIOR_AND_EXTERIOR
#undef INTERIOR
#define EXTERIOR
template<> void 
WilsonKernels<DomainWallVec5dImplF>::AsmDhopSiteDagExt(StencilView &st, DoubledGaugeFieldView &U,SiteHalfSpinor *buf,
							    int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
template<> void 
WilsonKernels<ZDomainWallVec5dImplF>::AsmDhopSiteDagExt(StencilView &st, DoubledGaugeFieldView &U,SiteHalfSpinor *buf,
							    int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

template<> void 
WilsonKernels<DomainWallVec5dImplFH>::AsmDhopSiteDagExt(StencilView &st, DoubledGaugeFieldView &U,SiteHalfSpinor *buf,
							    int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
template<> void 
WilsonKernels<ZDomainWallVec5dImplFH>::AsmDhopSiteDagExt(StencilView &st, DoubledGaugeFieldView &U,SiteHalfSpinor *buf,
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
WilsonKernels<WilsonImplD>::AsmDhopSite(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
template<> void 
WilsonKernels<ZWilsonImplD>::AsmDhopSite(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

template<> void 
WilsonKernels<WilsonImplDF>::AsmDhopSite(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
template<> void 
WilsonKernels<ZWilsonImplDF>::AsmDhopSite(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

#undef INTERIOR_AND_EXTERIOR
#define INTERIOR
#undef EXTERIOR
template<> void 
WilsonKernels<WilsonImplD>::AsmDhopSiteInt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
template<> void 
WilsonKernels<ZWilsonImplD>::AsmDhopSiteInt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

template<> void 
WilsonKernels<WilsonImplDF>::AsmDhopSiteInt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
template<> void 
WilsonKernels<ZWilsonImplDF>::AsmDhopSiteInt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

#undef INTERIOR_AND_EXTERIOR
#undef INTERIOR
#define EXTERIOR
template<> void 
WilsonKernels<WilsonImplD>::AsmDhopSiteExt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
template<> void 
WilsonKernels<ZWilsonImplD>::AsmDhopSiteExt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
      
template<> void 
WilsonKernels<WilsonImplDF>::AsmDhopSiteExt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
template<> void 
WilsonKernels<ZWilsonImplDF>::AsmDhopSiteExt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
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
WilsonKernels<WilsonImplD>::AsmDhopSiteDag(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
template<> void 
WilsonKernels<ZWilsonImplD>::AsmDhopSiteDag(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

template<> void 
WilsonKernels<WilsonImplDF>::AsmDhopSiteDag(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
template<> void 
WilsonKernels<ZWilsonImplDF>::AsmDhopSiteDag(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

#undef INTERIOR_AND_EXTERIOR
#define INTERIOR
#undef EXTERIOR
template<> void 
WilsonKernels<WilsonImplD>::AsmDhopSiteDagInt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
template<> void 
WilsonKernels<ZWilsonImplD>::AsmDhopSiteDagInt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

template<> void 
WilsonKernels<WilsonImplDF>::AsmDhopSiteDagInt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
template<> void 
WilsonKernels<ZWilsonImplDF>::AsmDhopSiteDagInt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

#undef INTERIOR_AND_EXTERIOR
#undef INTERIOR
#define EXTERIOR
template<> void 
WilsonKernels<WilsonImplD>::AsmDhopSiteDagExt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
template<> void 
WilsonKernels<ZWilsonImplD>::AsmDhopSiteDagExt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
				    
template<> void 
WilsonKernels<WilsonImplDF>::AsmDhopSiteDagExt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
template<> void 
WilsonKernels<ZWilsonImplDF>::AsmDhopSiteDagExt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
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
WilsonKernels<DomainWallVec5dImplD>::AsmDhopSite(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
							 int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
template<> void 
WilsonKernels<ZDomainWallVec5dImplD>::AsmDhopSite(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
							 int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

template<> void 
WilsonKernels<DomainWallVec5dImplDF>::AsmDhopSite(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
							 int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
template<> void 
WilsonKernels<ZDomainWallVec5dImplDF>::AsmDhopSite(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
							 int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

#undef INTERIOR_AND_EXTERIOR
#define INTERIOR
#undef EXTERIOR
template<> void 
WilsonKernels<DomainWallVec5dImplD>::AsmDhopSiteInt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
							 int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
template<> void 
WilsonKernels<ZDomainWallVec5dImplD>::AsmDhopSiteInt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
							 int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

template<> void 
WilsonKernels<DomainWallVec5dImplDF>::AsmDhopSiteInt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
							 int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
template<> void 
WilsonKernels<ZDomainWallVec5dImplDF>::AsmDhopSiteInt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
							 int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

#undef INTERIOR_AND_EXTERIOR
#undef INTERIOR
#define EXTERIOR
#undef  MULT_2SPIN
#define MULT_2SPIN(ptr,pf) MULT_ADDSUB_2SPIN_LSNOPF(ptr,pf)
template<> void 
WilsonKernels<DomainWallVec5dImplD>::AsmDhopSiteExt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
							 int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
template<> void 
WilsonKernels<ZDomainWallVec5dImplD>::AsmDhopSiteExt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
							 int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
				    
template<> void 
WilsonKernels<DomainWallVec5dImplDF>::AsmDhopSiteExt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
							 int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
template<> void 
WilsonKernels<ZDomainWallVec5dImplDF>::AsmDhopSiteExt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,
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
WilsonKernels<DomainWallVec5dImplD>::AsmDhopSiteDag(StencilView &st, DoubledGaugeFieldView &U,SiteHalfSpinor *buf,
							    int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
template<> void 
WilsonKernels<ZDomainWallVec5dImplD>::AsmDhopSiteDag(StencilView &st, DoubledGaugeFieldView &U,SiteHalfSpinor *buf,
							    int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

template<> void 
WilsonKernels<DomainWallVec5dImplDF>::AsmDhopSiteDag(StencilView &st, DoubledGaugeFieldView &U,SiteHalfSpinor *buf,
							    int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
template<> void 
WilsonKernels<ZDomainWallVec5dImplDF>::AsmDhopSiteDag(StencilView &st, DoubledGaugeFieldView &U,SiteHalfSpinor *buf,
							    int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

#undef INTERIOR_AND_EXTERIOR
#define INTERIOR
#undef EXTERIOR
template<> void 
WilsonKernels<DomainWallVec5dImplD>::AsmDhopSiteDagInt(StencilView &st, DoubledGaugeFieldView &U,SiteHalfSpinor *buf,
							    int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
template<> void 
WilsonKernels<ZDomainWallVec5dImplD>::AsmDhopSiteDagInt(StencilView &st, DoubledGaugeFieldView &U,SiteHalfSpinor *buf,
							    int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

template<> void 
WilsonKernels<DomainWallVec5dImplDF>::AsmDhopSiteDagInt(StencilView &st, DoubledGaugeFieldView &U,SiteHalfSpinor *buf,
							    int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
template<> void 
WilsonKernels<ZDomainWallVec5dImplDF>::AsmDhopSiteDagInt(StencilView &st, DoubledGaugeFieldView &U,SiteHalfSpinor *buf,
							    int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

#undef INTERIOR_AND_EXTERIOR
#undef INTERIOR
#define EXTERIOR
template<> void 
WilsonKernels<DomainWallVec5dImplD>::AsmDhopSiteDagExt(StencilView &st, DoubledGaugeFieldView &U,SiteHalfSpinor *buf,
							    int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
template<> void 
WilsonKernels<ZDomainWallVec5dImplD>::AsmDhopSiteDagExt(StencilView &st, DoubledGaugeFieldView &U,SiteHalfSpinor *buf,
							    int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

template<> void 
WilsonKernels<DomainWallVec5dImplDF>::AsmDhopSiteDagExt(StencilView &st, DoubledGaugeFieldView &U,SiteHalfSpinor *buf,
							    int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
template<> void 
WilsonKernels<ZDomainWallVec5dImplDF>::AsmDhopSiteDagExt(StencilView &st, DoubledGaugeFieldView &U,SiteHalfSpinor *buf,
							    int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

#undef COMPLEX_SIGNS
#undef MAYBEPERM
#undef MULT_2SPIN

#endif //AVX512
