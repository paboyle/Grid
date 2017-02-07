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
WilsonKernels<WilsonImplF>::AsmDhopSite(StencilImpl &st,LebesgueOrder & lo,DoubledGaugeField &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionField &in, FermionField &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

#undef INTERIOR_AND_EXTERIOR
#define INTERIOR
#undef EXTERIOR
template<> void 
WilsonKernels<WilsonImplF>::AsmDhopSiteInt(StencilImpl &st,LebesgueOrder & lo,DoubledGaugeField &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionField &in, FermionField &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

#undef INTERIOR_AND_EXTERIOR
#undef INTERIOR
#define EXTERIOR
template<> void 
WilsonKernels<WilsonImplF>::AsmDhopSiteExt(StencilImpl &st,LebesgueOrder & lo,DoubledGaugeField &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionField &in, FermionField &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
      
/////////////////////////////////////////////////////////////////
// XYZT vectorised, dag Kernel, single
/////////////////////////////////////////////////////////////////
#define KERNEL_DAG
#define INTERIOR_AND_EXTERIOR
#undef INTERIOR
#undef EXTERIOR
template<> void 
WilsonKernels<WilsonImplF>::AsmDhopSiteDag(StencilImpl &st,LebesgueOrder & lo,DoubledGaugeField &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionField &in, FermionField &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

#undef INTERIOR_AND_EXTERIOR
#define INTERIOR
#undef EXTERIOR
template<> void 
WilsonKernels<WilsonImplF>::AsmDhopSiteDagInt(StencilImpl &st,LebesgueOrder & lo,DoubledGaugeField &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionField &in, FermionField &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

#undef INTERIOR_AND_EXTERIOR
#undef INTERIOR
#define EXTERIOR
template<> void 
WilsonKernels<WilsonImplF>::AsmDhopSiteDagExt(StencilImpl &st,LebesgueOrder & lo,DoubledGaugeField &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionField &in, FermionField &out)
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
WilsonKernels<DomainWallVec5dImplF>::AsmDhopSite(StencilImpl &st,LebesgueOrder & lo,DoubledGaugeField &U, SiteHalfSpinor *buf,
							 int ss,int ssU,int Ls,int Ns,const FermionField &in, FermionField &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

#undef INTERIOR_AND_EXTERIOR
#define INTERIOR
#undef EXTERIOR
template<> void 
WilsonKernels<DomainWallVec5dImplF>::AsmDhopSiteInt(StencilImpl &st,LebesgueOrder & lo,DoubledGaugeField &U, SiteHalfSpinor *buf,
							 int ss,int ssU,int Ls,int Ns,const FermionField &in, FermionField &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

#undef INTERIOR_AND_EXTERIOR
#undef INTERIOR
#define EXTERIOR
#undef  MULT_2SPIN
#define MULT_2SPIN(ptr,pf) MULT_ADDSUB_2SPIN_LSNOPF(ptr,pf)
template<> void 
WilsonKernels<DomainWallVec5dImplF>::AsmDhopSiteExt(StencilImpl &st,LebesgueOrder & lo,DoubledGaugeField &U, SiteHalfSpinor *buf,
							 int ss,int ssU,int Ls,int Ns,const FermionField &in, FermionField &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
				    
/////////////////////////////////////////////////////////////////
// Ls vectorised, dag Kernel, single
/////////////////////////////////////////////////////////////////
#define KERNEL_DAG
#define INTERIOR_AND_EXTERIOR
#undef INTERIOR
#undef EXTERIOR
template<> void 
WilsonKernels<DomainWallVec5dImplF>::AsmDhopSiteDag(StencilImpl &st,LebesgueOrder & lo,DoubledGaugeField &U,SiteHalfSpinor *buf,
							    int ss,int ssU,int Ls,int Ns,const FermionField &in, FermionField &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

#undef INTERIOR_AND_EXTERIOR
#define INTERIOR
#undef EXTERIOR
template<> void 
WilsonKernels<DomainWallVec5dImplF>::AsmDhopSiteDagInt(StencilImpl &st,LebesgueOrder & lo,DoubledGaugeField &U,SiteHalfSpinor *buf,
							    int ss,int ssU,int Ls,int Ns,const FermionField &in, FermionField &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

#undef INTERIOR_AND_EXTERIOR
#undef INTERIOR
#define EXTERIOR
template<> void 
WilsonKernels<DomainWallVec5dImplF>::AsmDhopSiteDagExt(StencilImpl &st,LebesgueOrder & lo,DoubledGaugeField &U,SiteHalfSpinor *buf,
							    int ss,int ssU,int Ls,int Ns,const FermionField &in, FermionField &out)
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
WilsonKernels<WilsonImplD>::AsmDhopSite(StencilImpl &st,LebesgueOrder & lo,DoubledGaugeField &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionField &in, FermionField &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

#undef INTERIOR_AND_EXTERIOR
#define INTERIOR
#undef EXTERIOR
template<> void 
WilsonKernels<WilsonImplD>::AsmDhopSiteInt(StencilImpl &st,LebesgueOrder & lo,DoubledGaugeField &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionField &in, FermionField &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

#undef INTERIOR_AND_EXTERIOR
#undef INTERIOR
#define EXTERIOR
template<> void 
WilsonKernels<WilsonImplD>::AsmDhopSiteExt(StencilImpl &st,LebesgueOrder & lo,DoubledGaugeField &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionField &in, FermionField &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
      
/////////////////////////////////////////////////////////////////
// XYZT vectorised, dag Kernel, single
/////////////////////////////////////////////////////////////////
#define KERNEL_DAG
#define INTERIOR_AND_EXTERIOR
#undef INTERIOR
#undef EXTERIOR
template<> void 
WilsonKernels<WilsonImplD>::AsmDhopSiteDag(StencilImpl &st,LebesgueOrder & lo,DoubledGaugeField &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionField &in, FermionField &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

#undef INTERIOR_AND_EXTERIOR
#define INTERIOR
#undef EXTERIOR
template<> void 
WilsonKernels<WilsonImplD>::AsmDhopSiteDagInt(StencilImpl &st,LebesgueOrder & lo,DoubledGaugeField &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionField &in, FermionField &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

#undef INTERIOR_AND_EXTERIOR
#undef INTERIOR
#define EXTERIOR
template<> void 
WilsonKernels<WilsonImplD>::AsmDhopSiteDagExt(StencilImpl &st,LebesgueOrder & lo,DoubledGaugeField &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionField &in, FermionField &out)
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
WilsonKernels<DomainWallVec5dImplD>::AsmDhopSite(StencilImpl &st,LebesgueOrder & lo,DoubledGaugeField &U, SiteHalfSpinor *buf,
							 int ss,int ssU,int Ls,int Ns,const FermionField &in, FermionField &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

#undef INTERIOR_AND_EXTERIOR
#define INTERIOR
#undef EXTERIOR
template<> void 
WilsonKernels<DomainWallVec5dImplD>::AsmDhopSiteInt(StencilImpl &st,LebesgueOrder & lo,DoubledGaugeField &U, SiteHalfSpinor *buf,
							 int ss,int ssU,int Ls,int Ns,const FermionField &in, FermionField &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

#undef INTERIOR_AND_EXTERIOR
#undef INTERIOR
#define EXTERIOR
#undef  MULT_2SPIN
#define MULT_2SPIN(ptr,pf) MULT_ADDSUB_2SPIN_LSNOPF(ptr,pf)
template<> void 
WilsonKernels<DomainWallVec5dImplD>::AsmDhopSiteExt(StencilImpl &st,LebesgueOrder & lo,DoubledGaugeField &U, SiteHalfSpinor *buf,
							 int ss,int ssU,int Ls,int Ns,const FermionField &in, FermionField &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
				    
/////////////////////////////////////////////////////////////////
// Ls vectorised, dag Kernel, single
/////////////////////////////////////////////////////////////////
#define KERNEL_DAG
#define INTERIOR_AND_EXTERIOR
#undef INTERIOR
#undef EXTERIOR
template<> void 
WilsonKernels<DomainWallVec5dImplD>::AsmDhopSiteDag(StencilImpl &st,LebesgueOrder & lo,DoubledGaugeField &U,SiteHalfSpinor *buf,
							    int ss,int ssU,int Ls,int Ns,const FermionField &in, FermionField &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

#undef INTERIOR_AND_EXTERIOR
#define INTERIOR
#undef EXTERIOR
template<> void 
WilsonKernels<DomainWallVec5dImplD>::AsmDhopSiteDagInt(StencilImpl &st,LebesgueOrder & lo,DoubledGaugeField &U,SiteHalfSpinor *buf,
							    int ss,int ssU,int Ls,int Ns,const FermionField &in, FermionField &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

#undef INTERIOR_AND_EXTERIOR
#undef INTERIOR
#define EXTERIOR
template<> void 
WilsonKernels<DomainWallVec5dImplD>::AsmDhopSiteDagExt(StencilImpl &st,LebesgueOrder & lo,DoubledGaugeField &U,SiteHalfSpinor *buf,
							    int ss,int ssU,int Ls,int Ns,const FermionField &in, FermionField &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>

#undef COMPLEX_SIGNS
#undef MAYBEPERM
#undef MULT_2SPIN

#endif //AVX512
