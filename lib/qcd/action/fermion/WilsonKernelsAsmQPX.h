/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 



    Source file: ./lib/qcd/action/fermion/WilsonKernelsAsmQPX.h

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


#if defined(QPX) 

    ///////////////////////////////////////////////////////////
    // If we are QPX specialise the single precision routine
    ///////////////////////////////////////////////////////////

#include <simd/IBM_qpx.h>
#include <simd/IBM_qpx_single.h>
  
#define MAYBEPERM(A,perm) if (perm) { A ; }
#define MULT_2SPIN(ptr,pf) MULT_2SPIN_QPX(ptr,pf)
#define COMPLEX_SIGNS(isigns) 

#define INTERIOR_AND_EXTERIOR    
#undef  INTERIOR
#undef  EXTERIOR
  
/////////////////////////////////////////////////////////////////
// XYZT vectorised, undag Kernel, single
/////////////////////////////////////////////////////////////////
#undef KERNEL_DAG
template<> void 
WilsonKernels<WilsonImplF>::AsmDhopSite(StencilImpl &st,LebesgueOrder & lo,DoubledGaugeField &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionField &in, FermionField &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
      
/////////////////////////////////////////////////////////////////
// XYZT vectorised, dag Kernel, single
/////////////////////////////////////////////////////////////////
#define KERNEL_DAG
template<> void 
WilsonKernels<WilsonImplF>::AsmDhopSiteDag(StencilImpl &st,LebesgueOrder & lo,DoubledGaugeField &U,SiteHalfSpinor *buf,
						   int ss,int ssU,int Ls,int Ns,const FermionField &in, FermionField &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
				    
#undef MAYBEPERM
#undef MULT_2SPIN
#define MAYBEPERM(A,B) 
#define MULT_2SPIN(ptr,pf) MULT_2SPIN_QPX_LS(ptr,pf)
				    
/////////////////////////////////////////////////////////////////
// Ls vectorised, undag Kernel, single
/////////////////////////////////////////////////////////////////
#undef KERNEL_DAG
template<> void 
WilsonKernels<DomainWallVec5dImplF>::AsmDhopSite(StencilImpl &st,LebesgueOrder & lo,DoubledGaugeField &U, SiteHalfSpinor *buf,
							 int ss,int ssU,int Ls,int Ns,const FermionField &in, FermionField &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
				    
/////////////////////////////////////////////////////////////////
// Ls vectorised, dag Kernel, single
/////////////////////////////////////////////////////////////////
#define KERNEL_DAG
template<> void 
WilsonKernels<DomainWallVec5dImplF>::AsmDhopSiteDag(StencilImpl &st,LebesgueOrder & lo,DoubledGaugeField &U,SiteHalfSpinor *buf,
							    int ss,int ssU,int Ls,int Ns,const FermionField &in, FermionField &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
#undef MAYBEPERM
#undef MULT_2SPIN
	
///////////////////////////////////////////////////////////
// DP routines
///////////////////////////////////////////////////////////

#include <simd/IBM_qpx_double.h>
    
#define MAYBEPERM(A,perm) if (perm) { A ; }
#define MULT_2SPIN(ptr,pf) MULT_2SPIN_QPX(ptr,pf)

/////////////////////////////////////////////////////////////////
// XYZT Vectorised, undag Kernel, double
/////////////////////////////////////////////////////////////////
#undef KERNEL_DAG
template<> void 
WilsonKernels<WilsonImplD>::AsmDhopSite(StencilImpl &st,LebesgueOrder & lo,DoubledGaugeField &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionField &in, FermionField &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
/////////////////////////////////////////////////////////////////
      

/////////////////////////////////////////////////////////////////
// XYZT Vectorised, dag Kernel, double
/////////////////////////////////////////////////////////////////
#define KERNEL_DAG
template<> void 
WilsonKernels<WilsonImplD>::AsmDhopSiteDag(StencilImpl &st,LebesgueOrder & lo,DoubledGaugeField &U,SiteHalfSpinor *buf,
						   int ss,int ssU,int Ls,int Ns,const FermionField &in, FermionField &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
/////////////////////////////////////////////////////////////////

#undef MAYBEPERM
#undef MULT_2SPIN
#define MAYBEPERM(A,B) 
#define MULT_2SPIN(ptr,pf) MULT_2SPIN_QPX_LS(ptr,pf)
/////////////////////////////////////////////////////////////////
// Ls vectorised, undag Kernel, double
/////////////////////////////////////////////////////////////////
#undef KERNEL_DAG
template<> void 
WilsonKernels<DomainWallVec5dImplD>::AsmDhopSite(StencilImpl &st,LebesgueOrder & lo,DoubledGaugeField &U, SiteHalfSpinor *buf,
							 int ss,int ssU,int Ls,int Ns,const FermionField &in, FermionField &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
/////////////////////////////////////////////////////////////////
				    
/////////////////////////////////////////////////////////////////
// Ls vectorised, dag Kernel, double
/////////////////////////////////////////////////////////////////
#define KERNEL_DAG
template<> void 
WilsonKernels<DomainWallVec5dImplD>::AsmDhopSiteDag(StencilImpl &st,LebesgueOrder & lo,DoubledGaugeField &U,SiteHalfSpinor *buf,
							    int ss,int ssU,int Ls,int Ns,const FermionField &in, FermionField &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
/////////////////////////////////////////////////////////////////
	
#undef MAYBEPERM
#undef MULT_2SPIN

#endif 
