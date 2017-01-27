/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 



    Source file: ./lib/qcd/action/fermion/WilsonKernelsAsm.cc

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: paboyle <paboyle@ph.ed.ac.uk>
Author: Guido Cossu <guido.cossu@ed.ac.uk>

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

#include <Grid/Grid.h>


namespace Grid {
namespace QCD {
    
///////////////////////////////////////////////////////////
// Default to no assembler implementation
///////////////////////////////////////////////////////////
template<class Impl> void 
WilsonKernels<Impl >::DiracOptAsmDhopSite(StencilImpl &st,LebesgueOrder & lo,DoubledGaugeField &U,SiteHalfSpinor *buf,
					  int ss,int ssU,int Ls,int Ns,const FermionField &in, FermionField &out)
{
  assert(0);
}

template<class Impl> void 
WilsonKernels<Impl >::DiracOptAsmDhopSiteDag(StencilImpl &st,LebesgueOrder & lo,DoubledGaugeField &U,SiteHalfSpinor *buf,
					     int ss,int ssU,int Ls,int Ns,const FermionField &in, FermionField &out)
{
  assert(0);
}

#if defined(AVX512) 
#include <simd/Intel512wilson.h>

    ///////////////////////////////////////////////////////////
    // If we are AVX512 specialise the single precision routine
    ///////////////////////////////////////////////////////////

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
  
#define label(A)  ilabel(A)
#define ilabel(A) ".globl\n"  #A ":\n" 
  
#define MAYBEPERM(A,perm) if (perm) { A ; }
#define MULT_2SPIN(ptr,pf) MULT_ADDSUB_2SPIN(ptr,pf)
#define FX(A) WILSONASM_ ##A
#define COMPLEX_TYPE vComplexF
#define signs signsF
  
#undef KERNEL_DAG
template<> void 
WilsonKernels<WilsonImplF>::DiracOptAsmDhopSite(StencilImpl &st,LebesgueOrder & lo,DoubledGaugeField &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionField &in, FermionField &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
      
#define KERNEL_DAG
template<> void 
WilsonKernels<WilsonImplF>::DiracOptAsmDhopSiteDag(StencilImpl &st,LebesgueOrder & lo,DoubledGaugeField &U,SiteHalfSpinor *buf,
						   int ss,int ssU,int Ls,int Ns,const FermionField &in, FermionField &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
				    
#undef VMOVIDUP
#undef VMOVRDUP
#undef MAYBEPERM
#undef MULT_2SPIN
#undef FX 
#define FX(A) DWFASM_ ## A
#define MAYBEPERM(A,B) 
//#define VMOVIDUP(A,B,C)                                  VBCASTIDUPf(A,B,C)
//#define VMOVRDUP(A,B,C)                                  VBCASTRDUPf(A,B,C)
#define MULT_2SPIN(ptr,pf) MULT_ADDSUB_2SPIN_LS(ptr,pf)
				    
#undef KERNEL_DAG
template<> void 
WilsonKernels<DomainWallVec5dImplF>::DiracOptAsmDhopSite(StencilImpl &st,LebesgueOrder & lo,DoubledGaugeField &U, SiteHalfSpinor *buf,
							 int ss,int ssU,int Ls,int Ns,const FermionField &in, FermionField &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
				    
#define KERNEL_DAG
template<> void 
WilsonKernels<DomainWallVec5dImplF>::DiracOptAsmDhopSiteDag(StencilImpl &st,LebesgueOrder & lo,DoubledGaugeField &U,SiteHalfSpinor *buf,
							    int ss,int ssU,int Ls,int Ns,const FermionField &in, FermionField &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
#undef COMPLEX_TYPE
#undef signs
#undef VMOVRDUP
#undef MAYBEPERM
#undef MULT_2SPIN
#undef FX 
	
///////////////////////////////////////////////////////////
// If we are AVX512 specialise the double precision routine
///////////////////////////////////////////////////////////

#include <simd/Intel512double.h>
    
static Vector<vComplexD> signsD;
#define signs signsD
static int signInitD = setupSigns(signsD);
    
#define MAYBEPERM(A,perm) if (perm) { A ; }
#define MULT_2SPIN(ptr,pf) MULT_ADDSUB_2SPIN(ptr,pf)
#define FX(A) WILSONASM_ ##A
#define COMPLEX_TYPE vComplexD
  
#undef KERNEL_DAG
template<> void 
WilsonKernels<WilsonImplD>::DiracOptAsmDhopSite(StencilImpl &st,LebesgueOrder & lo,DoubledGaugeField &U, SiteHalfSpinor *buf,
						int ss,int ssU,int Ls,int Ns,const FermionField &in, FermionField &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
      
#define KERNEL_DAG
template<> void 
WilsonKernels<WilsonImplD>::DiracOptAsmDhopSiteDag(StencilImpl &st,LebesgueOrder & lo,DoubledGaugeField &U,SiteHalfSpinor *buf,
						   int ss,int ssU,int Ls,int Ns,const FermionField &in, FermionField &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
				    
#undef VMOVIDUP
#undef VMOVRDUP
#undef MAYBEPERM
#undef MULT_2SPIN
#undef FX 
#define FX(A) DWFASM_ ## A
#define MAYBEPERM(A,B) 
//#define VMOVIDUP(A,B,C)                                  VBCASTIDUPd(A,B,C)
//#define VMOVRDUP(A,B,C)                                  VBCASTRDUPd(A,B,C)
#define MULT_2SPIN(ptr,pf) MULT_ADDSUB_2SPIN_LS(ptr,pf)
				    
#undef KERNEL_DAG
template<> void 
WilsonKernels<DomainWallVec5dImplD>::DiracOptAsmDhopSite(StencilImpl &st,LebesgueOrder & lo,DoubledGaugeField &U, SiteHalfSpinor *buf,
							 int ss,int ssU,int Ls,int Ns,const FermionField &in, FermionField &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
				    
#define KERNEL_DAG
template<> void 
WilsonKernels<DomainWallVec5dImplD>::DiracOptAsmDhopSiteDag(StencilImpl &st,LebesgueOrder & lo,DoubledGaugeField &U,SiteHalfSpinor *buf,
							    int ss,int ssU,int Ls,int Ns,const FermionField &in, FermionField &out)
#include <qcd/action/fermion/WilsonKernelsAsmBody.h>
	
#undef COMPLEX_TYPE
#undef signs
#undef VMOVRDUP
#undef MAYBEPERM
#undef MULT_2SPIN
#undef FX 

#endif //AVX512

#define INSTANTIATE_ASM(A)\
template void WilsonKernels<A>::DiracOptAsmDhopSite(StencilImpl &st,LebesgueOrder & lo,DoubledGaugeField &U, SiteHalfSpinor *buf,\
                                  int ss,int ssU,int Ls,int Ns,const FermionField &in, FermionField &out);\
 \
template void WilsonKernels<A>::DiracOptAsmDhopSiteDag(StencilImpl &st,LebesgueOrder & lo,DoubledGaugeField &U, SiteHalfSpinor *buf,\
                                  int ss,int ssU,int Ls,int Ns,const FermionField &in, FermionField &out);\

INSTANTIATE_ASM(WilsonImplF);
INSTANTIATE_ASM(WilsonImplD);
INSTANTIATE_ASM(ZWilsonImplF);
INSTANTIATE_ASM(ZWilsonImplD);
INSTANTIATE_ASM(GparityWilsonImplF);
INSTANTIATE_ASM(GparityWilsonImplD);
INSTANTIATE_ASM(DomainWallVec5dImplF);
INSTANTIATE_ASM(DomainWallVec5dImplD);
INSTANTIATE_ASM(ZDomainWallVec5dImplF);
INSTANTIATE_ASM(ZDomainWallVec5dImplD);

}}

