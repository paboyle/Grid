/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/action/fermion/StaggerdKernelsHand.cc

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
#include <Grid.h>

#define REGISTER

#define LOAD_CHI(b)				\
  const SiteSpinor & ref (b[offset]);		\
  Chi_0=ref()()(0);				\
  Chi_1=ref()()(1);				\
  Chi_2=ref()()(2);


// To splat or not to splat depends on the implementation
#define MULT(A,UChi)				\
  auto & ref(U._odata[sU](A));			\
  Impl::loadLinkElement(U_00,ref()(0,0));	\
  Impl::loadLinkElement(U_10,ref()(1,0));	\
  Impl::loadLinkElement(U_20,ref()(2,0));	\
  Impl::loadLinkElement(U_01,ref()(0,1));	\
  Impl::loadLinkElement(U_11,ref()(1,1));	\
  Impl::loadLinkElement(U_21,ref()(2,1));	\
  Impl::loadLinkElement(U_02,ref()(0,2));	\
  Impl::loadLinkElement(U_12,ref()(1,2));	\
  Impl::loadLinkElement(U_22,ref()(2,2));	\
  UChi ## _0  = U_00*Chi_0;			\
  UChi ## _1  = U_10*Chi_0;			\
  UChi ## _2  = U_20*Chi_0;			\
  UChi ## _0 += U_01*Chi_1;			\
  UChi ## _1 += U_11*Chi_1;			\
  UChi ## _2 += U_21*Chi_1;			\
  UChi ## _0 += U_02*Chi_2;			\
  UChi ## _1 += U_12*Chi_2;			\
  UChi ## _2 += U_22*Chi_2;

#define MULT_ADD(A,UChi)			\
  auto & ref(U._odata[sU](A));			\
  Impl::loadLinkElement(U_00,ref()(0,0));	\
  Impl::loadLinkElement(U_10,ref()(1,0));	\
  Impl::loadLinkElement(U_20,ref()(2,0));	\
  Impl::loadLinkElement(U_01,ref()(0,1));	\
  Impl::loadLinkElement(U_11,ref()(1,1));	\
  Impl::loadLinkElement(U_21,ref()(2,1));	\
  Impl::loadLinkElement(U_02,ref()(0,2));	\
  Impl::loadLinkElement(U_12,ref()(1,2));	\
  Impl::loadLinkElement(U_22,ref()(2,2));	\
  UChi ## _0 += U_00*Chi_0;			\
  UChi ## _1 += U_10*Chi_0;			\
  UChi ## _2 += U_20*Chi_0;			\
  UChi ## _0 += U_01*Chi_1;			\
  UChi ## _1 += U_11*Chi_1;			\
  UChi ## _2 += U_21*Chi_1;			\
  UChi ## _0 += U_02*Chi_2;			\
  UChi ## _1 += U_12*Chi_2;			\
  UChi ## _2 += U_22*Chi_2;


#define PERMUTE_DIR(dir)			\
  permute##dir(Chi_0,Chi_0);			\
  permute##dir(Chi_1,Chi_1);			\
  permute##dir(Chi_2,Chi_2);

NAMESPACE_BEGIN(Grid);

template <class Impl>
void StaggeredKernels<Impl>::DhopSiteHand(StencilImpl &st, LebesgueOrder &lo, DoubledGaugeField &U,DoubledGaugeField &UUU,
					  SiteSpinor *buf, int LLs,
					  int sU, const FermionField &in, FermionField &out, int dag) 
{
  SiteSpinor naik; 
  SiteSpinor naive;
  int oneLink  =0;
  int threeLink=1;
  Real scale(1.0);
  
  if(dag) scale = -1.0;
  
  for(int s=0;s<LLs;s++){
    int sF=s+LLs*sU;
    DhopSiteDepthHand(st,lo,U,buf,sF,sU,in,naive,oneLink);
    DhopSiteDepthHand(st,lo,UUU,buf,sF,sU,in,naik,threeLink);
    out._odata[sF] =scale*(naive+naik);
  }
}

template <class Impl>
void StaggeredKernels<Impl>::DhopSiteDepthHand(StencilImpl &st, LebesgueOrder &lo, DoubledGaugeField &U,
					       SiteSpinor *buf, int sF,
					       int sU, const FermionField &in, SiteSpinor &out,int threeLink) 
{
  typedef typename Simd::scalar_type S;
  typedef typename Simd::vector_type V;

  REGISTER Simd even_0; // 12 regs on knc
  REGISTER Simd even_1;
  REGISTER Simd even_2;
  REGISTER Simd odd_0; // 12 regs on knc
  REGISTER Simd odd_1;
  REGISTER Simd odd_2;

  REGISTER Simd Chi_0;    // two spinor; 6 regs
  REGISTER Simd Chi_1;
  REGISTER Simd Chi_2;

  REGISTER Simd U_00;  // two rows of U matrix
  REGISTER Simd U_10;
  REGISTER Simd U_20;  
  REGISTER Simd U_01;
  REGISTER Simd U_11;
  REGISTER Simd U_21;  // 2 reg left.
  REGISTER Simd U_02;
  REGISTER Simd U_12;
  REGISTER Simd U_22; 

  int skew = 0;
  if (threeLink) skew=8;

  int offset,local,perm, ptype;
  StencilEntry *SE;

  // Xp
  SE=st.GetEntry(ptype,Xp+skew,sF);
  offset = SE->_offset;
  local  = SE->_is_local;
  perm   = SE->_permute;
  
  if ( local ) {
    LOAD_CHI(in._odata);
    if ( perm) {
      PERMUTE_DIR(3); // T==0, Z==1, Y==2, Z==3 expect 1,2,2,2 simd layout etc...
    }
  } else { 
    LOAD_CHI(buf);
  }
  {
    MULT(Xp,even);
  }
  
  // Yp
  SE=st.GetEntry(ptype,Yp+skew,sF);
  offset = SE->_offset;
  local  = SE->_is_local;
  perm   = SE->_permute;
  
  if ( local ) {
    LOAD_CHI(in._odata);
    if ( perm) {
      PERMUTE_DIR(2); // T==0, Z==1, Y==2, Z==3 expect 1,2,2,2 simd layout etc...
    }
  } else { 
    LOAD_CHI(buf);
  }
  {
    MULT(Yp,odd);
  }


  // Zp
  SE=st.GetEntry(ptype,Zp+skew,sF);
  offset = SE->_offset;
  local  = SE->_is_local;
  perm   = SE->_permute;
  
  if ( local ) {
    LOAD_CHI(in._odata);
    if ( perm) {
      PERMUTE_DIR(1); // T==0, Z==1, Y==2, Z==3 expect 1,2,2,2 simd layout etc...
    }
  } else { 
    LOAD_CHI(buf);
  }
  {
    MULT_ADD(Zp,even);
  }

  // Tp
  SE=st.GetEntry(ptype,Tp+skew,sF);
  offset = SE->_offset;
  local  = SE->_is_local;
  perm   = SE->_permute;
  
  if ( local ) {
    LOAD_CHI(in._odata);
    if ( perm) {
      PERMUTE_DIR(0); // T==0, Z==1, Y==2, Z==3 expect 1,2,2,2 simd layout etc...
    }
  } else { 
    LOAD_CHI(buf);
  }
  {
    MULT_ADD(Tp,odd);
  }
  
  // Xm
  SE=st.GetEntry(ptype,Xm+skew,sF);
  offset = SE->_offset;
  local  = SE->_is_local;
  perm   = SE->_permute;
  
  if ( local ) {
    LOAD_CHI(in._odata);
    if ( perm) {
      PERMUTE_DIR(3); // T==0, Z==1, Y==2, Z==3 expect 1,2,2,2 simd layout etc...
    }
  } else { 
    LOAD_CHI(buf);
  }
  {
    MULT_ADD(Xm,even);
  }
  
  
  // Ym
  SE=st.GetEntry(ptype,Ym+skew,sF);
  offset = SE->_offset;
  local  = SE->_is_local;
  perm   = SE->_permute;
  
  if ( local ) {
    LOAD_CHI(in._odata);
    if ( perm) {
      PERMUTE_DIR(2); // T==0, Z==1, Y==2, Z==3 expect 1,2,2,2 simd layout etc...
    }
  } else { 
    LOAD_CHI(buf);
  }
  {
    MULT_ADD(Ym,odd);
  }

  // Zm
  SE=st.GetEntry(ptype,Zm+skew,sF);
  offset = SE->_offset;
  local  = SE->_is_local;
  perm   = SE->_permute;
  
  if ( local ) {
    LOAD_CHI(in._odata);
    if ( perm) {
      PERMUTE_DIR(1); // T==0, Z==1, Y==2, Z==3 expect 1,2,2,2 simd layout etc...
    }
  } else { 
    LOAD_CHI(buf);
  }
  {
    MULT_ADD(Zm,even);
  }

  // Tm
  SE=st.GetEntry(ptype,Tm+skew,sF);
  offset = SE->_offset;
  local  = SE->_is_local;
  perm   = SE->_permute;
  
  if ( local ) {
    LOAD_CHI(in._odata);
    if ( perm) {
      PERMUTE_DIR(0); // T==0, Z==1, Y==2, Z==3 expect 1,2,2,2 simd layout etc...
    }
  } else { 
    LOAD_CHI(buf);
  }
  {
    MULT_ADD(Tm,odd);
  }

  vstream(out()()(0),even_0+odd_0);
  vstream(out()()(1),even_1+odd_1);
  vstream(out()()(2),even_2+odd_2);

}

#define DHOP_SITE_HAND_INSTANTIATE(IMPL)				\
  template void StaggeredKernels<IMPL>::DhopSiteHand(StencilImpl &st, LebesgueOrder &lo, \
						     DoubledGaugeField &U,DoubledGaugeField &UUU, \
						     SiteSpinor *buf, int LLs, \
						     int sU, const FermionField &in, FermionField &out, int dag);

#define DHOP_SITE_DEPTH_HAND_INSTANTIATE(IMPL)				\
  template void StaggeredKernels<IMPL>::DhopSiteDepthHand(StencilImpl &st, LebesgueOrder &lo, DoubledGaugeField &U, \
							  SiteSpinor *buf, int sF, \
							  int sU, const FermionField &in, SiteSpinor &out,int threeLink) ;
DHOP_SITE_HAND_INSTANTIATE(StaggeredImplD);
DHOP_SITE_HAND_INSTANTIATE(StaggeredImplF);
DHOP_SITE_HAND_INSTANTIATE(StaggeredVec5dImplD);
DHOP_SITE_HAND_INSTANTIATE(StaggeredVec5dImplF);

DHOP_SITE_DEPTH_HAND_INSTANTIATE(StaggeredImplD);
DHOP_SITE_DEPTH_HAND_INSTANTIATE(StaggeredImplF);
DHOP_SITE_DEPTH_HAND_INSTANTIATE(StaggeredVec5dImplD);
DHOP_SITE_DEPTH_HAND_INSTANTIATE(StaggeredVec5dImplF);

NAMESPACE_END(Grid);
