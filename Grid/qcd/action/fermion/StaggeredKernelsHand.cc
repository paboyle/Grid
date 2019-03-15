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


#define LOAD_CHI(b)		\
  const SiteSpinor & ref (b[offset]);	\
    Chi_0=ref()()(0);\
    Chi_1=ref()()(1);\
    Chi_2=ref()()(2);


// To splat or not to splat depends on the implementation
#define MULT(A,UChi)				\
  auto & ref(U._odata[sU](A));			\
   Impl::loadLinkElement(U_00,ref()(0,0));      \
   Impl::loadLinkElement(U_10,ref()(1,0));      \
   Impl::loadLinkElement(U_20,ref()(2,0));      \
   Impl::loadLinkElement(U_01,ref()(0,1));      \
   Impl::loadLinkElement(U_11,ref()(1,1));      \
   Impl::loadLinkElement(U_21,ref()(2,1));      \
   Impl::loadLinkElement(U_02,ref()(0,2));     \
   Impl::loadLinkElement(U_12,ref()(1,2));     \
   Impl::loadLinkElement(U_22,ref()(2,2));     \
    UChi ## _0  = U_00*Chi_0;	       \
    UChi ## _1  = U_10*Chi_0;\
    UChi ## _2  = U_20*Chi_0;\
    UChi ## _0 += U_01*Chi_1;\
    UChi ## _1 += U_11*Chi_1;\
    UChi ## _2 += U_21*Chi_1;\
    UChi ## _0 += U_02*Chi_2;\
    UChi ## _1 += U_12*Chi_2;\
    UChi ## _2 += U_22*Chi_2;

#define MULT_ADD(U,A,UChi)			\
  auto & ref(U._odata[sU](A));			\
   Impl::loadLinkElement(U_00,ref()(0,0));      \
   Impl::loadLinkElement(U_10,ref()(1,0));      \
   Impl::loadLinkElement(U_20,ref()(2,0));      \
   Impl::loadLinkElement(U_01,ref()(0,1));      \
   Impl::loadLinkElement(U_11,ref()(1,1));      \
   Impl::loadLinkElement(U_21,ref()(2,1));      \
   Impl::loadLinkElement(U_02,ref()(0,2));     \
   Impl::loadLinkElement(U_12,ref()(1,2));     \
   Impl::loadLinkElement(U_22,ref()(2,2));     \
    UChi ## _0 += U_00*Chi_0;	       \
    UChi ## _1 += U_10*Chi_0;\
    UChi ## _2 += U_20*Chi_0;\
    UChi ## _0 += U_01*Chi_1;\
    UChi ## _1 += U_11*Chi_1;\
    UChi ## _2 += U_21*Chi_1;\
    UChi ## _0 += U_02*Chi_2;\
    UChi ## _1 += U_12*Chi_2;\
    UChi ## _2 += U_22*Chi_2;


#define PERMUTE_DIR(dir)			\
  permute##dir(Chi_0,Chi_0);			\
  permute##dir(Chi_1,Chi_1);			\
  permute##dir(Chi_2,Chi_2);


#define HAND_STENCIL_LEG_BASE(Dir,Perm,skew)	\
  SE=st.GetEntry(ptype,Dir+skew,sF);	\
  offset = SE->_offset;			\
  local  = SE->_is_local;		\
  perm   = SE->_permute;		\
  if ( local ) {						\
    LOAD_CHI(in._odata);					\
    if ( perm) {						\
      PERMUTE_DIR(Perm);					\
    }								\
  } else {							\
    LOAD_CHI(buf);						\
  }								

#define HAND_STENCIL_LEG_BEGIN(Dir,Perm,skew,even)		\
  HAND_STENCIL_LEG_BASE(Dir,Perm,skew)				\
  {								\
    MULT(Dir,even);						\
  }

#define HAND_STENCIL_LEG(U,Dir,Perm,skew,even)			\
  HAND_STENCIL_LEG_BASE(Dir,Perm,skew)				\
  {								\
    MULT_ADD(U,Dir,even);					\
  }



#define HAND_STENCIL_LEG_INT(U,Dir,Perm,skew,even)	\
  SE=st.GetEntry(ptype,Dir+skew,sF);			\
  offset = SE->_offset;					\
  local  = SE->_is_local;				\
  perm   = SE->_permute;				\
  if ( local ) {					\
    LOAD_CHI(in._odata);				\
    if ( perm) {					\
      PERMUTE_DIR(Perm);				\
    }							\
  } else if ( st.same_node[Dir] ) {			\
    LOAD_CHI(buf);					\
  }							\
  if (SE->_is_local || st.same_node[Dir] ) {		\
    MULT_ADD(U,Dir,even);				\
  }

#define HAND_STENCIL_LEG_EXT(U,Dir,Perm,skew,even)	\
  SE=st.GetEntry(ptype,Dir+skew,sF);			\
  offset = SE->_offset;					\
  local  = SE->_is_local;				\
  perm   = SE->_permute;				\
  if ((!SE->_is_local) && (!st.same_node[Dir]) ) {		\
    nmu++;							\
    { LOAD_CHI(buf);	  }					\
    { MULT_ADD(U,Dir,even); }					\
  }								

namespace Grid {
namespace QCD {


template <class Impl>
void StaggeredKernels<Impl>::DhopSiteHand(StencilImpl &st, LebesgueOrder &lo, 
					  DoubledGaugeField &U,DoubledGaugeField &UUU,
					  SiteSpinor *buf, int LLs, int sU, 
					  const FermionField &in, FermionField &out,int dag) 
{
  typedef typename Simd::scalar_type S;
  typedef typename Simd::vector_type V;

  Simd even_0; // 12 regs on knc
  Simd even_1;
  Simd even_2;
  Simd odd_0; // 12 regs on knc
  Simd odd_1;
  Simd odd_2;

  Simd Chi_0;    // two spinor; 6 regs
  Simd Chi_1;
  Simd Chi_2;
  
  Simd U_00;  // two rows of U matrix
  Simd U_10;
  Simd U_20;  
  Simd U_01;
  Simd U_11;
  Simd U_21;  // 2 reg left.
  Simd U_02;
  Simd U_12;
  Simd U_22; 

  SiteSpinor result;
  int offset,local,perm, ptype;

  StencilEntry *SE;
  int skew;

  for(int s=0;s<LLs;s++){
    int sF=s+LLs*sU;

    skew = 0;
    HAND_STENCIL_LEG_BEGIN(Xp,3,skew,even);  
    HAND_STENCIL_LEG_BEGIN(Yp,2,skew,odd);   
    HAND_STENCIL_LEG      (U,Zp,1,skew,even);  
    HAND_STENCIL_LEG      (U,Tp,0,skew,odd);  
    HAND_STENCIL_LEG      (U,Xm,3,skew,even);  
    HAND_STENCIL_LEG      (U,Ym,2,skew,odd);   
    HAND_STENCIL_LEG      (U,Zm,1,skew,even);  
    HAND_STENCIL_LEG      (U,Tm,0,skew,odd);  
    skew = 8;
    HAND_STENCIL_LEG(UUU,Xp,3,skew,even);  
    HAND_STENCIL_LEG(UUU,Yp,2,skew,odd);   
    HAND_STENCIL_LEG(UUU,Zp,1,skew,even);  
    HAND_STENCIL_LEG(UUU,Tp,0,skew,odd);  
    HAND_STENCIL_LEG(UUU,Xm,3,skew,even);  
    HAND_STENCIL_LEG(UUU,Ym,2,skew,odd);   
    HAND_STENCIL_LEG(UUU,Zm,1,skew,even);  
    HAND_STENCIL_LEG(UUU,Tm,0,skew,odd);  
    
    if ( dag ) {
      result()()(0) = - even_0 - odd_0;
      result()()(1) = - even_1 - odd_1;
      result()()(2) = - even_2 - odd_2;
    } else { 
      result()()(0) = even_0 + odd_0;
      result()()(1) = even_1 + odd_1;
      result()()(2) = even_2 + odd_2;
    }
    vstream(out._odata[sF],result);
  }
}


template <class Impl>
void StaggeredKernels<Impl>::DhopSiteHandInt(StencilImpl &st, LebesgueOrder &lo, 
					     DoubledGaugeField &U, DoubledGaugeField &UUU,
					     SiteSpinor *buf, int LLs, int sU, 
					     const FermionField &in, FermionField &out,int dag) 
{
  typedef typename Simd::scalar_type S;
  typedef typename Simd::vector_type V;

  Simd even_0; // 12 regs on knc
  Simd even_1;
  Simd even_2;
  Simd odd_0; // 12 regs on knc
  Simd odd_1;
  Simd odd_2;

  Simd Chi_0;    // two spinor; 6 regs
  Simd Chi_1;
  Simd Chi_2;
  
  Simd U_00;  // two rows of U matrix
  Simd U_10;
  Simd U_20;  
  Simd U_01;
  Simd U_11;
  Simd U_21;  // 2 reg left.
  Simd U_02;
  Simd U_12;
  Simd U_22; 

  SiteSpinor result;
  int offset,local,perm, ptype;

  StencilEntry *SE;
  int skew;

  for(int s=0;s<LLs;s++){
    int sF=s+LLs*sU;

    even_0 = zero;    even_1 = zero;    even_2 = zero;
     odd_0 = zero;     odd_1 = zero;     odd_2 = zero;

    skew = 0;
    HAND_STENCIL_LEG_INT(U,Xp,3,skew,even);  
    HAND_STENCIL_LEG_INT(U,Yp,2,skew,odd);   
    HAND_STENCIL_LEG_INT(U,Zp,1,skew,even);  
    HAND_STENCIL_LEG_INT(U,Tp,0,skew,odd);  
    HAND_STENCIL_LEG_INT(U,Xm,3,skew,even);  
    HAND_STENCIL_LEG_INT(U,Ym,2,skew,odd);   
    HAND_STENCIL_LEG_INT(U,Zm,1,skew,even);  
    HAND_STENCIL_LEG_INT(U,Tm,0,skew,odd);  
    skew = 8;
    HAND_STENCIL_LEG_INT(UUU,Xp,3,skew,even);  
    HAND_STENCIL_LEG_INT(UUU,Yp,2,skew,odd);   
    HAND_STENCIL_LEG_INT(UUU,Zp,1,skew,even);  
    HAND_STENCIL_LEG_INT(UUU,Tp,0,skew,odd);  
    HAND_STENCIL_LEG_INT(UUU,Xm,3,skew,even);  
    HAND_STENCIL_LEG_INT(UUU,Ym,2,skew,odd);   
    HAND_STENCIL_LEG_INT(UUU,Zm,1,skew,even);  
    HAND_STENCIL_LEG_INT(UUU,Tm,0,skew,odd);  

    // Assume every site must be connected to at least one interior point. No 1^4 subvols.
    if ( dag ) {
      result()()(0) = - even_0 - odd_0;
      result()()(1) = - even_1 - odd_1;
      result()()(2) = - even_2 - odd_2;
    } else { 
      result()()(0) = even_0 + odd_0;
      result()()(1) = even_1 + odd_1;
      result()()(2) = even_2 + odd_2;
    }
    vstream(out._odata[sF],result);
  }
}


template <class Impl>
void StaggeredKernels<Impl>::DhopSiteHandExt(StencilImpl &st, LebesgueOrder &lo, 
					     DoubledGaugeField &U, DoubledGaugeField &UUU,
					     SiteSpinor *buf, int LLs, int sU, 
					     const FermionField &in, FermionField &out,int dag) 
{
  typedef typename Simd::scalar_type S;
  typedef typename Simd::vector_type V;

  Simd even_0; // 12 regs on knc
  Simd even_1;
  Simd even_2;
  Simd odd_0; // 12 regs on knc
  Simd odd_1;
  Simd odd_2;

  Simd Chi_0;    // two spinor; 6 regs
  Simd Chi_1;
  Simd Chi_2;
  
  Simd U_00;  // two rows of U matrix
  Simd U_10;
  Simd U_20;  
  Simd U_01;
  Simd U_11;
  Simd U_21;  // 2 reg left.
  Simd U_02;
  Simd U_12;
  Simd U_22; 

  SiteSpinor result;
  int offset,local,perm, ptype;

  StencilEntry *SE;
  int skew;

  for(int s=0;s<LLs;s++){
    int sF=s+LLs*sU;

    even_0 = zero;    even_1 = zero;    even_2 = zero;
     odd_0 = zero;     odd_1 = zero;     odd_2 = zero;
    int nmu=0;
    skew = 0;
    HAND_STENCIL_LEG_EXT(U,Xp,3,skew,even);  
    HAND_STENCIL_LEG_EXT(U,Yp,2,skew,odd);   
    HAND_STENCIL_LEG_EXT(U,Zp,1,skew,even);  
    HAND_STENCIL_LEG_EXT(U,Tp,0,skew,odd);  
    HAND_STENCIL_LEG_EXT(U,Xm,3,skew,even);  
    HAND_STENCIL_LEG_EXT(U,Ym,2,skew,odd);   
    HAND_STENCIL_LEG_EXT(U,Zm,1,skew,even);  
    HAND_STENCIL_LEG_EXT(U,Tm,0,skew,odd);  
    skew = 8;
    HAND_STENCIL_LEG_EXT(UUU,Xp,3,skew,even);  
    HAND_STENCIL_LEG_EXT(UUU,Yp,2,skew,odd);   
    HAND_STENCIL_LEG_EXT(UUU,Zp,1,skew,even);  
    HAND_STENCIL_LEG_EXT(UUU,Tp,0,skew,odd);  
    HAND_STENCIL_LEG_EXT(UUU,Xm,3,skew,even);  
    HAND_STENCIL_LEG_EXT(UUU,Ym,2,skew,odd);   
    HAND_STENCIL_LEG_EXT(UUU,Zm,1,skew,even);  
    HAND_STENCIL_LEG_EXT(UUU,Tm,0,skew,odd);  

    // Add sum of all exterior connected stencil legs
    if ( nmu ) { 
      if ( dag ) {
	result()()(0) = - even_0 - odd_0;
	result()()(1) = - even_1 - odd_1;
	result()()(2) = - even_2 - odd_2;
      } else { 
	result()()(0) = even_0 + odd_0;
	result()()(1) = even_1 + odd_1;
	result()()(2) = even_2 + odd_2;
      }
      out._odata[sF] = out._odata[sF] + result;
    }
  }
}


#define DHOP_SITE_HAND_INSTANTIATE(IMPL)				\
  template void StaggeredKernels<IMPL>::DhopSiteHand(StencilImpl &st, LebesgueOrder &lo, \
						     DoubledGaugeField &U,DoubledGaugeField &UUU, \
						     SiteSpinor *buf, int LLs, int sU, \
						     const FermionField &in, FermionField &out, int dag); \
									\
  template void StaggeredKernels<IMPL>::DhopSiteHandInt(StencilImpl &st, LebesgueOrder &lo, \
						     DoubledGaugeField &U,DoubledGaugeField &UUU, \
						     SiteSpinor *buf, int LLs, int sU, \
						     const FermionField &in, FermionField &out, int dag); \
									\
  template void StaggeredKernels<IMPL>::DhopSiteHandExt(StencilImpl &st, LebesgueOrder &lo, \
						     DoubledGaugeField &U,DoubledGaugeField &UUU, \
						     SiteSpinor *buf, int LLs, int sU, \
						     const FermionField &in, FermionField &out, int dag); \

DHOP_SITE_HAND_INSTANTIATE(StaggeredImplD);
DHOP_SITE_HAND_INSTANTIATE(StaggeredImplF);
DHOP_SITE_HAND_INSTANTIATE(StaggeredVec5dImplD);
DHOP_SITE_HAND_INSTANTIATE(StaggeredVec5dImplF);


}
}

