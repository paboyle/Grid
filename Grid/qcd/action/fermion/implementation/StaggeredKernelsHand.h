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
#include <Grid/Grid.h>

#pragma once

NAMESPACE_BEGIN(Grid);

#ifdef GRID_SIMT

#define LOAD_CHI(ptype,b)			\
  const SiteSpinor & ref (b[offset]);				\
  Chi_0=coalescedReadPermute<ptype>(ref()()(0),perm,lane);	\
  Chi_1=coalescedReadPermute<ptype>(ref()()(1),perm,lane);	\
  Chi_2=coalescedReadPermute<ptype>(ref()()(2),perm,lane);

#define LOAD_CHI_COMMS(b)		\
  const SiteSpinor & ref (b[offset]);	\
  Chi_0=coalescedRead(ref()()(0),lane);	\
  Chi_1=coalescedRead(ref()()(1),lane);	\
  Chi_2=coalescedRead(ref()()(2),lane);

#define PERMUTE_DIR(dir)	;
#else
#define LOAD_CHI(ptype,b)      LOAD_CHI_COMMS(b)

#define LOAD_CHI_COMMS(b)		\
  const SiteSpinor & ref (b[offset]);	\
  Chi_0=ref()()(0);			\
  Chi_1=ref()()(1);			\
  Chi_2=ref()()(2);

#define PERMUTE_DIR(dir)			\
  permute##dir(Chi_0,Chi_0);			\
  permute##dir(Chi_1,Chi_1);			\
  permute##dir(Chi_2,Chi_2);

#endif


// To splat or not to splat depends on the implementation
#define MULT(A,UChi)				\
  auto & ref(U[sU](A));			\
    U_00=coalescedRead(ref()(0,0),lane);				\
    U_10=coalescedRead(ref()(1,0),lane);				\
    U_20=coalescedRead(ref()(2,0),lane);				\
    U_01=coalescedRead(ref()(0,1),lane);				\
    U_11=coalescedRead(ref()(1,1),lane);				\
    U_21=coalescedRead(ref()(2,1),lane);				\
    U_02=coalescedRead(ref()(0,2),lane);				\
    U_12=coalescedRead(ref()(1,2),lane);				\
    U_22=coalescedRead(ref()(2,2),lane);				\
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
  auto & ref(U[sU](A));			\
    U_00=coalescedRead(ref()(0,0),lane);				\
    U_10=coalescedRead(ref()(1,0),lane);				\
    U_20=coalescedRead(ref()(2,0),lane);				\
    U_01=coalescedRead(ref()(0,1),lane);				\
    U_11=coalescedRead(ref()(1,1),lane);				\
    U_21=coalescedRead(ref()(2,1),lane);				\
    U_02=coalescedRead(ref()(0,2),lane);				\
    U_12=coalescedRead(ref()(1,2),lane);				\
    U_22=coalescedRead(ref()(2,2),lane);				\
    UChi ## _0 += U_00*Chi_0;	       \
    UChi ## _1 += U_10*Chi_0;\
    UChi ## _2 += U_20*Chi_0;\
    UChi ## _0 += U_01*Chi_1;\
    UChi ## _1 += U_11*Chi_1;\
    UChi ## _2 += U_21*Chi_1;\
    UChi ## _0 += U_02*Chi_2;\
    UChi ## _1 += U_12*Chi_2;\
    UChi ## _2 += U_22*Chi_2;


#define HAND_STENCIL_LEG_BASE(Dir,Perm,skew)	\
  SE=st.GetEntry(ptype,Dir+skew,sF);	\
  offset = SE->_offset;			\
  local  = SE->_is_local;		\
  perm   = SE->_permute;		\
  if ( local ) {						\
    LOAD_CHI(Perm,in);						\
    if ( perm) {						\
      PERMUTE_DIR(Perm);					\
    }								\
  } else {							\
    LOAD_CHI_COMMS(buf);					\
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
    LOAD_CHI(Perm,in);					\
    if ( perm) {					\
      PERMUTE_DIR(Perm);				\
    }							\
  } else if ( st.same_node[Dir] ) {			\
    LOAD_CHI_COMMS(buf);				\
  }							\
  if (local || st.same_node[Dir] ) {		\
    MULT_ADD(U,Dir,even);				\
  }

#define HAND_STENCIL_LEG_EXT(U,Dir,Perm,skew,even)	\
  SE=st.GetEntry(ptype,Dir+skew,sF);			\
  offset = SE->_offset;					\
  local  = SE->_is_local;				\
  if ((!local) && (!st.same_node[Dir]) ) {		\
    nmu++;							\
    { LOAD_CHI_COMMS(buf);	  }				\
    { MULT_ADD(U,Dir,even); }					\
  }								

#define HAND_DECLARATIONS(Simd) \
  Simd even_0;			\
  Simd even_1;			\
  Simd even_2;			\
  Simd odd_0;			\
  Simd odd_1;			\
  Simd odd_2;		        \
		      		\
  Simd Chi_0;			\
  Simd Chi_1;			\
  Simd Chi_2;			\
				\
  Simd U_00;			\
  Simd U_10;			\
  Simd U_20;			\
  Simd U_01;			\
  Simd U_11;			\
  Simd U_21;			\
  Simd U_02;			\
  Simd U_12;			\
  Simd U_22;			
  

template <class Impl>
template <int Naik> accelerator_inline
void StaggeredKernels<Impl>::DhopSiteHand(StencilView &st,
					  DoubledGaugeFieldView &U,DoubledGaugeFieldView &UUU,
					  SiteSpinor *buf, int sF, int sU, 
					  const FermionFieldView &in, FermionFieldView &out,int dag) 
{
  typedef typename Simd::scalar_type S;
  typedef typename Simd::vector_type V;


  const int Nsimd = SiteHalfSpinor::Nsimd();
  const int lane=acceleratorSIMTlane(Nsimd);
  typedef decltype( coalescedRead( in[0]()()(0) )) Simt;
  HAND_DECLARATIONS(Simt);

  typedef decltype( coalescedRead( in[0] )) calcSiteSpinor;
  calcSiteSpinor result;
  int offset,local,perm, ptype;

  StencilEntry *SE;
  int skew;

  //  for(int s=0;s<LLs;s++){
  //    int sF=s+LLs*sU;
  {

    skew = 0;
    HAND_STENCIL_LEG_BEGIN(Xp,3,skew,even);  
    HAND_STENCIL_LEG_BEGIN(Yp,2,skew,odd);   
    HAND_STENCIL_LEG      (U,Zp,1,skew,even);  
    HAND_STENCIL_LEG      (U,Tp,0,skew,odd);  
    HAND_STENCIL_LEG      (U,Xm,3,skew,even);  
    HAND_STENCIL_LEG      (U,Ym,2,skew,odd);   
    HAND_STENCIL_LEG      (U,Zm,1,skew,even);  
    HAND_STENCIL_LEG      (U,Tm,0,skew,odd);  
    if (Naik) {
    skew = 8;
    HAND_STENCIL_LEG(UUU,Xp,3,skew,even);  
    HAND_STENCIL_LEG(UUU,Yp,2,skew,odd);   
    HAND_STENCIL_LEG(UUU,Zp,1,skew,even);  
    HAND_STENCIL_LEG(UUU,Tp,0,skew,odd);  
    HAND_STENCIL_LEG(UUU,Xm,3,skew,even);  
    HAND_STENCIL_LEG(UUU,Ym,2,skew,odd);   
    HAND_STENCIL_LEG(UUU,Zm,1,skew,even);  
    HAND_STENCIL_LEG(UUU,Tm,0,skew,odd);  
    }    
    if ( dag ) {
      result()()(0) = - even_0 - odd_0;
      result()()(1) = - even_1 - odd_1;
      result()()(2) = - even_2 - odd_2;
    } else { 
      result()()(0) = even_0 + odd_0;
      result()()(1) = even_1 + odd_1;
      result()()(2) = even_2 + odd_2;
    }
    coalescedWrite(out[sF],result);
  }
}


template <class Impl>
template <int Naik> accelerator_inline
void StaggeredKernels<Impl>::DhopSiteHandInt(StencilView &st, 
					     DoubledGaugeFieldView &U, DoubledGaugeFieldView &UUU,
					     SiteSpinor *buf, int sF, int sU, 
					     const FermionFieldView &in, FermionFieldView &out,int dag) 
{
  typedef typename Simd::scalar_type S;
  typedef typename Simd::vector_type V;

  const int Nsimd = SiteHalfSpinor::Nsimd();
  const int lane=acceleratorSIMTlane(Nsimd);
  typedef decltype( coalescedRead( in[0]()()(0) )) Simt;
  HAND_DECLARATIONS(Simt);

  typedef decltype( coalescedRead( in[0] )) calcSiteSpinor;
  calcSiteSpinor result;
  int offset, ptype, local, perm;

  StencilEntry *SE;
  int skew;

  //  for(int s=0;s<LLs;s++){
  //    int sF=s+LLs*sU;
  {

    zeroit(even_0);    zeroit(even_1);    zeroit(even_2);
    zeroit(odd_0);    zeroit(odd_1);    zeroit(odd_2);

    skew = 0;
    HAND_STENCIL_LEG_INT(U,Xp,3,skew,even);  
    HAND_STENCIL_LEG_INT(U,Yp,2,skew,odd);   
    HAND_STENCIL_LEG_INT(U,Zp,1,skew,even);  
    HAND_STENCIL_LEG_INT(U,Tp,0,skew,odd);  
    HAND_STENCIL_LEG_INT(U,Xm,3,skew,even);  
    HAND_STENCIL_LEG_INT(U,Ym,2,skew,odd);   
    HAND_STENCIL_LEG_INT(U,Zm,1,skew,even);  
    HAND_STENCIL_LEG_INT(U,Tm,0,skew,odd);  
    if (Naik) {
    skew = 8;
    HAND_STENCIL_LEG_INT(UUU,Xp,3,skew,even);  
    HAND_STENCIL_LEG_INT(UUU,Yp,2,skew,odd);   
    HAND_STENCIL_LEG_INT(UUU,Zp,1,skew,even);  
    HAND_STENCIL_LEG_INT(UUU,Tp,0,skew,odd);  
    HAND_STENCIL_LEG_INT(UUU,Xm,3,skew,even);  
    HAND_STENCIL_LEG_INT(UUU,Ym,2,skew,odd);   
    HAND_STENCIL_LEG_INT(UUU,Zm,1,skew,even);  
    HAND_STENCIL_LEG_INT(UUU,Tm,0,skew,odd);  
    }
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
    coalescedWrite(out[sF],result);
  }
}


template <class Impl>
template <int Naik> accelerator_inline
void StaggeredKernels<Impl>::DhopSiteHandExt(StencilView &st,
					     DoubledGaugeFieldView &U, DoubledGaugeFieldView &UUU,
					     SiteSpinor *buf, int sF, int sU, 
					     const FermionFieldView &in, FermionFieldView &out,int dag) 
{
  typedef typename Simd::scalar_type S;
  typedef typename Simd::vector_type V;

  const int Nsimd = SiteHalfSpinor::Nsimd();
  const int lane=acceleratorSIMTlane(Nsimd);
  typedef decltype( coalescedRead( in[0]()()(0) )) Simt;
  HAND_DECLARATIONS(Simt);

  typedef decltype( coalescedRead( in[0] )) calcSiteSpinor;
  calcSiteSpinor result;
  int offset, ptype, local;

  StencilEntry *SE;
  int skew;

  //  for(int s=0;s<LLs;s++){
  //    int sF=s+LLs*sU;
  {

    zeroit(even_0);    zeroit(even_1);    zeroit(even_2);
    zeroit(odd_0);    zeroit(odd_1);    zeroit(odd_2);
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
    if (Naik) {
    skew = 8;
    HAND_STENCIL_LEG_EXT(UUU,Xp,3,skew,even);  
    HAND_STENCIL_LEG_EXT(UUU,Yp,2,skew,odd);   
    HAND_STENCIL_LEG_EXT(UUU,Zp,1,skew,even);  
    HAND_STENCIL_LEG_EXT(UUU,Tp,0,skew,odd);  
    HAND_STENCIL_LEG_EXT(UUU,Xm,3,skew,even);  
    HAND_STENCIL_LEG_EXT(UUU,Ym,2,skew,odd);   
    HAND_STENCIL_LEG_EXT(UUU,Zm,1,skew,even);  
    HAND_STENCIL_LEG_EXT(UUU,Tm,0,skew,odd);  
    }
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
      coalescedWrite(out[sF] , out(sF)+ result);
    }
  }
}

/*
#define DHOP_SITE_HAND_INSTANTIATE(IMPL)				\
  template void StaggeredKernels<IMPL>::DhopSiteHand(StencilImpl &st, LebesgueOrder &lo, \
						     DoubledGaugeFieldView &U,DoubledGaugeFieldView &UUU, \
						     SiteSpinor *buf, int LLs, int sU, \
						     const FermionFieldView &in, FermionFieldView &out, int dag); \
									\
  template void StaggeredKernels<IMPL>::DhopSiteHandInt(StencilImpl &st, LebesgueOrder &lo, \
						     DoubledGaugeFieldView &U,DoubledGaugeFieldView &UUU, \
						     SiteSpinor *buf, int LLs, int sU, \
						     const FermionFieldView &in, FermionFieldView &out, int dag); \
									\
  template void StaggeredKernels<IMPL>::DhopSiteHandExt(StencilImpl &st, LebesgueOrder &lo, \
						     DoubledGaugeFieldView &U,DoubledGaugeFieldView &UUU, \
						     SiteSpinor *buf, int LLs, int sU, \
						     const FermionFieldView &in, FermionFieldView &out, int dag); \
*/
#undef LOAD_CHI
#undef HAND_DECLARATIONS

NAMESPACE_END(Grid);


