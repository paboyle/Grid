/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/action/fermion/WilsonKernels.cc

Copyright (C) 2015

Author: Azusa Yamaguchi, Peter Boyle

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
#include <Grid/qcd/action/fermion/FermionCore.h>

namespace Grid {
namespace QCD {

int StaggeredKernelsStatic::Opt= StaggeredKernelsStatic::OptGeneric;
int StaggeredKernelsStatic::Comms = StaggeredKernelsStatic::CommsAndCompute;

#define GENERIC_STENCIL_LEG(U,Dir,skew,multLink)		\
  SE = st.GetEntry(ptype, Dir+skew, sF);			\
  if (SE->_is_local ) {						\
    if (SE->_permute) {						\
      chi_p = &chi;						\
      permute(chi,  in._odata[SE->_offset], ptype);		\
    } else {							\
      chi_p = &in._odata[SE->_offset];				\
    }								\
  } else {							\
    chi_p = &buf[SE->_offset];					\
  }								\
  multLink(Uchi, U._odata[sU], *chi_p, Dir);			

#define GENERIC_STENCIL_LEG_INT(U,Dir,skew,multLink)		\
  SE = st.GetEntry(ptype, Dir+skew, sF);			\
  if (SE->_is_local ) {						\
    if (SE->_permute) {						\
      chi_p = &chi;						\
      permute(chi,  in._odata[SE->_offset], ptype);		\
    } else {							\
      chi_p = &in._odata[SE->_offset];				\
    }								\
  } else if ( st.same_node[Dir] ) {				\
    chi_p = &buf[SE->_offset];					\
  }								\
  if (SE->_is_local || st.same_node[Dir] ) {			\
    multLink(Uchi, U._odata[sU], *chi_p, Dir);			\
  }

#define GENERIC_STENCIL_LEG_EXT(U,Dir,skew,multLink)		\
  SE = st.GetEntry(ptype, Dir+skew, sF);			\
  if ((!SE->_is_local) && (!st.same_node[Dir]) ) {		\
    nmu++;							\
    chi_p = &buf[SE->_offset];					\
    multLink(Uchi, U._odata[sU], *chi_p, Dir);			\
  }

template <class Impl>
StaggeredKernels<Impl>::StaggeredKernels(const ImplParams &p) : Base(p){};

////////////////////////////////////////////////////////////////////////////////////
// Generic implementation; move to different file?
// Int, Ext, Int+Ext cases for comms overlap
////////////////////////////////////////////////////////////////////////////////////
template <class Impl>
void StaggeredKernels<Impl>::DhopSiteGeneric(StencilImpl &st, LebesgueOrder &lo, 
					     DoubledGaugeField &U, DoubledGaugeField &UUU,
					     SiteSpinor *buf, int LLs, int sU, 
					     const FermionField &in, FermionField &out, int dag) {
  const SiteSpinor *chi_p;
  SiteSpinor chi;
  SiteSpinor Uchi;
  StencilEntry *SE;
  int ptype;
  int skew;

  for(int s=0;s<LLs;s++){
    int sF=LLs*sU+s;
    skew = 0;
    GENERIC_STENCIL_LEG(U,Xp,skew,Impl::multLink);
    GENERIC_STENCIL_LEG(U,Yp,skew,Impl::multLinkAdd);
    GENERIC_STENCIL_LEG(U,Zp,skew,Impl::multLinkAdd);
    GENERIC_STENCIL_LEG(U,Tp,skew,Impl::multLinkAdd);
    GENERIC_STENCIL_LEG(U,Xm,skew,Impl::multLinkAdd);
    GENERIC_STENCIL_LEG(U,Ym,skew,Impl::multLinkAdd);
    GENERIC_STENCIL_LEG(U,Zm,skew,Impl::multLinkAdd);
    GENERIC_STENCIL_LEG(U,Tm,skew,Impl::multLinkAdd);
    skew=8;
    GENERIC_STENCIL_LEG(UUU,Xp,skew,Impl::multLinkAdd);
    GENERIC_STENCIL_LEG(UUU,Yp,skew,Impl::multLinkAdd);
    GENERIC_STENCIL_LEG(UUU,Zp,skew,Impl::multLinkAdd);
    GENERIC_STENCIL_LEG(UUU,Tp,skew,Impl::multLinkAdd);
    GENERIC_STENCIL_LEG(UUU,Xm,skew,Impl::multLinkAdd);
    GENERIC_STENCIL_LEG(UUU,Ym,skew,Impl::multLinkAdd);
    GENERIC_STENCIL_LEG(UUU,Zm,skew,Impl::multLinkAdd);
    GENERIC_STENCIL_LEG(UUU,Tm,skew,Impl::multLinkAdd);
    if ( dag ) { 
      Uchi = - Uchi;
    } 
    vstream(out._odata[sF], Uchi);
  }
};

  ///////////////////////////////////////////////////
  // Only contributions from interior of our node
  ///////////////////////////////////////////////////
template <class Impl>
void StaggeredKernels<Impl>::DhopSiteGenericInt(StencilImpl &st, LebesgueOrder &lo, 
						DoubledGaugeField &U, DoubledGaugeField &UUU,
						SiteSpinor *buf, int LLs, int sU, 
						const FermionField &in, FermionField &out,int dag) {
  const SiteSpinor *chi_p;
  SiteSpinor chi;
  SiteSpinor Uchi;
  StencilEntry *SE;
  int ptype;
  int skew ;

  for(int s=0;s<LLs;s++){
    int sF=LLs*sU+s;
    skew = 0;
    Uchi=zero;
    GENERIC_STENCIL_LEG_INT(U,Xp,skew,Impl::multLinkAdd);
    GENERIC_STENCIL_LEG_INT(U,Yp,skew,Impl::multLinkAdd);
    GENERIC_STENCIL_LEG_INT(U,Zp,skew,Impl::multLinkAdd);
    GENERIC_STENCIL_LEG_INT(U,Tp,skew,Impl::multLinkAdd);
    GENERIC_STENCIL_LEG_INT(U,Xm,skew,Impl::multLinkAdd);
    GENERIC_STENCIL_LEG_INT(U,Ym,skew,Impl::multLinkAdd);
    GENERIC_STENCIL_LEG_INT(U,Zm,skew,Impl::multLinkAdd);
    GENERIC_STENCIL_LEG_INT(U,Tm,skew,Impl::multLinkAdd);
    skew=8;
    GENERIC_STENCIL_LEG_INT(UUU,Xp,skew,Impl::multLinkAdd);
    GENERIC_STENCIL_LEG_INT(UUU,Yp,skew,Impl::multLinkAdd);
    GENERIC_STENCIL_LEG_INT(UUU,Zp,skew,Impl::multLinkAdd);
    GENERIC_STENCIL_LEG_INT(UUU,Tp,skew,Impl::multLinkAdd);
    GENERIC_STENCIL_LEG_INT(UUU,Xm,skew,Impl::multLinkAdd);
    GENERIC_STENCIL_LEG_INT(UUU,Ym,skew,Impl::multLinkAdd);
    GENERIC_STENCIL_LEG_INT(UUU,Zm,skew,Impl::multLinkAdd);
    GENERIC_STENCIL_LEG_INT(UUU,Tm,skew,Impl::multLinkAdd);
    if ( dag ) {
      Uchi = - Uchi;
    }
    vstream(out._odata[sF], Uchi);
  }
};


  ///////////////////////////////////////////////////
  // Only contributions from exterior of our node
  ///////////////////////////////////////////////////
template <class Impl>
void StaggeredKernels<Impl>::DhopSiteGenericExt(StencilImpl &st, LebesgueOrder &lo, 
						DoubledGaugeField &U, DoubledGaugeField &UUU,
						SiteSpinor *buf, int LLs, int sU,
						const FermionField &in, FermionField &out,int dag) {
  const SiteSpinor *chi_p;
  SiteSpinor chi;
  SiteSpinor Uchi;
  StencilEntry *SE;
  int ptype;
  int nmu=0;
  int skew ;

  for(int s=0;s<LLs;s++){
    int sF=LLs*sU+s;
    skew = 0;
    Uchi=zero;
    GENERIC_STENCIL_LEG_EXT(U,Xp,skew,Impl::multLinkAdd);
    GENERIC_STENCIL_LEG_EXT(U,Yp,skew,Impl::multLinkAdd);
    GENERIC_STENCIL_LEG_EXT(U,Zp,skew,Impl::multLinkAdd);
    GENERIC_STENCIL_LEG_EXT(U,Tp,skew,Impl::multLinkAdd);
    GENERIC_STENCIL_LEG_EXT(U,Xm,skew,Impl::multLinkAdd);
    GENERIC_STENCIL_LEG_EXT(U,Ym,skew,Impl::multLinkAdd);
    GENERIC_STENCIL_LEG_EXT(U,Zm,skew,Impl::multLinkAdd);
    GENERIC_STENCIL_LEG_EXT(U,Tm,skew,Impl::multLinkAdd);
    skew=8;
    GENERIC_STENCIL_LEG_EXT(UUU,Xp,skew,Impl::multLinkAdd);
    GENERIC_STENCIL_LEG_EXT(UUU,Yp,skew,Impl::multLinkAdd);
    GENERIC_STENCIL_LEG_EXT(UUU,Zp,skew,Impl::multLinkAdd);
    GENERIC_STENCIL_LEG_EXT(UUU,Tp,skew,Impl::multLinkAdd);
    GENERIC_STENCIL_LEG_EXT(UUU,Xm,skew,Impl::multLinkAdd);
    GENERIC_STENCIL_LEG_EXT(UUU,Ym,skew,Impl::multLinkAdd);
    GENERIC_STENCIL_LEG_EXT(UUU,Zm,skew,Impl::multLinkAdd);
    GENERIC_STENCIL_LEG_EXT(UUU,Tm,skew,Impl::multLinkAdd);

    if ( nmu ) { 
      if ( dag ) { 
	out._odata[sF] = out._odata[sF] - Uchi;
      } else { 
	out._odata[sF] = out._odata[sF] + Uchi;
      }
    }
  }
};

////////////////////////////////////////////////////////////////////////////////////
// Driving / wrapping routine to select right kernel
////////////////////////////////////////////////////////////////////////////////////

template <class Impl>
void StaggeredKernels<Impl>::DhopSiteDag(StencilImpl &st, LebesgueOrder &lo, DoubledGaugeField &U, DoubledGaugeField &UUU,
					 SiteSpinor *buf, int LLs, int sU,
					 const FermionField &in, FermionField &out,
					 int interior,int exterior)
{
  int dag=1;
  DhopSite(st,lo,U,UUU,buf,LLs,sU,in,out,dag,interior,exterior);
};

template <class Impl>
void StaggeredKernels<Impl>::DhopSite(StencilImpl &st, LebesgueOrder &lo, DoubledGaugeField &U, DoubledGaugeField &UUU,
				      SiteSpinor *buf, int LLs, int sU,
				      const FermionField &in, FermionField &out,
				      int interior,int exterior)
{
  int dag=0;
  DhopSite(st,lo,U,UUU,buf,LLs,sU,in,out,dag,interior,exterior);
};

template <class Impl>
void StaggeredKernels<Impl>::DhopSite(StencilImpl &st, LebesgueOrder &lo, DoubledGaugeField &U, DoubledGaugeField &UUU,
				      SiteSpinor *buf, int LLs,
				      int sU, const FermionField &in, FermionField &out,
				      int dag,int interior,int exterior) 
{
  switch(Opt) {
#ifdef AVX512
  case OptInlineAsm:
    if ( interior && exterior ) {
      DhopSiteAsm(st,lo,U,UUU,buf,LLs,sU,in,out,dag);
    } else { 
      std::cout << GridLogError << "Cannot overlap comms and compute with Staggered assembly"<<std::endl;
      assert(0);
    }
    break;
#endif
  case OptHandUnroll:
    if ( interior && exterior ) {
      DhopSiteHand   (st,lo,U,UUU,buf,LLs,sU,in,out,dag);
    } else if ( interior ) {
      DhopSiteHandInt(st,lo,U,UUU,buf,LLs,sU,in,out,dag);
    } else if ( exterior ) {
      DhopSiteHandExt(st,lo,U,UUU,buf,LLs,sU,in,out,dag);
    }
    break;
  case OptGeneric:
    if ( interior && exterior ) {
      DhopSiteGeneric   (st,lo,U,UUU,buf,LLs,sU,in,out,dag);
    } else if ( interior ) {
      DhopSiteGenericInt(st,lo,U,UUU,buf,LLs,sU,in,out,dag);
    } else if ( exterior ) {
      DhopSiteGenericExt(st,lo,U,UUU,buf,LLs,sU,in,out,dag);
    }
    break;
  default:
    std::cout<<"Oops Opt = "<<Opt<<std::endl;
    assert(0);
    break;
  }
};

template <class Impl>
void StaggeredKernels<Impl>::DhopDir( StencilImpl &st, DoubledGaugeField &U,  DoubledGaugeField &UUU, SiteSpinor *buf, int sF,
				      int sU, const FermionField &in, FermionField &out, int dir, int disp) 
{
  // Disp should be either +1,-1,+3,-3
  // What about "dag" ?
  // Because we work out pU . dS/dU 
  // U
  assert(0);
}

FermOpStaggeredTemplateInstantiate(StaggeredKernels);
FermOpStaggeredVec5dTemplateInstantiate(StaggeredKernels);

}}

