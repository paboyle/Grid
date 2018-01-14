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

NAMESPACE_BEGIN(Grid);

int StaggeredKernelsStatic::Opt= StaggeredKernelsStatic::OptGeneric;

template <class Impl>
StaggeredKernels<Impl>::StaggeredKernels(const ImplParams &p) : Base(p){};

////////////////////////////////////////////
// Generic implementation; move to different file?
////////////////////////////////////////////

template <class Impl>
void StaggeredKernels<Impl>::DhopSiteDepth(StencilImpl &st, LebesgueOrder &lo, DoubledGaugeField &U,
					   SiteSpinor *buf, int sF,
					   int sU, const FermionField &in, SiteSpinor &out,int threeLink) {
  const SiteSpinor *chi_p;
  SiteSpinor chi;
  SiteSpinor Uchi;
  StencilEntry *SE;
  int ptype;
  int skew = 0;
  if (threeLink) skew=8;
  ///////////////////////////
  // Xp
  ///////////////////////////

  SE = st.GetEntry(ptype, Xp+skew, sF);
  if (SE->_is_local) {
    if (SE->_permute) {
      chi_p = &chi;
      permute(chi,  in._odata[SE->_offset], ptype);
    } else {
      chi_p = &in._odata[SE->_offset];
    }
  } else {
    chi_p = &buf[SE->_offset];
  }
  Impl::multLink(Uchi, U._odata[sU], *chi_p, Xp);

  ///////////////////////////
  // Yp
  ///////////////////////////
  SE = st.GetEntry(ptype, Yp+skew, sF);
  if (SE->_is_local) {
    if (SE->_permute) {
      chi_p = &chi;
      permute(chi,  in._odata[SE->_offset], ptype);
    } else {
      chi_p = &in._odata[SE->_offset];
    }
  } else {
    chi_p = &buf[SE->_offset];
  }
  Impl::multLinkAdd(Uchi, U._odata[sU], *chi_p, Yp);

  ///////////////////////////
  // Zp
  ///////////////////////////
  SE = st.GetEntry(ptype, Zp+skew, sF);
  if (SE->_is_local) {
    if (SE->_permute) {
      chi_p = &chi;
      permute(chi,  in._odata[SE->_offset], ptype);
    } else {
      chi_p = &in._odata[SE->_offset];
    }
  } else {
    chi_p = &buf[SE->_offset];
  }
  Impl::multLinkAdd(Uchi, U._odata[sU], *chi_p, Zp);

  ///////////////////////////
  // Tp
  ///////////////////////////
  SE = st.GetEntry(ptype, Tp+skew, sF);
  if (SE->_is_local) {
    if (SE->_permute) {
      chi_p = &chi;
      permute(chi,  in._odata[SE->_offset], ptype);
    } else {
      chi_p = &in._odata[SE->_offset];
    }
  } else {
    chi_p = &buf[SE->_offset];
  }
  Impl::multLinkAdd(Uchi, U._odata[sU], *chi_p, Tp);

  ///////////////////////////
  // Xm
  ///////////////////////////
  SE = st.GetEntry(ptype, Xm+skew, sF);
  if (SE->_is_local) {
    if (SE->_permute) {
      chi_p = &chi;
      permute(chi,  in._odata[SE->_offset], ptype);
    } else {
      chi_p = &in._odata[SE->_offset];
    }
  } else {
    chi_p = &buf[SE->_offset];
  }
  Impl::multLinkAdd(Uchi, U._odata[sU], *chi_p, Xm);

  ///////////////////////////
  // Ym
  ///////////////////////////
  SE = st.GetEntry(ptype, Ym+skew, sF);
  if (SE->_is_local) {
    if (SE->_permute) {
      chi_p = &chi;
      permute(chi,  in._odata[SE->_offset], ptype);
    } else {
      chi_p = &in._odata[SE->_offset];
    }
  } else {
    chi_p = &buf[SE->_offset];
  }
  Impl::multLinkAdd(Uchi, U._odata[sU], *chi_p, Ym);

  ///////////////////////////
  // Zm
  ///////////////////////////
  SE = st.GetEntry(ptype, Zm+skew, sF);
  if (SE->_is_local) {
    if (SE->_permute) {
      chi_p = &chi;
      permute(chi,  in._odata[SE->_offset], ptype);
    } else {
      chi_p = &in._odata[SE->_offset];
    }
  } else {
    chi_p = &buf[SE->_offset];
  }
  Impl::multLinkAdd(Uchi, U._odata[sU], *chi_p, Zm);

  ///////////////////////////
  // Tm
  ///////////////////////////
  SE = st.GetEntry(ptype, Tm+skew, sF);
  if (SE->_is_local) {
    if (SE->_permute) {
      chi_p = &chi;
      permute(chi,  in._odata[SE->_offset], ptype);
    } else {
      chi_p = &in._odata[SE->_offset];
    }
  } else {
    chi_p = &buf[SE->_offset];
  }
  Impl::multLinkAdd(Uchi, U._odata[sU], *chi_p, Tm);

  vstream(out, Uchi);
};

template <class Impl>
void StaggeredKernels<Impl>::DhopSiteDag(StencilImpl &st, LebesgueOrder &lo, DoubledGaugeField &U, DoubledGaugeField &UUU,
					 SiteSpinor *buf, int LLs, int sU,
					 const FermionField &in, FermionField &out) {
  SiteSpinor naik;
  SiteSpinor naive;
  int oneLink  =0;
  int threeLink=1;
  int dag=1;
  switch(Opt) {
#ifdef AVX512
    //FIXME; move the sign into the Asm routine
  case OptInlineAsm:
    DhopSiteAsm(st,lo,U,UUU,buf,LLs,sU,in,out);
    for(int s=0;s<LLs;s++) {
      int sF=s+LLs*sU;
      out._odata[sF]=-out._odata[sF];
    }
    break;
#endif
  case OptHandUnroll:
    DhopSiteHand(st,lo,U,UUU,buf,LLs,sU,in,out,dag);
    break;
  case OptGeneric:
    for(int s=0;s<LLs;s++){
      int sF=s+LLs*sU;
      DhopSiteDepth(st,lo,U,buf,sF,sU,in,naive,oneLink);
      DhopSiteDepth(st,lo,UUU,buf,sF,sU,in,naik,threeLink);
      out._odata[sF] =-naive-naik; 
    }
    break;
  default:
    std::cout<<"Oops Opt = "<<Opt<<std::endl;
    assert(0);
    break;
  }
};

template <class Impl>
void StaggeredKernels<Impl>::DhopSite(StencilImpl &st, LebesgueOrder &lo, DoubledGaugeField &U, DoubledGaugeField &UUU,
				      SiteSpinor *buf, int LLs,
				      int sU, const FermionField &in, FermionField &out) 
{
  int oneLink  =0;
  int threeLink=1;
  SiteSpinor naik;
  SiteSpinor naive;
  int dag=0;
  switch(Opt) {
#ifdef AVX512
  case OptInlineAsm:
    DhopSiteAsm(st,lo,U,UUU,buf,LLs,sU,in,out);
    break;
#endif
  case OptHandUnroll:
    DhopSiteHand(st,lo,U,UUU,buf,LLs,sU,in,out,dag);
    break;
  case OptGeneric:
    for(int s=0;s<LLs;s++){
      int sF=LLs*sU+s;
      //      assert(sF<in._odata.size());
      //      assert(sU< U._odata.size());
      //      assert(sF>=0);      assert(sU>=0);
      DhopSiteDepth(st,lo,U,buf,sF,sU,in,naive,oneLink);
      DhopSiteDepth(st,lo,UUU,buf,sF,sU,in,naik,threeLink);
      out._odata[sF] =naive+naik;
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

NAMESPACE_END(Grid);

