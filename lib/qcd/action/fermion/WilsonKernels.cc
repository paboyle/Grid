    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/action/fermion/WilsonKernels.cc

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: Peter Boyle <peterboyle@Peters-MacBook-Pro-2.local>
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
namespace Grid {
namespace QCD {

template<class Impl> 
WilsonKernels<Impl>::WilsonKernels(const ImplParams &p): Base(p) {};

  // Need controls to do interior, exterior, or both
template<class Impl> 
void WilsonKernels<Impl>::DiracOptDhopSiteDag(StencilImpl &st,DoubledGaugeField &U,
					   std::vector<SiteHalfSpinor,alignedAllocator<SiteHalfSpinor> >  &buf,
					   int sF,int sU,const FermionField &in, FermionField &out)
{
  SiteHalfSpinor  tmp;    
  SiteHalfSpinor  chi;    
  SiteHalfSpinor *chi_p;
  SiteHalfSpinor Uchi;
  SiteSpinor result;
  StencilEntry *SE;
  int ptype;

  ///////////////////////////
  // Xp
  ///////////////////////////
  SE=st.GetEntry(ptype,Xp,sF);

  if (SE->_is_local ) { 
    chi_p = &chi;
    if ( SE->_permute ) {
      spProjXp(tmp,in._odata[SE->_offset]);
      permute(chi,tmp,ptype);
    } else {
      spProjXp(chi,in._odata[SE->_offset]);
    }
  } else { 
    chi_p=&buf[SE->_offset];
  }
  
  Impl::multLink(Uchi,U._odata[sU],*chi_p,Xp,SE,st);
  spReconXp(result,Uchi);
    
  ///////////////////////////
  // Yp
  ///////////////////////////
  SE=st.GetEntry(ptype,Yp,sF);

  if ( SE->_is_local ) { 
    chi_p = &chi;
    if ( SE->_permute ) {
      spProjYp(tmp,in._odata[SE->_offset]);
      permute(chi,tmp,ptype);
    } else {
      spProjYp(chi,in._odata[SE->_offset]);
    }
  } else { 
    chi_p=&buf[SE->_offset];
  }

  Impl::multLink(Uchi,U._odata[sU],*chi_p,Yp,SE,st);
  accumReconYp(result,Uchi);

  ///////////////////////////
  // Zp
  ///////////////////////////
  SE=st.GetEntry(ptype,Zp,sF);

  if ( SE->_is_local ) { 
    chi_p = &chi;
    if ( SE->_permute ) {
      spProjZp(tmp,in._odata[SE->_offset]);
      permute(chi,tmp,ptype);
    } else {
      spProjZp(chi,in._odata[SE->_offset]);
    }
  } else { 
    chi_p=&buf[SE->_offset];
  }

  Impl::multLink(Uchi,U._odata[sU],*chi_p,Zp,SE,st);
  accumReconZp(result,Uchi);

  ///////////////////////////
  // Tp
  ///////////////////////////
  SE=st.GetEntry(ptype,Tp,sF);

  if ( SE->_is_local ) {
    chi_p = &chi;
    if ( SE->_permute ) {
      spProjTp(tmp,in._odata[SE->_offset]);
      permute(chi,tmp,ptype);
    } else {
      spProjTp(chi,in._odata[SE->_offset]);
    }
  } else {
    chi_p=&buf[SE->_offset];
  }

  Impl::multLink(Uchi,U._odata[sU],*chi_p,Tp,SE,st);
  accumReconTp(result,Uchi);

  ///////////////////////////
  // Xm
  ///////////////////////////
  SE=st.GetEntry(ptype,Xm,sF);

  if ( SE->_is_local ) {
    chi_p = &chi;
    if ( SE->_permute ) {
      spProjXm(tmp,in._odata[SE->_offset]);
      permute(chi,tmp,ptype);
    } else {
      spProjXm(chi,in._odata[SE->_offset]);
    }
  } else {
    chi_p=&buf[SE->_offset];
  }

  Impl::multLink(Uchi,U._odata[sU],*chi_p,Xm,SE,st);
  accumReconXm(result,Uchi);
  
  ///////////////////////////
  // Ym
  ///////////////////////////
  SE=st.GetEntry(ptype,Ym,sF);

  if ( SE->_is_local ) {
    chi_p = &chi;
    if ( SE->_permute ) {
      spProjYm(tmp,in._odata[SE->_offset]);
      permute(chi,tmp,ptype);
    } else {
      spProjYm(chi,in._odata[SE->_offset]);
    }
  } else {
    chi_p=&buf[SE->_offset];
  }

  Impl::multLink(Uchi,U._odata[sU],*chi_p,Ym,SE,st);
  accumReconYm(result,Uchi);
  
  ///////////////////////////
  // Zm
  ///////////////////////////
  SE=st.GetEntry(ptype,Zm,sF);

  if ( SE->_is_local ) {
    chi_p = &chi;
    if ( SE->_permute ) {
      spProjZm(tmp,in._odata[SE->_offset]);
      permute(chi,tmp,ptype);
    } else {
      spProjZm(chi,in._odata[SE->_offset]);
    }
  } else {
    chi_p=&buf[SE->_offset];
  }

  Impl::multLink(Uchi,U._odata[sU],*chi_p,Zm,SE,st);
  accumReconZm(result,Uchi);

  ///////////////////////////
  // Tm
  ///////////////////////////
  SE=st.GetEntry(ptype,Tm,sF);

  if ( SE->_is_local ) {
    chi_p = &chi;
    if ( SE->_permute ) {
      spProjTm(tmp,in._odata[SE->_offset]);
      permute(chi,tmp,ptype);
    } else { 
      spProjTm(chi,in._odata[SE->_offset]);
    }
  } else {
    chi_p=&buf[SE->_offset];
  }

  Impl::multLink(Uchi,U._odata[sU],*chi_p,Tm,SE,st);
  accumReconTm(result,Uchi);

  vstream(out._odata[sF],result);
};


  // Need controls to do interior, exterior, or both
template<class Impl> 
void WilsonKernels<Impl>::DiracOptDhopSite(StencilImpl &st,DoubledGaugeField &U,
					   std::vector<SiteHalfSpinor,alignedAllocator<SiteHalfSpinor> >  &buf,
					   int sF,int sU,const FermionField &in, FermionField &out)
{
  SiteHalfSpinor  tmp;    
  SiteHalfSpinor  chi;    
  SiteHalfSpinor *chi_p;    
  SiteHalfSpinor Uchi;
  SiteSpinor result;
  StencilEntry *SE;
  int ptype;

  ///////////////////////////
  // Xp
  ///////////////////////////
  SE=st.GetEntry(ptype,Xm,sF);

  if ( SE->_is_local ) { 
    chi_p = &chi;
    if ( SE->_permute ) {
      spProjXp(tmp,in._odata[SE->_offset]);
      permute(chi,tmp,ptype);
    } else {
      spProjXp(chi,in._odata[SE->_offset]);
    }
  } else { 
    chi_p=&buf[SE->_offset];
  }
  
  Impl::multLink(Uchi,U._odata[sU],*chi_p,Xm,SE,st);
  spReconXp(result,Uchi);
    
  ///////////////////////////
  // Yp
  ///////////////////////////
  SE=st.GetEntry(ptype,Ym,sF);

  if ( SE->_is_local ) { 
    chi_p = &chi;
    if ( SE->_permute ) {
      spProjYp(tmp,in._odata[SE->_offset]);
      permute(chi,tmp,ptype);
    } else {
      spProjYp(chi,in._odata[SE->_offset]);
    }
  } else { 
    chi_p=&buf[SE->_offset];
  }

  Impl::multLink(Uchi,U._odata[sU],*chi_p,Ym,SE,st);
  accumReconYp(result,Uchi);

  ///////////////////////////
  // Zp
  ///////////////////////////
  SE=st.GetEntry(ptype,Zm,sF);

  if ( SE->_is_local ) { 
    chi_p = &chi;
    if ( SE->_permute ) {
      spProjZp(tmp,in._odata[SE->_offset]);
      permute(chi,tmp,ptype);
    } else {
      spProjZp(chi,in._odata[SE->_offset]);
    }
  } else { 
    chi_p=&buf[SE->_offset];
  }

  Impl::multLink(Uchi,U._odata[sU],*chi_p,Zm,SE,st);
  accumReconZp(result,Uchi);

  ///////////////////////////
  // Tp
  ///////////////////////////
  SE=st.GetEntry(ptype,Tm,sF);

  if ( SE->_is_local ) {
    chi_p = &chi;
    if ( SE->_permute ) {
      spProjTp(tmp,in._odata[SE->_offset]);
      permute(chi,tmp,ptype);
    } else {
      spProjTp(chi,in._odata[SE->_offset]);
    }
  } else {
    chi_p=&buf[SE->_offset];
  }

  Impl::multLink(Uchi,U._odata[sU],*chi_p,Tm,SE,st);
  accumReconTp(result,Uchi);

  ///////////////////////////
  // Xm
  ///////////////////////////
  SE=st.GetEntry(ptype,Xp,sF);

  if ( SE->_is_local ) {
    chi_p = &chi;
    if ( SE->_permute ) {
      spProjXm(tmp,in._odata[SE->_offset]);
      permute(chi,tmp,ptype);
    } else {
      spProjXm(chi,in._odata[SE->_offset]);
    }
  } else {
    chi_p=&buf[SE->_offset];
  }

  Impl::multLink(Uchi,U._odata[sU],*chi_p,Xp,SE,st);
  accumReconXm(result,Uchi);

  ///////////////////////////
  // Ym
  ///////////////////////////
  SE=st.GetEntry(ptype,Yp,sF);

  if ( SE->_is_local ) {
    chi_p = &chi;
    if ( SE->_permute ) {
      spProjYm(tmp,in._odata[SE->_offset]);
      permute(chi,tmp,ptype);
    } else {
      spProjYm(chi,in._odata[SE->_offset]);
    }
  } else {
    chi_p=&buf[SE->_offset];
  }

  Impl::multLink(Uchi,U._odata[sU],*chi_p,Yp,SE,st);
  accumReconYm(result,Uchi);
  
  ///////////////////////////
  // Zm
  ///////////////////////////
  SE=st.GetEntry(ptype,Zp,sF);

  if ( SE->_is_local ) {
    chi_p = &chi;
    if ( SE->_permute ) {
      spProjZm(tmp,in._odata[SE->_offset]);
      permute(chi,tmp,ptype);
    } else {
      spProjZm(chi,in._odata[SE->_offset]);
    }
  } else {
    chi_p=&buf[SE->_offset];
  }

  Impl::multLink(Uchi,U._odata[sU],*chi_p,Zp,SE,st);
  accumReconZm(result,Uchi);

  ///////////////////////////
  // Tm
  ///////////////////////////
  SE=st.GetEntry(ptype,Tp,sF);

  if ( SE->_is_local ) {
    chi_p = &chi;
    if ( SE->_permute ) {
      spProjTm(tmp,in._odata[SE->_offset]);
      permute(chi,tmp,ptype);
    } else { 
      spProjTm(chi,in._odata[SE->_offset]);
    }
  } else {
    chi_p=&buf[SE->_offset];
  }

  Impl::multLink(Uchi,U._odata[sU],*chi_p,Tp,SE,st);
  accumReconTm(result,Uchi);

  vstream(out._odata[sF],result);
};

template<class Impl> 
void WilsonKernels<Impl>::DiracOptDhopDir(StencilImpl &st,DoubledGaugeField &U,
					  std::vector<SiteHalfSpinor,alignedAllocator<SiteHalfSpinor> >  &buf,
					  int sF,int sU,const FermionField &in, FermionField &out,int dir,int gamma)
{
  SiteHalfSpinor  tmp;    
  SiteHalfSpinor  chi;    
  SiteSpinor   result;
  SiteHalfSpinor Uchi;
  StencilEntry *SE;
  int ptype;

  SE=st.GetEntry(ptype,dir,sF);

  // Xp
  if(gamma==Xp){
    if (  SE->_is_local && SE->_permute ) {
      spProjXp(tmp,in._odata[SE->_offset]);
      permute(chi,tmp,ptype);
    } else if ( SE->_is_local ) {
      spProjXp(chi,in._odata[SE->_offset]);
    } else { 
      chi=buf[SE->_offset];
    }
    Impl::multLink(Uchi,U._odata[sU],chi,dir,SE,st);
    spReconXp(result,Uchi);
  }

  // Yp
  if ( gamma==Yp ){
    if (  SE->_is_local && SE->_permute ) {
      spProjYp(tmp,in._odata[SE->_offset]);
      permute(chi,tmp,ptype);
    } else if ( SE->_is_local ) {
      spProjYp(chi,in._odata[SE->_offset]);
    } else { 
      chi=buf[SE->_offset];
    }
    Impl::multLink(Uchi,U._odata[sU],chi,dir,SE,st);
    spReconYp(result,Uchi);
  }
  
  // Zp
  if ( gamma ==Zp ){
    if (  SE->_is_local && SE->_permute ) {
      spProjZp(tmp,in._odata[SE->_offset]);
      permute(chi,tmp,ptype);
    } else if ( SE->_is_local ) {
      spProjZp(chi,in._odata[SE->_offset]);
    } else { 
      chi=buf[SE->_offset];
    }
    Impl::multLink(Uchi,U._odata[sU],chi,dir,SE,st);
    spReconZp(result,Uchi);
  }
  
  // Tp
  if ( gamma ==Tp ){
    if (  SE->_is_local && SE->_permute ) {
      spProjTp(tmp,in._odata[SE->_offset]);
      permute(chi,tmp,ptype);
    } else if ( SE->_is_local ) {
      spProjTp(chi,in._odata[SE->_offset]);
    } else { 
      chi=buf[SE->_offset];
    }
    Impl::multLink(Uchi,U._odata[sU],chi,dir,SE,st);
    spReconTp(result,Uchi);
  }

  // Xm
  if ( gamma==Xm ){
    if (  SE->_is_local && SE->_permute ) {
      spProjXm(tmp,in._odata[SE->_offset]);
      permute(chi,tmp,ptype);
    } else if ( SE->_is_local ) {
      spProjXm(chi,in._odata[SE->_offset]);
    } else { 
      chi=buf[SE->_offset];
    }
    Impl::multLink(Uchi,U._odata[sU],chi,dir,SE,st);
    spReconXm(result,Uchi);
  }

  // Ym
  if ( gamma == Ym ){
    if (  SE->_is_local && SE->_permute ) {
      spProjYm(tmp,in._odata[SE->_offset]);
      permute(chi,tmp,ptype);
    } else if ( SE->_is_local ) {
      spProjYm(chi,in._odata[SE->_offset]);
    } else { 
      chi=buf[SE->_offset];
    }
    Impl::multLink(Uchi,U._odata[sU],chi,dir,SE,st);
    spReconYm(result,Uchi);
  }

  // Zm
  if ( gamma == Zm ){
    if (  SE->_is_local && SE->_permute ) {
      spProjZm(tmp,in._odata[SE->_offset]);
      permute(chi,tmp,ptype);
    } else if ( SE->_is_local ) {
      spProjZm(chi,in._odata[SE->_offset]);
    } else { 
      chi=buf[SE->_offset];
    }
    Impl::multLink(Uchi,U._odata[sU],chi,dir,SE,st);
    spReconZm(result,Uchi);
  }
  
  // Tm
  if ( gamma==Tm ) {
    if (  SE->_is_local && SE->_permute ) {
      spProjTm(tmp,in._odata[SE->_offset]);
      permute(chi,tmp,ptype);
    } else if ( SE->_is_local ) {
      spProjTm(chi,in._odata[SE->_offset]);
    } else { 
      chi=buf[SE->_offset];
    }
    Impl::multLink(Uchi,U._odata[sU],chi,dir,SE,st);
    spReconTm(result,Uchi);
  }

  vstream(out._odata[sF],result);
}

#if ( ! defined(AVX512) )
template<class Impl> 
void WilsonKernels<Impl>::DiracOptAsmDhopSite(StencilImpl &st,DoubledGaugeField &U,
					      std::vector<SiteHalfSpinor,alignedAllocator<SiteHalfSpinor> >  &buf,
					      int sF,int sU,const FermionField &in, FermionField &out)
{
  DiracOptDhopSite(st,U,buf,sF,sU,in,out); // will template override for Wilson Nc=3
}
#endif

  FermOpTemplateInstantiate(WilsonKernels);
template class WilsonKernels<DomainWallRedBlack5dImplF>;		
template class WilsonKernels<DomainWallRedBlack5dImplD>;

}}
