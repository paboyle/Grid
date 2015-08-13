#include <Grid.h>
namespace Grid {
namespace QCD {

template<class Impl> 
void WilsonKernels<Impl>::DiracOptDhopSite(CartesianStencil &st,DoubledGaugeField &U,
						  std::vector<SiteHalfSpinor,alignedAllocator<SiteHalfSpinor> >  &buf,
						  int sF,int sU,const FermionField &in, FermionField &out)
{
  SiteHalfSpinor  tmp;    
  SiteHalfSpinor  chi;    
  SiteHalfSpinor Uchi;
  SiteSpinor result;
  StencilEntry *SE;
  int ptype;

  // Xp
  SE=st.GetEntry(ptype,Xp,sF);
  if ( SE->_is_local && SE->_permute ) {
    spProjXp(tmp,in._odata[SE->_offset]);
    permute(chi,tmp,ptype);
  } else if ( SE->_is_local ) {
    spProjXp(chi,in._odata[SE->_offset]);
  } else { 
    chi=buf[SE->_offset];
  }
  Impl::multLink(Uchi,U._odata[sU],chi,Xp,SE);
  spReconXp(result,Uchi);
    
  // Yp
  SE=st.GetEntry(ptype,Yp,sF);
  if ( SE->_is_local && SE->_permute ) {
    spProjYp(tmp,in._odata[SE->_offset]);
    permute(chi,tmp,ptype);
  } else if ( SE->_is_local ) {
    spProjYp(chi,in._odata[SE->_offset]);
  } else { 
    chi=buf[SE->_offset];
  }
  Impl::multLink(Uchi,U._odata[sU],chi,Yp,SE);
  accumReconYp(result,Uchi);

  // Zp
  SE=st.GetEntry(ptype,Zp,sF);
  if (  SE->_is_local && SE->_permute ) {
    spProjZp(tmp,in._odata[SE->_offset]);
    permute(chi,tmp,ptype);
  } else if ( SE->_is_local ) {
    spProjZp(chi,in._odata[SE->_offset]);
  } else { 
    chi=buf[SE->_offset];
  }
  Impl::multLink(Uchi,U._odata[sU],chi,Zp,SE);
  accumReconZp(result,Uchi);

  // Tp
  SE=st.GetEntry(ptype,Tp,sF);
  if (  SE->_is_local && SE->_permute ) {
    spProjTp(tmp,in._odata[SE->_offset]);
    permute(chi,tmp,ptype);
  } else if ( SE->_is_local ) {
    spProjTp(chi,in._odata[SE->_offset]);
  } else { 
    chi=buf[SE->_offset];
  }
  Impl::multLink(Uchi,U._odata[sU],chi,Tp,SE);
  accumReconTp(result,Uchi);

  // Xm
  SE=st.GetEntry(ptype,Xm,sF);
  if (  SE->_is_local && SE->_permute ) {
    spProjXm(tmp,in._odata[SE->_offset]);
    permute(chi,tmp,ptype);
  } else if ( SE->_is_local ) {
    spProjXm(chi,in._odata[SE->_offset]);
  } else { 
    chi=buf[SE->_offset];
  }
  Impl::multLink(Uchi,U._odata[sU],chi,Xm,SE);
  accumReconXm(result,Uchi);
  
  // Ym
  SE=st.GetEntry(ptype,Ym,sF);
  if (  SE->_is_local && SE->_permute ) {
    spProjYm(tmp,in._odata[SE->_offset]);
    permute(chi,tmp,ptype);
  } else if ( SE->_is_local ) {
    spProjYm(chi,in._odata[SE->_offset]);
  } else { 
    chi=buf[SE->_offset];
  }
  Impl::multLink(Uchi,U._odata[sU],chi,Ym,SE);
  accumReconYm(result,Uchi);
  
  // Zm
  SE=st.GetEntry(ptype,Zm,sF);
  if (  SE->_is_local && SE->_permute ) {
    spProjZm(tmp,in._odata[SE->_offset]);
    permute(chi,tmp,ptype);
  } else if ( SE->_is_local ) {
    spProjZm(chi,in._odata[SE->_offset]);
  } else { 
    chi=buf[SE->_offset];
  }
  Impl::multLink(Uchi,U._odata[sU],chi,Zm,SE);
  accumReconZm(result,Uchi);

  // Tm
  SE=st.GetEntry(ptype,Tm,sF);
  if (  SE->_is_local && SE->_permute ) {
    spProjTm(tmp,in._odata[SE->_offset]);
    permute(chi,tmp,ptype);
  } else if ( SE->_is_local ) {
    spProjTm(chi,in._odata[SE->_offset]);
  } else { 
    chi=buf[SE->_offset];
  }
  Impl::multLink(Uchi,U._odata[sU],chi,Tm,SE);
  accumReconTm(result,Uchi);

  vstream(out._odata[sF],result*(-0.5));
};

template<class Impl> 
void WilsonKernels<Impl>::DiracOptDhopSiteDag(CartesianStencil &st,DoubledGaugeField &U,
					      std::vector<SiteHalfSpinor,alignedAllocator<SiteHalfSpinor> >  &buf,
					      int sF,int sU,const FermionField &in, FermionField &out)
{
  SiteHalfSpinor  tmp;    
  SiteHalfSpinor  chi;    
  SiteSpinor result;
  SiteHalfSpinor Uchi;
  StencilEntry *SE;
  int ptype;

  // Xp
  SE=st.GetEntry(ptype,Xm,sF);
  if (  SE->_is_local && SE->_permute ) {
    spProjXp(tmp,in._odata[SE->_offset]);
    permute(chi,tmp,ptype);
  } else if ( SE->_is_local ) {
    spProjXp(chi,in._odata[SE->_offset]);
  } else { 
    chi=buf[SE->_offset];
  }
  Impl::multLink(Uchi,U._odata[sU],chi,Xm,SE);
  spReconXp(result,Uchi);

  // Yp
  SE=st.GetEntry(ptype,Ym,sF);
  if (  SE->_is_local && SE->_permute ) {
    spProjYp(tmp,in._odata[SE->_offset]);
    permute(chi,tmp,ptype);
  } else if ( SE->_is_local ) {
    spProjYp(chi,in._odata[SE->_offset]);
  } else { 
    chi=buf[SE->_offset];
  }
  Impl::multLink(Uchi,U._odata[sU],chi,Ym,SE);
  accumReconYp(result,Uchi);
  
  // Zp
  SE=st.GetEntry(ptype,Zm,sF);
  if (  SE->_is_local && SE->_permute ) {
    spProjZp(tmp,in._odata[SE->_offset]);
    permute(chi,tmp,ptype);
  } else if ( SE->_is_local ) {
    spProjZp(chi,in._odata[SE->_offset]);
  } else { 
    chi=buf[SE->_offset];
  }
  Impl::multLink(Uchi,U._odata[sU],chi,Zm,SE);
  accumReconZp(result,Uchi);
  
  // Tp
  SE=st.GetEntry(ptype,Tm,sF);
  if (  SE->_is_local && SE->_permute ) {
    spProjTp(tmp,in._odata[SE->_offset]);
    permute(chi,tmp,ptype);
  } else if ( SE->_is_local ) {
    spProjTp(chi,in._odata[SE->_offset]);
  } else { 
    chi=buf[SE->_offset];
  }
  Impl::multLink(Uchi,U._odata[sU],chi,Tm,SE);
  accumReconTp(result,Uchi);
  
  // Xm
  SE=st.GetEntry(ptype,Xp,sF);
  if (  SE->_is_local && SE->_permute ) {
    spProjXm(tmp,in._odata[SE->_offset]);
    permute(chi,tmp,ptype);
  } else if ( SE->_is_local ) {
    spProjXm(chi,in._odata[SE->_offset]);
  } else { 
    chi=buf[SE->_offset];
  }
  Impl::multLink(Uchi,U._odata[sU],chi,Xp,SE);
  accumReconXm(result,Uchi);

  // Ym
  SE=st.GetEntry(ptype,Yp,sF);
  if (  SE->_is_local && SE->_permute ) {
    spProjYm(tmp,in._odata[SE->_offset]);
    permute(chi,tmp,ptype);
  } else if ( SE->_is_local ) {
    spProjYm(chi,in._odata[SE->_offset]);
  } else { 
    chi=buf[SE->_offset];
  }
  Impl::multLink(Uchi,U._odata[sU],chi,Yp,SE);
  accumReconYm(result,Uchi);

  // Zm
  SE=st.GetEntry(ptype,Zp,sF);
  if (  SE->_is_local && SE->_permute ) {
    spProjZm(tmp,in._odata[SE->_offset]);
    permute(chi,tmp,ptype);
  } else if ( SE->_is_local ) {
    spProjZm(chi,in._odata[SE->_offset]);
  } else { 
    chi=buf[SE->_offset];
  }
  Impl::multLink(Uchi,U._odata[sU],chi,Zp,SE);
  accumReconZm(result,Uchi);
    
  // Tm
  SE=st.GetEntry(ptype,Tp,sF);
  if (  SE->_is_local && SE->_permute ) {
    spProjTm(tmp,in._odata[SE->_offset]);
    permute(chi,tmp,ptype);
  } else if ( SE->_is_local ) {
    spProjTm(chi,in._odata[SE->_offset]);
  } else { 
    chi=buf[SE->_offset];
  }
  Impl::multLink(Uchi,U._odata[sU],chi,Tp,SE);
  accumReconTm(result,Uchi);
  
  vstream(out._odata[sF],result*(-0.5));
}

template<class Impl> 
void WilsonKernels<Impl>::DiracOptDhopDir(CartesianStencil &st,DoubledGaugeField &U,
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
    Impl::multLink(Uchi,U._odata[sU],chi,dir,SE);
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
    Impl::multLink(Uchi,U._odata[sU],chi,dir,SE);
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
    Impl::multLink(Uchi,U._odata[sU],chi,dir,SE);
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
    Impl::multLink(Uchi,U._odata[sU],chi,dir,SE);
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
    Impl::multLink(Uchi,U._odata[sU],chi,dir,SE);
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
    Impl::multLink(Uchi,U._odata[sU],chi,dir,SE);
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
    Impl::multLink(Uchi,U._odata[sU],chi,dir,SE);
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
    Impl::multLink(Uchi,U._odata[sU],chi,dir,SE);
    spReconTm(result,Uchi);
  }

  vstream(out._odata[sF],result*(-0.5));
}

  FermOpTemplateInstantiate(WilsonKernels);

}}
