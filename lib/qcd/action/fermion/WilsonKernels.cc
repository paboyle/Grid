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
  SiteSpinor result;
  SiteHalfSpinor Uchi;
  int offset,local,perm, ptype;

  // Xp
  int ss = sF;
  offset = st._offsets [Xp][ss];
  local  = st._is_local[Xp][ss];
  perm   = st._permute[Xp][ss];
  
  ptype  = st._permute_type[Xp];
  if ( local && perm ) {
    spProjXp(tmp,in._odata[offset]);
    permute(chi,tmp,ptype);
  } else if ( local ) {
    spProjXp(chi,in._odata[offset]);
  } else { 
    chi=buf[offset];
  }
  Impl::multLink(Uchi,U._odata[sU],chi,Xp);
  spReconXp(result,Uchi);
    
  // Yp
  offset = st._offsets [Yp][ss];
  local  = st._is_local[Yp][ss];
  perm   = st._permute[Yp][ss];
  ptype  = st._permute_type[Yp];
  if ( local && perm ) {
    spProjYp(tmp,in._odata[offset]);
    permute(chi,tmp,ptype);
  } else if ( local ) {
    spProjYp(chi,in._odata[offset]);
  } else { 
    chi=buf[offset];
  }
  Impl::multLink(Uchi,U._odata[sU],chi,Yp);
  accumReconYp(result,Uchi);

  // Zp
  offset = st._offsets [Zp][ss];
  local  = st._is_local[Zp][ss];
  perm   = st._permute[Zp][ss];
  ptype  = st._permute_type[Zp];
  if ( local && perm ) {
    spProjZp(tmp,in._odata[offset]);
    permute(chi,tmp,ptype);
  } else if ( local ) {
    spProjZp(chi,in._odata[offset]);
  } else { 
    chi=buf[offset];
  }
  Impl::multLink(Uchi,U._odata[sU],chi,Zp);
  accumReconZp(result,Uchi);

  // Tp
  offset = st._offsets [Tp][ss];
  local  = st._is_local[Tp][ss];
  perm   = st._permute[Tp][ss];
  ptype  = st._permute_type[Tp];
  if ( local && perm ) {
    spProjTp(tmp,in._odata[offset]);
    permute(chi,tmp,ptype);
  } else if ( local ) {
    spProjTp(chi,in._odata[offset]);
  } else { 
    chi=buf[offset];
  }
  Impl::multLink(Uchi,U._odata[sU],chi,Tp);
  accumReconTp(result,Uchi);

  // Xm
  offset = st._offsets [Xm][ss];
  local  = st._is_local[Xm][ss];
  perm   = st._permute[Xm][ss];
  ptype  = st._permute_type[Xm];
  
  if ( local && perm ) {
    spProjXm(tmp,in._odata[offset]);
    permute(chi,tmp,ptype);
  } else if ( local ) {
    spProjXm(chi,in._odata[offset]);
  } else { 
    chi=buf[offset];
  }
  Impl::multLink(Uchi,U._odata[sU],chi,Xm);
  accumReconXm(result,Uchi);
  
  // Ym
  offset = st._offsets [Ym][ss];
  local  = st._is_local[Ym][ss];
  perm   = st._permute[Ym][ss];
  ptype  = st._permute_type[Ym];
  
  if ( local && perm ) {
    spProjYm(tmp,in._odata[offset]);
    permute(chi,tmp,ptype);
  } else if ( local ) {
    spProjYm(chi,in._odata[offset]);
  } else { 
    chi=buf[offset];
  }
  Impl::multLink(Uchi,U._odata[sU],chi,Ym);
  accumReconYm(result,Uchi);
  
  // Zm
  offset = st._offsets [Zm][ss];
  local  = st._is_local[Zm][ss];
  perm   = st._permute[Zm][ss];
  ptype  = st._permute_type[Zm];
  if ( local && perm ) {
    spProjZm(tmp,in._odata[offset]);
    permute(chi,tmp,ptype);
  } else if ( local ) {
    spProjZm(chi,in._odata[offset]);
  } else { 
    chi=buf[offset];
  }
  Impl::multLink(Uchi,U._odata[sU],chi,Zm);
  accumReconZm(result,Uchi);

  // Tm
  offset = st._offsets [Tm][ss];
  local  = st._is_local[Tm][ss];
  perm   = st._permute[Tm][ss];
  ptype  = st._permute_type[Tm];
  if ( local && perm ) {
    spProjTm(tmp,in._odata[offset]);
    permute(chi,tmp,ptype);
  } else if ( local ) {
    spProjTm(chi,in._odata[offset]);
  } else { 
    chi=buf[offset];
  }
  Impl::multLink(Uchi,U._odata[sU],chi,Tm);
  accumReconTm(result,Uchi);

  vstream(out._odata[ss],result*(-0.5));
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
  int offset,local,perm, ptype;

  // Xp
  int ss=sF;
  offset = st._offsets [Xm][ss];
  local  = st._is_local[Xm][ss];
  perm   = st._permute[Xm][ss];

  ptype  = st._permute_type[Xm];
  if ( local && perm ) {
    spProjXp(tmp,in._odata[offset]);
    permute(chi,tmp,ptype);
  } else if ( local ) {
    spProjXp(chi,in._odata[offset]);
  } else { 
    chi=buf[offset];
  }
  Impl::multLink(Uchi,U._odata[sU],chi,Xm);
  spReconXp(result,Uchi);

  // Yp
  offset = st._offsets [Ym][ss];
  local  = st._is_local[Ym][ss];
  perm   = st._permute[Ym][ss];
  ptype  = st._permute_type[Ym];
  if ( local && perm ) {
    spProjYp(tmp,in._odata[offset]);
    permute(chi,tmp,ptype);
  } else if ( local ) {
    spProjYp(chi,in._odata[offset]);
  } else { 
    chi=buf[offset];
  }
  Impl::multLink(Uchi,U._odata[sU],chi,Ym);
  accumReconYp(result,Uchi);
  
  // Zp
  offset = st._offsets [Zm][ss];
  local  = st._is_local[Zm][ss];
  perm   = st._permute[Zm][ss];
  ptype  = st._permute_type[Zm];
  if ( local && perm ) {
    spProjZp(tmp,in._odata[offset]);
    permute(chi,tmp,ptype);
  } else if ( local ) {
    spProjZp(chi,in._odata[offset]);
  } else { 
    chi=buf[offset];
  }
  Impl::multLink(Uchi,U._odata[sU],chi,Zm);
  accumReconZp(result,Uchi);
  
  // Tp
  offset = st._offsets [Tm][ss];
  local  = st._is_local[Tm][ss];
  perm   = st._permute[Tm][ss];
  ptype  = st._permute_type[Tm];
  if ( local && perm ) {
    spProjTp(tmp,in._odata[offset]);
    permute(chi,tmp,ptype);
  } else if ( local ) {
    spProjTp(chi,in._odata[offset]);
  } else { 
    chi=buf[offset];
  }
  Impl::multLink(Uchi,U._odata[sU],chi,Tm);
  accumReconTp(result,Uchi);
  
  // Xm
  offset = st._offsets [Xp][ss];
  local  = st._is_local[Xp][ss];
  perm   = st._permute[Xp][ss];
  ptype  = st._permute_type[Xp];

  if ( local && perm ) {
    spProjXm(tmp,in._odata[offset]);
    permute(chi,tmp,ptype);
  } else if ( local ) {
    spProjXm(chi,in._odata[offset]);
  } else { 
    chi=buf[offset];
  }
  Impl::multLink(Uchi,U._odata[sU],chi,Xp);
  accumReconXm(result,Uchi);

  // Ym
  offset = st._offsets [Yp][ss];
  local  = st._is_local[Yp][ss];
  perm   = st._permute[Yp][ss];
  ptype  = st._permute_type[Yp];

  if ( local && perm ) {
    spProjYm(tmp,in._odata[offset]);
    permute(chi,tmp,ptype);
  } else if ( local ) {
    spProjYm(chi,in._odata[offset]);
  } else { 
    chi=buf[offset];
  }
  Impl::multLink(Uchi,U._odata[sU],chi,Yp);
  accumReconYm(result,Uchi);

  // Zm
  offset = st._offsets [Zp][ss];
  local  = st._is_local[Zp][ss];
  perm   = st._permute[Zp][ss];
  ptype  = st._permute_type[Zp];
  if ( local && perm ) {
    spProjZm(tmp,in._odata[offset]);
    permute(chi,tmp,ptype);
  } else if ( local ) {
    spProjZm(chi,in._odata[offset]);
  } else { 
    chi=buf[offset];
  }
  Impl::multLink(Uchi,U._odata[sU],chi,Zp);
  accumReconZm(result,Uchi);
    
  // Tm
  offset = st._offsets [Tp][ss];
  local  = st._is_local[Tp][ss];
  perm   = st._permute[Tp][ss];
  ptype  = st._permute_type[Tp];
  if ( local && perm ) {
    spProjTm(tmp,in._odata[offset]);
    permute(chi,tmp,ptype);
  } else if ( local ) {
    spProjTm(chi,in._odata[offset]);
  } else { 
    chi=buf[offset];
  }
  Impl::multLink(Uchi,U._odata[sU],chi,Tp);
  accumReconTm(result,Uchi);
  
  vstream(out._odata[ss],result*(-0.5));
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
  int offset,local,perm, ptype;
  int ss=sF;

  offset = st._offsets [dir][ss];
  local  = st._is_local[dir][ss];
  perm   = st._permute[dir][ss];
  ptype  = st._permute_type[dir];

  // Xp
  if(gamma==Xp){
    if ( local && perm ) {
      spProjXp(tmp,in._odata[offset]);
      permute(chi,tmp,ptype);
    } else if ( local ) {
      spProjXp(chi,in._odata[offset]);
    } else { 
      chi=buf[offset];
    }
    Impl::multLink(Uchi,U._odata[sU],chi,dir);
    spReconXp(result,Uchi);
  }

  // Yp
  if ( gamma==Yp ){
    if ( local && perm ) {
      spProjYp(tmp,in._odata[offset]);
      permute(chi,tmp,ptype);
    } else if ( local ) {
      spProjYp(chi,in._odata[offset]);
    } else { 
      chi=buf[offset];
    }
    Impl::multLink(Uchi,U._odata[sU],chi,dir);
    spReconYp(result,Uchi);
  }
  
  // Zp
  if ( gamma ==Zp ){
    if ( local && perm ) {
      spProjZp(tmp,in._odata[offset]);
      permute(chi,tmp,ptype);
    } else if ( local ) {
      spProjZp(chi,in._odata[offset]);
    } else { 
      chi=buf[offset];
    }
    Impl::multLink(Uchi,U._odata[sU],chi,dir);
    spReconZp(result,Uchi);
  }
  
  // Tp
  if ( gamma ==Tp ){
    if ( local && perm ) {
      spProjTp(tmp,in._odata[offset]);
      permute(chi,tmp,ptype);
    } else if ( local ) {
      spProjTp(chi,in._odata[offset]);
    } else { 
      chi=buf[offset];
    }
    Impl::multLink(Uchi,U._odata[sU],chi,dir);
    spReconTp(result,Uchi);
  }

  // Xm
  if ( gamma==Xm ){
    if ( local && perm ) {
      spProjXm(tmp,in._odata[offset]);
      permute(chi,tmp,ptype);
    } else if ( local ) {
      spProjXm(chi,in._odata[offset]);
    } else { 
      chi=buf[offset];
    }
    Impl::multLink(Uchi,U._odata[sU],chi,dir);
    spReconXm(result,Uchi);
  }

  // Ym
  if ( gamma == Ym ){
    if ( local && perm ) {
      spProjYm(tmp,in._odata[offset]);
      permute(chi,tmp,ptype);
    } else if ( local ) {
      spProjYm(chi,in._odata[offset]);
    } else { 
      chi=buf[offset];
    }
    Impl::multLink(Uchi,U._odata[sU],chi,dir);
    spReconYm(result,Uchi);
  }

  // Zm
  if ( gamma == Zm ){
    if ( local && perm ) {
      spProjZm(tmp,in._odata[offset]);
      permute(chi,tmp,ptype);
    } else if ( local ) {
      spProjZm(chi,in._odata[offset]);
    } else { 
      chi=buf[offset];
    }
    Impl::multLink(Uchi,U._odata[sU],chi,dir);
    spReconZm(result,Uchi);
  }
  
  // Tm
  if ( gamma==Tm ) {
    if ( local && perm ) {
      spProjTm(tmp,in._odata[offset]);
      permute(chi,tmp,ptype);
    } else if ( local ) {
      spProjTm(chi,in._odata[offset]);
    } else { 
      chi=buf[offset];
    }
    Impl::multLink(Uchi,U._odata[sU],chi,dir);
    spReconTm(result,Uchi);
  }

  vstream(out._odata[ss],result*(-0.5));
}

  FermOpTemplateInstantiate(WilsonKernels);

}}
