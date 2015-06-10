#include <Grid.h>

namespace Grid {
namespace QCD {

void DiracOpt::DhopSite(CartesianStencil &st,LatticeDoubledGaugeField &U,
			std::vector<vHalfSpinColourVector,alignedAllocator<vHalfSpinColourVector> >  &buf,
			int sF,int sU,const LatticeFermion &in, LatticeFermion &out)
{
    vHalfSpinColourVector  tmp;    
    vHalfSpinColourVector  chi;    
    vSpinColourVector result;
    vHalfSpinColourVector Uchi;
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
    mult(&Uchi(),&U._odata[sU](Xp),&chi());
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
    mult(&Uchi(),&U._odata[sU](Yp),&chi());
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
    mult(&Uchi(),&U._odata[sU](Zp),&chi());
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
    mult(&Uchi(),&U._odata[sU](Tp),&chi());
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
    mult(&Uchi(),&U._odata[sU](Xm),&chi());
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
    mult(&Uchi(),&U._odata[sU](Ym),&chi());
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
    mult(&Uchi(),&U._odata[sU](Zm),&chi());
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
    mult(&Uchi(),&U._odata[sU](Tm),&chi());
    accumReconTm(result,Uchi);

    vstream(out._odata[ss],result*(-0.5));
}

void DiracOpt::DhopSiteDag(CartesianStencil &st,LatticeDoubledGaugeField &U,
			   std::vector<vHalfSpinColourVector,alignedAllocator<vHalfSpinColourVector> >  &buf,
			   int sF,int sU,const LatticeFermion &in, LatticeFermion &out)
{
    vHalfSpinColourVector  tmp;    
    vHalfSpinColourVector  chi;    
    vSpinColourVector result;
    vHalfSpinColourVector Uchi;
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
    mult(&Uchi(),&U._odata[sU](Xm),&chi());
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
    mult(&Uchi(),&U._odata[sU](Ym),&chi());
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
    mult(&Uchi(),&U._odata[sU](Zm),&chi());
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
    mult(&Uchi(),&U._odata[sU](Tm),&chi());
    accumReconTp(result,Uchi);

    // Xm
    offset = st._offsets [Xp][ss];
    local  = st._is_local[Xp][ss];
    perm   = st._permute[Xp][ss];
    ptype  = st._permute_type[Xp];

    if ( local && perm ) 
    {
      spProjXm(tmp,in._odata[offset]);
      permute(chi,tmp,ptype);
    } else if ( local ) {
      spProjXm(chi,in._odata[offset]);
    } else { 
      chi=buf[offset];
    }
    mult(&Uchi(),&U._odata[sU](Xp),&chi());
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
    mult(&Uchi(),&U._odata[sU](Yp),&chi());
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
    mult(&Uchi(),&U._odata[sU](Zp),&chi());
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
    mult(&Uchi(),&U._odata[sU](Tp),&chi());
    accumReconTm(result,Uchi);

    vstream(out._odata[ss],result*(-0.5));
}

void DiracOpt::DhopDir(CartesianStencil &st,LatticeDoubledGaugeField &U,
			std::vector<vHalfSpinColourVector,alignedAllocator<vHalfSpinColourVector> >  &buf,
		       int sF,int sU,const LatticeFermion &in, LatticeFermion &out,int dirdisp)
{
    vHalfSpinColourVector  tmp;    
    vHalfSpinColourVector  chi;    
    vSpinColourVector result;
    vHalfSpinColourVector Uchi;
    int offset,local,perm, ptype;
    int ss=sF;
    
    offset = st._offsets [dirdisp][ss];
    local  = st._is_local[dirdisp][ss];
    perm   = st._permute[dirdisp][ss];
    ptype  = st._permute_type[dirdisp];

    // Xp
    if(dirdisp==Xp){
      if ( local && perm ) {
	spProjXp(tmp,in._odata[offset]);
	permute(chi,tmp,ptype);
      } else if ( local ) {
	spProjXp(chi,in._odata[offset]);
      } else { 
	chi=buf[offset];
      }
      mult(&Uchi(),&U._odata[sU](Xp),&chi());
      spReconXp(result,Uchi);
    }

    // Yp
    if ( dirdisp==Yp ){
      if ( local && perm ) {
	spProjYp(tmp,in._odata[offset]);
	permute(chi,tmp,ptype);
      } else if ( local ) {
	spProjYp(chi,in._odata[offset]);
      } else { 
	chi=buf[offset];
      }
      mult(&Uchi(),&U._odata[sU](Yp),&chi());
      spReconYp(result,Uchi);
    }

    // Zp
    if ( dirdisp ==Zp ){
      if ( local && perm ) {
	spProjZp(tmp,in._odata[offset]);
	permute(chi,tmp,ptype);
      } else if ( local ) {
	spProjZp(chi,in._odata[offset]);
      } else { 
	chi=buf[offset];
      }
      mult(&Uchi(),&U._odata[sU](Zp),&chi());
      spReconZp(result,Uchi);
    }

    // Tp
    if ( dirdisp ==Tp ){
      if ( local && perm ) {
	spProjTp(tmp,in._odata[offset]);
	permute(chi,tmp,ptype);
      } else if ( local ) {
	spProjTp(chi,in._odata[offset]);
      } else { 
	chi=buf[offset];
      }
      mult(&Uchi(),&U._odata[sU](Tp),&chi());
      spReconTp(result,Uchi);
    }

    // Xm
    if ( dirdisp==Xm ){
      if ( local && perm ) {
	spProjXm(tmp,in._odata[offset]);
	permute(chi,tmp,ptype);
      } else if ( local ) {
	spProjXm(chi,in._odata[offset]);
      } else { 
	chi=buf[offset];
      }
      mult(&Uchi(),&U._odata[sU](Xm),&chi());
      spReconXm(result,Uchi);
    }

    // Ym
    if ( dirdisp == Ym ){
      if ( local && perm ) {
	spProjYm(tmp,in._odata[offset]);
	permute(chi,tmp,ptype);
      } else if ( local ) {
	spProjYm(chi,in._odata[offset]);
      } else { 
	chi=buf[offset];
      }
      mult(&Uchi(),&U._odata[sU](Ym),&chi());
      spReconYm(result,Uchi);
    }

    // Zm
    if ( dirdisp == Zm ){
      if ( local && perm ) {
	spProjZm(tmp,in._odata[offset]);
	permute(chi,tmp,ptype);
      } else if ( local ) {
	spProjZm(chi,in._odata[offset]);
      } else { 
	chi=buf[offset];
      }
      mult(&Uchi(),&U._odata[sU](Zm),&chi());
      spReconZm(result,Uchi);
    }

    // Tm
    if ( dirdisp==Tm ) {
      if ( local && perm ) {
	spProjTm(tmp,in._odata[offset]);
	permute(chi,tmp,ptype);
      } else if ( local ) {
	spProjTm(chi,in._odata[offset]);
      } else { 
	chi=buf[offset];
      }
      mult(&Uchi(),&U._odata[sU](Tm),&chi());
      spReconTm(result,Uchi);
    }

    vstream(out._odata[ss],result*(-0.5));
}

}}
