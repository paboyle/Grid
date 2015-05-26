#include <Grid.h>

namespace Grid {
namespace QCD {

void DiracOpt::DhopSite(CartesianStencil &st,LatticeDoubledGaugeField &U,
			    std::vector<vHalfSpinColourVector,alignedAllocator<vHalfSpinColourVector> >  &buf,
			    int ss,const LatticeFermion &in, LatticeFermion &out)
{
    vHalfSpinColourVector  tmp;    
    vHalfSpinColourVector  chi;    
    vSpinColourVector result;
    vHalfSpinColourVector Uchi;
    int offset,local,perm, ptype;

    //#define VERBOSE( A)  if ( ss<10 ) { std::cout << "site " <<ss << " " #A " neigh " << offset << " perm "<< perm <<std::endl;}    

    // Xp
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
    mult(&Uchi(),&U._odata[ss](Xp),&chi());
    spReconXp(result,Uchi);

    //    std::cout << "XP_RECON"<<std::endl;
    //    std::cout << result()(0)(0) <<" "<<result()(0)(1) <<" "<<result()(0)(2) <<std::endl;
    //    std::cout << result()(1)(0) <<" "<<result()(1)(1) <<" "<<result()(1)(2) <<std::endl;
    //    std::cout << result()(2)(0) <<" "<<result()(2)(1) <<" "<<result()(2)(2) <<std::endl;
    //    std::cout << result()(3)(0) <<" "<<result()(3)(1) <<" "<<result()(3)(2) <<std::endl;

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
    mult(&Uchi(),&U._odata[ss](Yp),&chi());
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
    mult(&Uchi(),&U._odata[ss](Zp),&chi());
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
    mult(&Uchi(),&U._odata[ss](Tp),&chi());
    accumReconTp(result,Uchi);

    // Xm
    offset = st._offsets [Xm][ss];
    local  = st._is_local[Xm][ss];
    perm   = st._permute[Xm][ss];
    ptype  = st._permute_type[Xm];

    if ( local && perm ) 
    {
      spProjXm(tmp,in._odata[offset]);
      permute(chi,tmp,ptype);
    } else if ( local ) {
      spProjXm(chi,in._odata[offset]);
    } else { 
      chi=buf[offset];
    }
    mult(&Uchi(),&U._odata[ss](Xm),&chi());
    accumReconXm(result,Uchi);
    //  std::cout << "XM_RECON_ACCUM"<<std::endl;
    //    std::cout << result()(0)(0) <<" "<<result()(0)(1) <<" "<<result()(0)(2) <<std::endl;
    //    std::cout << result()(1)(0) <<" "<<result()(1)(1) <<" "<<result()(1)(2) <<std::endl;
    //    std::cout << result()(2)(0) <<" "<<result()(2)(1) <<" "<<result()(2)(2) <<std::endl;
    //    std::cout << result()(3)(0) <<" "<<result()(3)(1) <<" "<<result()(3)(2) <<std::endl;


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
    mult(&Uchi(),&U._odata[ss](Ym),&chi());
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
    mult(&Uchi(),&U._odata[ss](Zm),&chi());
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
    mult(&Uchi(),&U._odata[ss](Tm),&chi());
    accumReconTm(result,Uchi);

    vstream(out._odata[ss],result);
}

void DiracOpt::DhopSiteDag(CartesianStencil &st,LatticeDoubledGaugeField &U,
			       std::vector<vHalfSpinColourVector,alignedAllocator<vHalfSpinColourVector> >  &buf,
			       int ss,const LatticeFermion &in, LatticeFermion &out)
{
    vHalfSpinColourVector  tmp;    
    vHalfSpinColourVector  chi;    
    vSpinColourVector result;
    vHalfSpinColourVector Uchi;
    int offset,local,perm, ptype;

    // Xp
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
    mult(&Uchi(),&U._odata[ss](Xm),&chi());
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
    mult(&Uchi(),&U._odata[ss](Ym),&chi());
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
    mult(&Uchi(),&U._odata[ss](Zm),&chi());
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
    mult(&Uchi(),&U._odata[ss](Tm),&chi());
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
    mult(&Uchi(),&U._odata[ss](Xp),&chi());
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
    mult(&Uchi(),&U._odata[ss](Yp),&chi());
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
    mult(&Uchi(),&U._odata[ss](Zp),&chi());
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
    mult(&Uchi(),&U._odata[ss](Tp),&chi());
    accumReconTm(result,Uchi);

    vstream(out._odata[ss],result);
}
}}
