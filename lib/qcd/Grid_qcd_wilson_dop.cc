
#include <Grid.h>

namespace Grid {
namespace QCD {

const std::vector<int> WilsonMatrix::directions   ({0,1,2,3, 0, 1, 2, 3,0});
const std::vector<int> WilsonMatrix::displacements({1,1,1,1,-1,-1,-1,-1,0});

  // Should be in header?
const int WilsonMatrix::Xp = 0;
const int WilsonMatrix::Yp = 1;
const int WilsonMatrix::Zp = 2;
const int WilsonMatrix::Tp = 3;
const int WilsonMatrix::Xm = 4;
const int WilsonMatrix::Ym = 5;
const int WilsonMatrix::Zm = 6;
const int WilsonMatrix::Tm = 7;
  //const int WilsonMatrix::X0 = 8;

  class WilsonCompressor {
  public:
    int mu;
    int dag;

    WilsonCompressor(int _dag){
      mu=0;
      dag=_dag;
      assert((dag==0)||(dag==1));
    }
    void Point(int p) { 
      mu=p;
    };

    vHalfSpinColourVector operator () (const vSpinColourVector &in)
    {
      vHalfSpinColourVector ret;
      int mudag=mu;
      if (dag) {
	mudag=(mu+Nd)%(2*Nd);
      }
      switch(mudag) {
      case WilsonMatrix::Xp:
	spProjXp(ret,in);
	break;
      case WilsonMatrix::Yp:
	spProjYp(ret,in);
	break;
      case WilsonMatrix::Zp:
	spProjZp(ret,in);
	break;
      case WilsonMatrix::Tp:
	spProjTp(ret,in);
	break;
      case WilsonMatrix::Xm:
	spProjXm(ret,in);
	break;
      case WilsonMatrix::Ym:
	spProjYm(ret,in);
	break;
      case WilsonMatrix::Zm:
	spProjZm(ret,in);
	break;
      case WilsonMatrix::Tm:
	spProjTm(ret,in);
	break;
      default: 
	assert(0);
	break;
      }
      return ret;
    }
  };

  WilsonMatrix::WilsonMatrix(LatticeGaugeField &_Umu,double _mass)
  : Stencil(Umu._grid,npoint,0,directions,displacements),
    mass(_mass),
    Umu(_Umu._grid)
{
  // Allocate the required comms buffer
  grid = _Umu._grid;
  comm_buf.resize(Stencil._unified_buffer_size);
  DoubleStore(Umu,_Umu);
}
      
void WilsonMatrix::DoubleStore(LatticeDoubledGaugeField &Uds,const LatticeGaugeField &Umu)
{
  LatticeColourMatrix U(grid);

  for(int mu=0;mu<Nd;mu++){
    U = peekIndex<LorentzIndex>(Umu,mu);
    pokeIndex<LorentzIndex>(Uds,U,mu);
    U = adj(Cshift(U,mu,-1));
    pokeIndex<LorentzIndex>(Uds,U,mu+4);
  }
}

RealD WilsonMatrix::M(const LatticeFermion &in, LatticeFermion &out)
{
  this->Dhop(in,out,0);
  out = (4+mass)*in - 0.5*out  ; // FIXME : axpby_norm! fusion fun
  return norm2(out);
}
RealD WilsonMatrix::Mdag(const LatticeFermion &in, LatticeFermion &out)
{
  this->Dhop(in,out,1);
  out = (4+mass)*in - 0.5*out  ; // FIXME : axpby_norm! fusion fun
  return norm2(out);
}

void WilsonMatrix::Meooe(const LatticeFermion &in, LatticeFermion &out)
{
  this->Dhop(in,out,0);
  out = 0.5*out; // FIXME : scale factor in Dhop
}
void WilsonMatrix::MeooeDag(const LatticeFermion &in, LatticeFermion &out)
{
  this->Dhop(in,out,1);
  out = 0.5*out; // FIXME : scale factor in Dhop
}
void WilsonMatrix::Mooee(const LatticeFermion &in, LatticeFermion &out)
{
  out = (4.0+mass)*in;
  return ;
}
void WilsonMatrix::MooeeDag(const LatticeFermion &in, LatticeFermion &out)
{
  this->Mooee(in,out);
}
void WilsonMatrix::MooeeInv(const LatticeFermion &in, LatticeFermion &out)
{
  out = (1.0/(4.0+mass))*in;
  return ;
}
void WilsonMatrix::MooeeInvDag(const LatticeFermion &in, LatticeFermion &out)
{
  this->MooeeInv(in,out);
}

void WilsonMatrix::DhopSite(int ss,const LatticeFermion &in, LatticeFermion &out)
{
    vHalfSpinColourVector  tmp;    
    vHalfSpinColourVector  chi;    
    vSpinColourVector result;
    vHalfSpinColourVector Uchi;
    int offset,local,perm, ptype;

    //    int ss = Stencil._LebesgueReorder[sss];
    int ssu= ss;

    // Xp
    offset = Stencil._offsets [Xp][ss];
    local  = Stencil._is_local[Xp][ss];
    perm   = Stencil._permute[Xp][ss];

    ptype  = Stencil._permute_type[Xp];
    if ( local && perm ) {
      spProjXp(tmp,in._odata[offset]);
      permute(chi,tmp,ptype);
    } else if ( local ) {
      spProjXp(chi,in._odata[offset]);
    } else { 
      chi=comm_buf[offset];
    }
    mult(&Uchi(),&Umu._odata[ssu](Xp),&chi());
    spReconXp(result,Uchi);

    // Yp
    offset = Stencil._offsets [Yp][ss];
    local  = Stencil._is_local[Yp][ss];
    perm   = Stencil._permute[Yp][ss];
    ptype  = Stencil._permute_type[Yp];
    if ( local && perm ) {
      spProjYp(tmp,in._odata[offset]);
      permute(chi,tmp,ptype);
    } else if ( local ) {
      spProjYp(chi,in._odata[offset]);
    } else { 
      chi=comm_buf[offset];
    }
    mult(&Uchi(),&Umu._odata[ssu](Yp),&chi());
    accumReconYp(result,Uchi);

    // Zp
    offset = Stencil._offsets [Zp][ss];
    local  = Stencil._is_local[Zp][ss];
    perm   = Stencil._permute[Zp][ss];
    ptype  = Stencil._permute_type[Zp];
    if ( local && perm ) {
      spProjZp(tmp,in._odata[offset]);
      permute(chi,tmp,ptype);
    } else if ( local ) {
      spProjZp(chi,in._odata[offset]);
    } else { 
      chi=comm_buf[offset];
    }
    mult(&Uchi(),&Umu._odata[ssu](Zp),&chi());
    accumReconZp(result,Uchi);

    // Tp
    offset = Stencil._offsets [Tp][ss];
    local  = Stencil._is_local[Tp][ss];
    perm   = Stencil._permute[Tp][ss];
    ptype  = Stencil._permute_type[Tp];
    if ( local && perm ) {
      spProjTp(tmp,in._odata[offset]);
      permute(chi,tmp,ptype);
    } else if ( local ) {
      spProjTp(chi,in._odata[offset]);
    } else { 
      chi=comm_buf[offset];
    }
    mult(&Uchi(),&Umu._odata[ssu](Tp),&chi());
    accumReconTp(result,Uchi);

    // Xm
    offset = Stencil._offsets [Xm][ss];
    local  = Stencil._is_local[Xm][ss];
    perm   = Stencil._permute[Xm][ss];
    ptype  = Stencil._permute_type[Xm];

    if ( local && perm ) 
    {
      spProjXm(tmp,in._odata[offset]);
      permute(chi,tmp,ptype);
    } else if ( local ) {
      spProjXm(chi,in._odata[offset]);
    } else { 
      chi=comm_buf[offset];
    }
    mult(&Uchi(),&Umu._odata[ssu](Xm),&chi());
    accumReconXm(result,Uchi);


    // Ym
    offset = Stencil._offsets [Ym][ss];
    local  = Stencil._is_local[Ym][ss];
    perm   = Stencil._permute[Ym][ss];
    ptype  = Stencil._permute_type[Ym];

    if ( local && perm ) {
      spProjYm(tmp,in._odata[offset]);
      permute(chi,tmp,ptype);
    } else if ( local ) {
      spProjYm(chi,in._odata[offset]);
    } else { 
      chi=comm_buf[offset];
    }
    mult(&Uchi(),&Umu._odata[ssu](Ym),&chi());
    accumReconYm(result,Uchi);

    // Zm
    offset = Stencil._offsets [Zm][ss];
    local  = Stencil._is_local[Zm][ss];
    perm   = Stencil._permute[Zm][ss];
    ptype  = Stencil._permute_type[Zm];
    if ( local && perm ) {
      spProjZm(tmp,in._odata[offset]);
      permute(chi,tmp,ptype);
    } else if ( local ) {
      spProjZm(chi,in._odata[offset]);
    } else { 
      chi=comm_buf[offset];
    }
    mult(&Uchi(),&Umu._odata[ssu](Zm),&chi());
    accumReconZm(result,Uchi);

    // Tm
    offset = Stencil._offsets [Tm][ss];
    local  = Stencil._is_local[Tm][ss];
    perm   = Stencil._permute[Tm][ss];
    ptype  = Stencil._permute_type[Tm];
    if ( local && perm ) {
      spProjTm(tmp,in._odata[offset]);
      permute(chi,tmp,ptype);
    } else if ( local ) {
      spProjTm(chi,in._odata[offset]);
    } else { 
      chi=comm_buf[offset];
    }
    mult(&Uchi(),&Umu._odata[ssu](Tm),&chi());
    accumReconTm(result,Uchi);

    vstream(out._odata[ss],result);
}
void WilsonMatrix::DhopSiteDag(int ss,const LatticeFermion &in, LatticeFermion &out)
{
    vHalfSpinColourVector  tmp;    
    vHalfSpinColourVector  chi;    
    vSpinColourVector result;
    vHalfSpinColourVector Uchi;
    int offset,local,perm, ptype;

    int ssu= ss;

    // Xp
    offset = Stencil._offsets [Xm][ss];
    local  = Stencil._is_local[Xm][ss];
    perm   = Stencil._permute[Xm][ss];

    ptype  = Stencil._permute_type[Xm];
    if ( local && perm ) {
      spProjXp(tmp,in._odata[offset]);
      permute(chi,tmp,ptype);
    } else if ( local ) {
      spProjXp(chi,in._odata[offset]);
    } else { 
      chi=comm_buf[offset];
    }
    mult(&Uchi(),&Umu._odata[ssu](Xm),&chi());
    spReconXp(result,Uchi);

    // Yp
    offset = Stencil._offsets [Ym][ss];
    local  = Stencil._is_local[Ym][ss];
    perm   = Stencil._permute[Ym][ss];
    ptype  = Stencil._permute_type[Ym];
    if ( local && perm ) {
      spProjYp(tmp,in._odata[offset]);
      permute(chi,tmp,ptype);
    } else if ( local ) {
      spProjYp(chi,in._odata[offset]);
    } else { 
      chi=comm_buf[offset];
    }
    mult(&Uchi(),&Umu._odata[ssu](Ym),&chi());
    accumReconYp(result,Uchi);

    // Zp
    offset = Stencil._offsets [Zm][ss];
    local  = Stencil._is_local[Zm][ss];
    perm   = Stencil._permute[Zm][ss];
    ptype  = Stencil._permute_type[Zm];
    if ( local && perm ) {
      spProjZp(tmp,in._odata[offset]);
      permute(chi,tmp,ptype);
    } else if ( local ) {
      spProjZp(chi,in._odata[offset]);
    } else { 
      chi=comm_buf[offset];
    }
    mult(&Uchi(),&Umu._odata[ssu](Zm),&chi());
    accumReconZp(result,Uchi);

    // Tp
    offset = Stencil._offsets [Tm][ss];
    local  = Stencil._is_local[Tm][ss];
    perm   = Stencil._permute[Tm][ss];
    ptype  = Stencil._permute_type[Tm];
    if ( local && perm ) {
      spProjTp(tmp,in._odata[offset]);
      permute(chi,tmp,ptype);
    } else if ( local ) {
      spProjTp(chi,in._odata[offset]);
    } else { 
      chi=comm_buf[offset];
    }
    mult(&Uchi(),&Umu._odata[ssu](Tm),&chi());
    accumReconTp(result,Uchi);

    // Xm
    offset = Stencil._offsets [Xp][ss];
    local  = Stencil._is_local[Xp][ss];
    perm   = Stencil._permute[Xp][ss];
    ptype  = Stencil._permute_type[Xp];

    if ( local && perm ) 
    {
      spProjXm(tmp,in._odata[offset]);
      permute(chi,tmp,ptype);
    } else if ( local ) {
      spProjXm(chi,in._odata[offset]);
    } else { 
      chi=comm_buf[offset];
    }
    mult(&Uchi(),&Umu._odata[ssu](Xp),&chi());
    accumReconXm(result,Uchi);


    // Ym
    offset = Stencil._offsets [Yp][ss];
    local  = Stencil._is_local[Yp][ss];
    perm   = Stencil._permute[Yp][ss];
    ptype  = Stencil._permute_type[Yp];

    if ( local && perm ) {
      spProjYm(tmp,in._odata[offset]);
      permute(chi,tmp,ptype);
    } else if ( local ) {
      spProjYm(chi,in._odata[offset]);
    } else { 
      chi=comm_buf[offset];
    }
    mult(&Uchi(),&Umu._odata[ssu](Yp),&chi());
    accumReconYm(result,Uchi);

    // Zm
    offset = Stencil._offsets [Zp][ss];
    local  = Stencil._is_local[Zp][ss];
    perm   = Stencil._permute[Zp][ss];
    ptype  = Stencil._permute_type[Zp];
    if ( local && perm ) {
      spProjZm(tmp,in._odata[offset]);
      permute(chi,tmp,ptype);
    } else if ( local ) {
      spProjZm(chi,in._odata[offset]);
    } else { 
      chi=comm_buf[offset];
    }
    mult(&Uchi(),&Umu._odata[ssu](Zp),&chi());
    accumReconZm(result,Uchi);

    // Tm
    offset = Stencil._offsets [Tp][ss];
    local  = Stencil._is_local[Tp][ss];
    perm   = Stencil._permute[Tp][ss];
    ptype  = Stencil._permute_type[Tp];
    if ( local && perm ) {
      spProjTm(tmp,in._odata[offset]);
      permute(chi,tmp,ptype);
    } else if ( local ) {
      spProjTm(chi,in._odata[offset]);
    } else { 
      chi=comm_buf[offset];
    }
    mult(&Uchi(),&Umu._odata[ssu](Tp),&chi());
    accumReconTm(result,Uchi);

    vstream(out._odata[ss],result);
}

void WilsonMatrix::Dhop(const LatticeFermion &in, LatticeFermion &out,int dag)
{
  assert((dag==0) ||(dag==1));

  WilsonCompressor compressor(dag);
  Stencil.HaloExchange<vSpinColourVector,vHalfSpinColourVector,WilsonCompressor>(in,comm_buf,compressor);

  if ( dag ) {
PARALLEL_FOR_LOOP
    for(int sss=0;sss<grid->oSites();sss++){
      DhopSiteDag(sss,in,out);
    }
  } else {
PARALLEL_FOR_LOOP
    for(int sss=0;sss<grid->oSites();sss++){
      DhopSite(sss,in,out);
    }
  }
}


}
}

