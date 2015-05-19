
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
  Dhop(in,out,0);
  out = (4+mass)*in - 0.5*out  ; // FIXME : axpby_norm! fusion fun
  return norm2(out);
}
RealD WilsonMatrix::Mdag(const LatticeFermion &in, LatticeFermion &out)
{
  Dhop(in,out,1);
  out = (4+mass)*in - 0.5*out  ; // FIXME : axpby_norm! fusion fun
  return norm2(out);
}

void WilsonMatrix::Meooe(const LatticeFermion &in, LatticeFermion &out)
{
  Dhop(in,out,0);
  out = 0.5*out; // FIXME : scale factor in Dhop
}
void WilsonMatrix::MeooeDag(const LatticeFermion &in, LatticeFermion &out)
{
  Dhop(in,out,1);
}
void WilsonMatrix::Mooee(const LatticeFermion &in, LatticeFermion &out)
{
  out = (4.0+mass)*in;
  return ;
}
void WilsonMatrix::MooeeInv(const LatticeFermion &in, LatticeFermion &out)
{
  out = (1.0/(4.0+mass))*in;
  return ;
}
void WilsonMatrix::MooeeDag(const LatticeFermion &in, LatticeFermion &out)
{
  out = (1.0/(4.0+mass))*in;
  return ;
}
void WilsonMatrix::MooeeInvDag(const LatticeFermion &in, LatticeFermion &out)
{
  out = (1.0/(4.0+mass))*in;
  return ;
}

void WilsonMatrix::Dhop(const LatticeFermion &in, LatticeFermion &out,int dag)
{
  assert((dag==0) ||(dag==1));

  WilsonCompressor compressor(dag);
  Stencil.HaloExchange<vSpinColourVector,vHalfSpinColourVector,WilsonCompressor>(in,comm_buf,compressor);

PARALLEL_FOR_LOOP
  for(int sss=0;sss<grid->oSites();sss++){

    vHalfSpinColourVector  tmp;    
    vHalfSpinColourVector  chi;    
    vSpinColourVector result;
    vHalfSpinColourVector Uchi;
    int offset,local,perm, ptype;

    //    int ss = Stencil._LebesgueReorder[sss];
    int ss = sss;
    int ssu= ss;

    // Xp
    int idx = (Xp+dag*4)%8;
    offset = Stencil._offsets [idx][ss];
    local  = Stencil._is_local[idx][ss];
    perm   = Stencil._permute[idx][ss];

    ptype  = Stencil._permute_type[idx];
    if ( local && perm ) {
      spProjXp(tmp,in._odata[offset]);
      permute(chi,tmp,ptype);
    } else if ( local ) {
      spProjXp(chi,in._odata[offset]);
    } else { 
      chi=comm_buf[offset];
    }
    mult(&Uchi(),&Umu._odata[ssu](idx),&chi());
    spReconXp(result,Uchi);

    // Yp
    idx = (Yp+dag*4)%8;
    offset = Stencil._offsets [idx][ss];
    local  = Stencil._is_local[idx][ss];
    perm   = Stencil._permute[idx][ss];
    ptype  = Stencil._permute_type[idx];

    if ( local && perm ) {
      spProjYp(tmp,in._odata[offset]);
      permute(chi,tmp,ptype);
    } else if ( local ) {
      spProjYp(chi,in._odata[offset]);
    } else { 
      chi=comm_buf[offset];
    }
    mult(&Uchi(),&Umu._odata[ssu](idx),&chi());
    accumReconYp(result,Uchi);

    // Zp
    idx = (Zp+dag*4)%8;
    offset = Stencil._offsets [idx][ss];
    local  = Stencil._is_local[idx][ss];
    perm   = Stencil._permute[idx][ss];
    ptype  = Stencil._permute_type[idx];
    if ( local && perm ) {
      spProjZp(tmp,in._odata[offset]);
      permute(chi,tmp,ptype);
    } else if ( local ) {
      spProjZp(chi,in._odata[offset]);
    } else { 
      chi=comm_buf[offset];
    }
    mult(&Uchi(),&Umu._odata[ssu](idx),&chi());
    accumReconZp(result,Uchi);

    // Tp
    idx = (Tp+dag*4)%8;
    offset = Stencil._offsets [idx][ss];
    local  = Stencil._is_local[idx][ss];
    perm   = Stencil._permute[idx][ss];
    ptype  = Stencil._permute_type[idx];
    if ( local && perm ) {
      spProjTp(tmp,in._odata[offset]);
      permute(chi,tmp,ptype);
    } else if ( local ) {
      spProjTp(chi,in._odata[offset]);
    } else { 
      chi=comm_buf[offset];
    }
    mult(&Uchi(),&Umu._odata[ssu](idx),&chi());
    accumReconTp(result,Uchi);

    // Xm
    idx = (Xm+dag*4)%8;
    offset = Stencil._offsets [idx][ss];
    local  = Stencil._is_local[idx][ss];
    perm   = Stencil._permute[idx][ss];
    ptype  = Stencil._permute_type[idx];

    if ( local && perm ) 
    {
      spProjXm(tmp,in._odata[offset]);
      permute(chi,tmp,ptype);
    } else if ( local ) {
      spProjXm(chi,in._odata[offset]);
    } else { 
      chi=comm_buf[offset];
    }
    mult(&Uchi(),&Umu._odata[ssu](idx),&chi());
    accumReconXm(result,Uchi);


    // Ym
    idx = (Ym+dag*4)%8;
    offset = Stencil._offsets [idx][ss];
    local  = Stencil._is_local[idx][ss];
    perm   = Stencil._permute[idx][ss];
    ptype  = Stencil._permute_type[idx];

    if ( local && perm ) {
      spProjYm(tmp,in._odata[offset]);
      permute(chi,tmp,ptype);
    } else if ( local ) {
      spProjYm(chi,in._odata[offset]);
    } else { 
      chi=comm_buf[offset];
    }
    mult(&Uchi(),&Umu._odata[ssu](idx),&chi());
    accumReconYm(result,Uchi);

    // Zm
    idx = (Zm+dag*4)%8;
    offset = Stencil._offsets [idx][ss];
    local  = Stencil._is_local[idx][ss];
    perm   = Stencil._permute[idx][ss];
    ptype  = Stencil._permute_type[idx];
    if ( local && perm ) {
      spProjZm(tmp,in._odata[offset]);
      permute(chi,tmp,ptype);
    } else if ( local ) {
      spProjZm(chi,in._odata[offset]);
    } else { 
      chi=comm_buf[offset];
    }
    mult(&Uchi(),&Umu._odata[ssu](idx),&chi());
    accumReconZm(result,Uchi);

    // Tm
    idx = (Tm+dag*4)%8;
    offset = Stencil._offsets [idx][ss];
    local  = Stencil._is_local[idx][ss];
    perm   = Stencil._permute[idx][ss];
    ptype  = Stencil._permute_type[idx];
    if ( local && perm ) {
      spProjTm(tmp,in._odata[offset]);
      permute(chi,tmp,ptype);
    } else if ( local ) {
      spProjTm(chi,in._odata[offset]);
    } else { 
      chi=comm_buf[offset];
    }
    mult(&Uchi(),&Umu._odata[ssu](idx),&chi());
    accumReconTm(result,Uchi);

    vstream(out._odata[ss],result);
  }

}





}}

