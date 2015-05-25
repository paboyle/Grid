
#include <Grid.h>

namespace Grid {
namespace QCD {

const std::vector<int> WilsonMatrix::directions   ({0,1,2,3, 0, 1, 2, 3});
const std::vector<int> WilsonMatrix::displacements({1,1,1,1,-1,-1,-1,-1});

  // Should be in header?
const int WilsonMatrix::Xp = 0;
const int WilsonMatrix::Yp = 1;
const int WilsonMatrix::Zp = 2;
const int WilsonMatrix::Tp = 3;
const int WilsonMatrix::Xm = 4;
const int WilsonMatrix::Ym = 5;
const int WilsonMatrix::Zm = 6;
const int WilsonMatrix::Tm = 7;

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

  WilsonMatrix::WilsonMatrix(LatticeGaugeField &_Umu,GridCartesian &Fgrid,GridRedBlackCartesian &Hgrid, double _mass)  : 
    CheckerBoardedSparseMatrixBase<LatticeFermion>(&Fgrid,&Hgrid),
    Stencil    (  _grid,npoint,Even,directions,displacements),
    StencilEven(_cbgrid,npoint,Even,directions,displacements), // source is Even
    StencilOdd (_cbgrid,npoint,Odd ,directions,displacements), // source is Odd
    mass(_mass),
    Umu(_grid),
    UmuEven(_cbgrid),
    UmuOdd (_cbgrid)
  {
    // Allocate the required comms buffer
    comm_buf.resize(Stencil._unified_buffer_size); // this is always big enough to contain EO
    
    DoubleStore(Umu,_Umu);
    pickCheckerboard(Even,UmuEven,Umu);
    pickCheckerboard(Odd ,UmuOdd,Umu);
  }
      
void WilsonMatrix::DoubleStore(LatticeDoubledGaugeField &Uds,const LatticeGaugeField &Umu)
{
  LatticeColourMatrix U(_grid);

  for(int mu=0;mu<Nd;mu++){
    U = peekIndex<LorentzIndex>(Umu,mu);
    pokeIndex<LorentzIndex>(Uds,U,mu);
    U = adj(Cshift(U,mu,-1));
    pokeIndex<LorentzIndex>(Uds,U,mu+4);
  }
}

RealD WilsonMatrix::M(const LatticeFermion &in, LatticeFermion &out)
{
  out.checkerboard=in.checkerboard;
  Dhop(in,out,DaggerNo);
  out = (4+mass)*in - 0.5*out  ; // FIXME : axpby_norm! fusion fun
  return norm2(out);
}
RealD WilsonMatrix::Mdag(const LatticeFermion &in, LatticeFermion &out)
{
  out.checkerboard=in.checkerboard;
  Dhop(in,out,DaggerYes);
  out = (4+mass)*in - 0.5*out  ; // FIXME : axpby_norm! fusion fun
  return norm2(out);
}

void WilsonMatrix::Meooe(const LatticeFermion &in, LatticeFermion &out)
{
  if ( in.checkerboard == Odd ) {
    DhopEO(in,out,DaggerNo);
  } else {
    DhopOE(in,out,DaggerNo);
  }
  out = (-0.5)*out; // FIXME : scale factor in Dhop
}
void WilsonMatrix::MeooeDag(const LatticeFermion &in, LatticeFermion &out)
{
  if ( in.checkerboard == Odd ) {
    DhopEO(in,out,DaggerYes);
  } else {
    DhopOE(in,out,DaggerYes);
  }
  out = (-0.5)*out; // FIXME : scale factor in Dhop
}
void WilsonMatrix::Mooee(const LatticeFermion &in, LatticeFermion &out)
{
  out.checkerboard = in.checkerboard;
  out = (4.0+mass)*in;
  return ;
}
void WilsonMatrix::MooeeDag(const LatticeFermion &in, LatticeFermion &out)
{
  out.checkerboard = in.checkerboard;
  Mooee(in,out);
}
void WilsonMatrix::MooeeInv(const LatticeFermion &in, LatticeFermion &out)
{
  out.checkerboard = in.checkerboard;
  out = (1.0/(4.0+mass))*in;
  return ;
}
void WilsonMatrix::MooeeInvDag(const LatticeFermion &in, LatticeFermion &out)
{
  out.checkerboard = in.checkerboard;
  MooeeInv(in,out);
}

void WilsonMatrix::DhopSite(CartesianStencil &st,LatticeDoubledGaugeField &U,
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

void WilsonMatrix::DhopSiteDag(CartesianStencil &st,LatticeDoubledGaugeField &U,
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

void WilsonMatrix::DhopInternal(CartesianStencil & st,LatticeDoubledGaugeField & U,
				const LatticeFermion &in, LatticeFermion &out,int dag)
{
  assert((dag==DaggerNo) ||(dag==DaggerYes));
  WilsonCompressor compressor(dag);

  st.HaloExchange<vSpinColourVector,vHalfSpinColourVector,WilsonCompressor>(in,comm_buf,compressor);

  if ( dag == DaggerYes ) {
PARALLEL_FOR_LOOP
    for(int sss=0;sss<in._grid->oSites();sss++){
      DhopSiteDag(st,U,comm_buf,sss,in,out);
    }
  } else {
PARALLEL_FOR_LOOP
    for(int sss=0;sss<in._grid->oSites();sss++){
      DhopSite(st,U,comm_buf,sss,in,out);
    }
  }
}
void WilsonMatrix::DhopOE(const LatticeFermion &in, LatticeFermion &out,int dag)
{
  conformable(in._grid,_cbgrid);    // verifies half grid
  conformable(in._grid,out._grid); // drops the cb check

  assert(in.checkerboard==Even);
  out.checkerboard = Odd;

  DhopInternal(StencilEven,UmuOdd,in,out,dag);
}
void WilsonMatrix::DhopEO(const LatticeFermion &in, LatticeFermion &out,int dag)
{
  conformable(in._grid,_cbgrid);    // verifies half grid
  conformable(in._grid,out._grid); // drops the cb check

  assert(in.checkerboard==Odd);
  out.checkerboard = Even;

  DhopInternal(StencilOdd,UmuEven,in,out,dag);
}
void WilsonMatrix::Dhop(const LatticeFermion &in, LatticeFermion &out,int dag)
{
  conformable(in._grid,_grid); // verifies full grid
  conformable(in._grid,out._grid);

  out.checkerboard = in.checkerboard;

  DhopInternal(Stencil,Umu,in,out,dag);
}

}}



