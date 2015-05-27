#include <Grid.h>

namespace Grid {
namespace QCD {

const std::vector<int> WilsonMatrix::directions   ({0,1,2,3, 0, 1, 2, 3});
const std::vector<int> WilsonMatrix::displacements({1,1,1,1,-1,-1,-1,-1});

  int WilsonMatrix::HandOptDslash;

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
      case Xp:
	spProjXp(ret,in);
	break;
      case Yp:
	spProjYp(ret,in);
	break;
      case Zp:
	spProjZp(ret,in);
	break;
      case Tp:
	spProjTp(ret,in);
	break;
      case Xm:
	spProjXm(ret,in);
	break;
      case Ym:
	spProjYm(ret,in);
	break;
      case Zm:
	spProjZm(ret,in);
	break;
      case Tm:
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

void WilsonMatrix::DhopInternal(CartesianStencil & st,LatticeDoubledGaugeField & U,
				const LatticeFermion &in, LatticeFermion &out,int dag)
{
  assert((dag==DaggerNo) ||(dag==DaggerYes));
  WilsonCompressor compressor(dag);
  st.HaloExchange<vSpinColourVector,vHalfSpinColourVector,WilsonCompressor>(in,comm_buf,compressor);

  if ( dag == DaggerYes ) {
    if( HandOptDslash ) {
PARALLEL_FOR_LOOP
      for(int sss=0;sss<in._grid->oSites();sss++){
        DiracOptHand::DhopSiteDag(st,U,comm_buf,sss,in,out);
      }
    } else { 
PARALLEL_FOR_LOOP
      for(int sss=0;sss<in._grid->oSites();sss++){
        DiracOpt::DhopSiteDag(st,U,comm_buf,sss,in,out);
      }
    }
  } else {
    if( HandOptDslash ) {
PARALLEL_FOR_LOOP
      for(int sss=0;sss<in._grid->oSites();sss++){
        DiracOptHand::DhopSite(st,U,comm_buf,sss,in,out);
      }
    } else { 
PARALLEL_FOR_LOOP
      for(int sss=0;sss<in._grid->oSites();sss++){
        DiracOpt::DhopSite(st,U,comm_buf,sss,in,out);
      }
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



