#include <Grid.h>

namespace Grid {
namespace QCD {

const std::vector<int> WilsonFermion::directions   ({0,1,2,3, 0, 1, 2, 3});
const std::vector<int> WilsonFermion::displacements({1,1,1,1,-1,-1,-1,-1});

int WilsonFermion::HandOptDslash;

WilsonFermion::WilsonFermion(LatticeGaugeField &_Umu,
			     GridCartesian         &Fgrid,
			     GridRedBlackCartesian &Hgrid, 
			     RealD _mass) :
  _grid(&Fgrid),
  _cbgrid(&Hgrid),
  Stencil    (&Fgrid,npoint,Even,directions,displacements),
  StencilEven(&Hgrid,npoint,Even,directions,displacements), // source is Even
  StencilOdd (&Hgrid,npoint,Odd ,directions,displacements), // source is Odd
  mass(_mass),
  Umu(&Fgrid),
  UmuEven(&Hgrid),
  UmuOdd (&Hgrid)
{
  // Allocate the required comms buffer
  comm_buf.resize(Stencil._unified_buffer_size); // this is always big enough to contain EO
  DoubleStore(Umu,_Umu);
  pickCheckerboard(Even,UmuEven,Umu);
  pickCheckerboard(Odd ,UmuOdd,Umu);
}
      
void WilsonFermion::DoubleStore(LatticeDoubledGaugeField &Uds,const LatticeGaugeField &Umu)
{
  conformable(Uds._grid,GaugeGrid());
  conformable(Umu._grid,GaugeGrid());
  LatticeColourMatrix U(GaugeGrid());
  for(int mu=0;mu<Nd;mu++){
    U = PeekIndex<LorentzIndex>(Umu,mu);
    PokeIndex<LorentzIndex>(Uds,U,mu);
    U = adj(Cshift(U,mu,-1));
    PokeIndex<LorentzIndex>(Uds,U,mu+4);
  }
}

RealD WilsonFermion::M(const LatticeFermion &in, LatticeFermion &out)
{
  out.checkerboard=in.checkerboard;
  Dhop(in,out,DaggerNo);
  return axpy_norm(out,4+mass,in,out);
}
RealD WilsonFermion::Mdag(const LatticeFermion &in, LatticeFermion &out)
{
  out.checkerboard=in.checkerboard;
  Dhop(in,out,DaggerYes);
  return axpy_norm(out,4+mass,in,out);
}

void WilsonFermion::Meooe(const LatticeFermion &in, LatticeFermion &out)
{
  if ( in.checkerboard == Odd ) {
    DhopEO(in,out,DaggerNo);
  } else {
    DhopOE(in,out,DaggerNo);
  }
}
void WilsonFermion::MeooeDag(const LatticeFermion &in, LatticeFermion &out)
{
  if ( in.checkerboard == Odd ) {
    DhopEO(in,out,DaggerYes);
  } else {
    DhopOE(in,out,DaggerYes);
  }
}
void WilsonFermion::Mooee(const LatticeFermion &in, LatticeFermion &out)
{
  out.checkerboard = in.checkerboard;
  out = (4.0+mass)*in;
  return ;
}
void WilsonFermion::MooeeDag(const LatticeFermion &in, LatticeFermion &out)
{
  out.checkerboard = in.checkerboard;
  Mooee(in,out);
}
void WilsonFermion::MooeeInv(const LatticeFermion &in, LatticeFermion &out)
{
  out.checkerboard = in.checkerboard;
  out = (1.0/(4.0+mass))*in;
  return ;
}
void WilsonFermion::MooeeInvDag(const LatticeFermion &in, LatticeFermion &out)
{
  out.checkerboard = in.checkerboard;
  MooeeInv(in,out);
}
void WilsonFermion::Mdir (const LatticeFermion &in, LatticeFermion &out,int dir,int disp)
{
  DhopDir(in,out,dir,disp);
}
void WilsonFermion::DhopDir(const LatticeFermion &in, LatticeFermion &out,int dir,int disp){

  WilsonCompressor compressor(DaggerNo);

  Stencil.HaloExchange<vSpinColourVector,vHalfSpinColourVector,WilsonCompressor>(in,comm_buf,compressor);
  
  assert( (disp==1)||(disp==-1) );

  int skip = (disp==1) ? 0 : 1;

  int dirdisp = dir+skip*4;

PARALLEL_FOR_LOOP
  for(int sss=0;sss<in._grid->oSites();sss++){
    DiracOptDhopDir(Stencil,Umu,comm_buf,sss,sss,in,out,dirdisp,dirdisp);
  }
 
};
void WilsonFermion::DhopDirDisp(const LatticeFermion &in, LatticeFermion &out,int dirdisp,int gamma,int dag)
{

  WilsonCompressor compressor(dag);

  Stencil.HaloExchange<vSpinColourVector,vHalfSpinColourVector,WilsonCompressor>(in,comm_buf,compressor);

PARALLEL_FOR_LOOP
  for(int sss=0;sss<in._grid->oSites();sss++){
    DiracOptDhopDir(Stencil,Umu,comm_buf,sss,sss,in,out,dirdisp,gamma);
  }
 
};

void WilsonFermion::DhopInternal(CartesianStencil & st,LatticeDoubledGaugeField & U,
				const LatticeFermion &in, LatticeFermion &out,int dag)
{
  assert((dag==DaggerNo) ||(dag==DaggerYes));
  WilsonCompressor compressor(dag);
  st.HaloExchange<vSpinColourVector,vHalfSpinColourVector,WilsonCompressor>(in,comm_buf,compressor);

  if ( dag == DaggerYes ) {
    if( HandOptDslash ) {
PARALLEL_FOR_LOOP
      for(int sss=0;sss<in._grid->oSites();sss++){
        DiracOptHandDhopSiteDag(st,U,comm_buf,sss,sss,in,out);
      }
    } else { 
PARALLEL_FOR_LOOP
      for(int sss=0;sss<in._grid->oSites();sss++){
        DiracOptDhopSiteDag(st,U,comm_buf,sss,sss,in,out);
      }
    }
  } else {
    if( HandOptDslash ) {
PARALLEL_FOR_LOOP
      for(int sss=0;sss<in._grid->oSites();sss++){
        DiracOptHandDhopSite(st,U,comm_buf,sss,sss,in,out);
      }
    } else { 
PARALLEL_FOR_LOOP
      for(int sss=0;sss<in._grid->oSites();sss++){
        DiracOptDhopSite(st,U,comm_buf,sss,sss,in,out);
      }
    }
  }
}
void WilsonFermion::DhopOE(const LatticeFermion &in, LatticeFermion &out,int dag)
{
  conformable(in._grid,_cbgrid);    // verifies half grid
  conformable(in._grid,out._grid); // drops the cb check

  assert(in.checkerboard==Even);
  out.checkerboard = Odd;

  DhopInternal(StencilEven,UmuOdd,in,out,dag);
}
void WilsonFermion::DhopEO(const LatticeFermion &in, LatticeFermion &out,int dag)
{
  conformable(in._grid,_cbgrid);    // verifies half grid
  conformable(in._grid,out._grid); // drops the cb check

  assert(in.checkerboard==Odd);
  out.checkerboard = Even;

  DhopInternal(StencilOdd,UmuEven,in,out,dag);
}
void WilsonFermion::Dhop(const LatticeFermion &in, LatticeFermion &out,int dag)
{
  conformable(in._grid,_grid); // verifies full grid
  conformable(in._grid,out._grid);

  out.checkerboard = in.checkerboard;

  DhopInternal(Stencil,Umu,in,out,dag);
}

void WilsonFermion::DerivInternal(CartesianStencil & st,LatticeDoubledGaugeField & U,
				  LatticeGaugeField &mat,const LatticeFermion &A,const LatticeFermion &B,int dag)
{
  assert((dag==DaggerNo) ||(dag==DaggerYes));

  WilsonCompressor compressor(dag);

  LatticeColourMatrix tmp(B._grid);
  LatticeFermion Btilde(B._grid);

  st.HaloExchange<vSpinColourVector,vHalfSpinColourVector,WilsonCompressor>(B,comm_buf,compressor);
  
  for(int mu=0;mu<Nd;mu++){
      
    ////////////////////////////////////////////////////////////////////////
    // Flip gamma (1+g)<->(1-g) if dag
    ////////////////////////////////////////////////////////////////////////
    int gamma = mu;
    if ( dag ) gamma+= Nd;

    ////////////////////////
    // Call the single hop
    ////////////////////////
PARALLEL_FOR_LOOP
    for(int sss=0;sss<B._grid->oSites();sss++){
      DiracOptDhopDir(st,U,comm_buf,sss,sss,B,Btilde,mu,gamma);
    }

    //////////////////////////////////////////////////
    // spin trace outer product
    //////////////////////////////////////////////////
    tmp = TraceIndex<SpinIndex>(outerProduct(Btilde,A)); 
    PokeIndex<LorentzIndex>(mat,tmp,mu);

  }
}
void WilsonFermion::DhopDeriv(LatticeGaugeField &mat,const LatticeFermion &U,const LatticeFermion &V,int dag)
{
  conformable(U._grid,_grid);  
  conformable(U._grid,V._grid);
  conformable(U._grid,mat._grid);

  mat.checkerboard = U.checkerboard;

  DerivInternal(Stencil,Umu,mat,U,V,dag);
}
void WilsonFermion::DhopDerivOE(LatticeGaugeField &mat,const LatticeFermion &U,const LatticeFermion &V,int dag)
{
  conformable(U._grid,_cbgrid);  
  conformable(U._grid,V._grid);
  conformable(U._grid,mat._grid);

  assert(V.checkerboard==Even);
  assert(U.checkerboard==Odd);
  mat.checkerboard = Odd;

  DerivInternal(StencilEven,UmuOdd,mat,U,V,dag);
}
void WilsonFermion::DhopDerivEO(LatticeGaugeField &mat,const LatticeFermion &U,const LatticeFermion &V,int dag)
{
  conformable(U._grid,_cbgrid);  
  conformable(U._grid,V._grid);
  conformable(U._grid,mat._grid);

  assert(V.checkerboard==Odd);
  assert(U.checkerboard==Even);
  mat.checkerboard = Even;

  DerivInternal(StencilOdd,UmuEven,mat,U,V,dag);
}

}}



