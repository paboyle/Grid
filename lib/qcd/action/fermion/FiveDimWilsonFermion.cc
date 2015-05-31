#include <Grid.h>

namespace Grid {
namespace QCD {
  
  // S-direction is INNERMOST and takes no part in the parity.
  const std::vector<int> FiveDimWilsonFermion::directions   ({1,2,3,4, 1, 2, 3, 4});
  const std::vector<int> FiveDimWilsonFermion::displacements({1,1,1,1,-1,-1,-1,-1});

  int FiveDimWilsonFermion::HandOptDslash;

  // 5d lattice for DWF.
  FiveDimWilsonFermion::FiveDimWilsonFermion(LatticeGaugeField &_Umu,
					   GridCartesian         &FiveDimGrid,
					   GridRedBlackCartesian &FiveDimRedBlackGrid,
					   GridCartesian         &FourDimGrid,
					   GridRedBlackCartesian &FourDimRedBlackGrid,
					   double _mass) :
  _FiveDimGrid(&FiveDimGrid),
  _FiveDimRedBlackGrid(&FiveDimRedBlackGrid),
  _FourDimGrid(&FourDimGrid),
  _FourDimRedBlackGrid(&FourDimRedBlackGrid),
  Stencil    (_FiveDimGrid,npoint,Even,directions,displacements),
  StencilEven(_FiveDimRedBlackGrid,npoint,Even,directions,displacements), // source is Even
  StencilOdd (_FiveDimRedBlackGrid,npoint,Odd ,directions,displacements), // source is Odd
  mass(_mass),
  Umu(_FourDimGrid),
  UmuEven(_FourDimRedBlackGrid),
  UmuOdd (_FourDimRedBlackGrid),
  Lebesgue(_FourDimGrid),
  LebesgueEvenOdd(_FourDimRedBlackGrid)
{
  // some assertions
  assert(FiveDimGrid._ndimension==5);
  assert(FourDimGrid._ndimension==4);
  
  assert(FiveDimRedBlackGrid._ndimension==5);
  assert(FourDimRedBlackGrid._ndimension==4);

  assert(FiveDimRedBlackGrid._checker_dim==1);

  // Dimension zero of the five-d is the Ls direction
  Ls=FiveDimGrid._fdimensions[0];
  assert(FiveDimRedBlackGrid._fdimensions[0]==Ls);
  assert(FiveDimRedBlackGrid._processors[0] ==1);
  assert(FiveDimRedBlackGrid._simd_layout[0]==1);
  assert(FiveDimGrid._processors[0]         ==1);
  assert(FiveDimGrid._simd_layout[0]        ==1);

  // Other dimensions must match the decomposition of the four-D fields 
  for(int d=0;d<4;d++){
    assert(FourDimRedBlackGrid._fdimensions[d]  ==FourDimGrid._fdimensions[d]);
    assert(FiveDimRedBlackGrid._fdimensions[d+1]==FourDimGrid._fdimensions[d]);

    assert(FourDimRedBlackGrid._processors[d]   ==FourDimGrid._processors[d]);
    assert(FiveDimRedBlackGrid._processors[d+1] ==FourDimGrid._processors[d]);

    assert(FourDimRedBlackGrid._simd_layout[d]  ==FourDimGrid._simd_layout[d]);
    assert(FiveDimRedBlackGrid._simd_layout[d+1]==FourDimGrid._simd_layout[d]);

    assert(FiveDimGrid._fdimensions[d+1]        ==FourDimGrid._fdimensions[d]);
    assert(FiveDimGrid._processors[d+1]         ==FourDimGrid._processors[d]);
    assert(FiveDimGrid._simd_layout[d+1]        ==FourDimGrid._simd_layout[d]);
  }

  // Allocate the required comms buffer
  comm_buf.resize(Stencil._unified_buffer_size); // this is always big enough to contain EO
  
  DoubleStore(Umu,_Umu);
  pickCheckerboard(Even,UmuEven,Umu);
  pickCheckerboard(Odd ,UmuOdd,Umu);
}
void FiveDimWilsonFermion::DoubleStore(LatticeDoubledGaugeField &Uds,const LatticeGaugeField &Umu)
{
  conformable(Uds._grid,GaugeGrid());
  conformable(Umu._grid,GaugeGrid());
  LatticeColourMatrix U(GaugeGrid());
  for(int mu=0;mu<Nd;mu++){
    U = peekIndex<LorentzIndex>(Umu,mu);
    pokeIndex<LorentzIndex>(Uds,U,mu);
    U = adj(Cshift(U,mu,-1));
    pokeIndex<LorentzIndex>(Uds,U,mu+4);
  }
}

RealD FiveDimWilsonFermion::M(const LatticeFermion &in, LatticeFermion &out)
{
  out.checkerboard=in.checkerboard;
  Dhop(in,out,DaggerNo);
  return axpy_norm(out,5.0-M5,in,out);
}
RealD FiveDimWilsonFermion::Mdag(const LatticeFermion &in, LatticeFermion &out)
{
  out.checkerboard=in.checkerboard;
  Dhop(in,out,DaggerYes);
  return axpy_norm(out,5.0-M5,in,out);
}
void FiveDimWilsonFermion::Meooe(const LatticeFermion &in, LatticeFermion &out)
{
  if ( in.checkerboard == Odd ) {
    DhopEO(in,out,DaggerNo);
  } else {
    DhopOE(in,out,DaggerNo);
  }
}
void FiveDimWilsonFermion::MeooeDag(const LatticeFermion &in, LatticeFermion &out)
{
  if ( in.checkerboard == Odd ) {
    DhopEO(in,out,DaggerYes);
  } else {
    DhopOE(in,out,DaggerYes);
  }
}
void FiveDimWilsonFermion::Mooee(const LatticeFermion &in, LatticeFermion &out)
{
  out.checkerboard = in.checkerboard;
  out = (5.0-M5)*in;
  return ;
}
void FiveDimWilsonFermion::MooeeDag(const LatticeFermion &in, LatticeFermion &out)
{
  out.checkerboard = in.checkerboard;
  Mooee(in,out);
}
void FiveDimWilsonFermion::MooeeInv(const LatticeFermion &in, LatticeFermion &out)
{
  out.checkerboard = in.checkerboard;
  out = (1.0/(5.0-M5))*in;
  return ;
}
void FiveDimWilsonFermion::MooeeInvDag(const LatticeFermion &in, LatticeFermion &out)
{
  out.checkerboard = in.checkerboard;
  MooeeInv(in,out);
}
void FiveDimWilsonFermion::DhopInternal(CartesianStencil & st, LebesgueOrder &lo,
					LatticeDoubledGaugeField & U,
					const LatticeFermion &in, LatticeFermion &out,int dag)
{
  assert((dag==DaggerNo) ||(dag==DaggerYes));

  WilsonCompressor compressor(dag);

  st.HaloExchange<vSpinColourVector,vHalfSpinColourVector,WilsonCompressor>(in,comm_buf,compressor);
  
  // Dhop takes the 4d grid from U, and makes a 5d index for fermion
  // Not loop ordering and data layout.
  // Designed to create 
  // - per thread reuse in L1 cache for U
  // - 8 linear access unit stride streams per thread for Fermion for hw prefetchable.
  if ( dag == DaggerYes ) {
    if( HandOptDslash ) {
      for(int ss=0;ss<U._grid->oSites();ss++){
	int sU=lo.Reorder(ss);
PARALLEL_FOR_LOOP
	for(int s=0;s<Ls;s++){
	  int sF = s+Ls*sU;
	  DiracOptHand::DhopSiteDag(st,U,comm_buf,sF,sU,in,out);
	}
      }
    } else { 
      for(int ss=0;ss<U._grid->oSites();ss++){
	int sU=lo.Reorder(ss);
PARALLEL_FOR_LOOP
	for(int s=0;s<Ls;s++){
	  int sF = s+Ls*sU;
	  DiracOpt::DhopSiteDag(st,U,comm_buf,sF,sU,in,out);
	}
      }
    }
  } else {
    if( HandOptDslash ) {

PARALLEL_FOR_LOOP
      for(int ss=0;ss<U._grid->oSites();ss++){
	int sU=lo.Reorder(ss);
	for(int s=0;s<Ls;s++){
	  int sF = s+Ls*sU;
	  DiracOptHand::DhopSite(st,U,comm_buf,sF,sU,in,out);
	}
      }

    } else { 
      for(int ss=0;ss<U._grid->oSites();ss++){
	int sU=lo.Reorder(ss);
PARALLEL_FOR_LOOP
	for(int s=0;s<Ls;s++){
	  int sF = s+Ls*sU; 
	  DiracOpt::DhopSite(st,U,comm_buf,sF,sU,in,out);
	}
      }
    }
  }
}
void FiveDimWilsonFermion::DhopOE(const LatticeFermion &in, LatticeFermion &out,int dag)
{
  conformable(in._grid,FermionRedBlackGrid());    // verifies half grid
  conformable(in._grid,out._grid); // drops the cb check

  assert(in.checkerboard==Even);
  out.checkerboard = Odd;

  DhopInternal(StencilEven,LebesgueEvenOdd,UmuOdd,in,out,dag);
}
void FiveDimWilsonFermion::DhopEO(const LatticeFermion &in, LatticeFermion &out,int dag)
{
  conformable(in._grid,FermionRedBlackGrid());    // verifies half grid
  conformable(in._grid,out._grid); // drops the cb check

  assert(in.checkerboard==Odd);
  out.checkerboard = Even;

  DhopInternal(StencilOdd,LebesgueEvenOdd,UmuEven,in,out,dag);
}
void FiveDimWilsonFermion::Dhop(const LatticeFermion &in, LatticeFermion &out,int dag)
{
  conformable(in._grid,FermionGrid()); // verifies full grid
  conformable(in._grid,out._grid);

  out.checkerboard = in.checkerboard;

  DhopInternal(Stencil,Lebesgue,Umu,in,out,dag);
}

}}



