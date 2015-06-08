#include <Grid.h>

namespace Grid {
namespace QCD {
  
  // S-direction is INNERMOST and takes no part in the parity.
  const std::vector<int> WilsonFermion5D::directions   ({1,2,3,4, 1, 2, 3, 4});
  const std::vector<int> WilsonFermion5D::displacements({1,1,1,1,-1,-1,-1,-1});

  int WilsonFermion5D::HandOptDslash;

  // 5d lattice for DWF.
  WilsonFermion5D::WilsonFermion5D(LatticeGaugeField &_Umu,
					   GridCartesian         &FiveDimGrid,
					   GridRedBlackCartesian &FiveDimRedBlackGrid,
					   GridCartesian         &FourDimGrid,
					   GridRedBlackCartesian &FourDimRedBlackGrid,
					   RealD _M5) :
  _FiveDimGrid(&FiveDimGrid),
  _FiveDimRedBlackGrid(&FiveDimRedBlackGrid),
  _FourDimGrid(&FourDimGrid),
  _FourDimRedBlackGrid(&FourDimRedBlackGrid),
  Stencil    (_FiveDimGrid,npoint,Even,directions,displacements),
  StencilEven(_FiveDimRedBlackGrid,npoint,Even,directions,displacements), // source is Even
  StencilOdd (_FiveDimRedBlackGrid,npoint,Odd ,directions,displacements), // source is Odd
  M5(_M5),
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
void WilsonFermion5D::DoubleStore(LatticeDoubledGaugeField &Uds,const LatticeGaugeField &Umu)
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
void WilsonFermion5D::DhopDir(const LatticeFermion &in, LatticeFermion &out,int dir,int disp)
{
  assert( (disp==1)||(disp==-1) );

  WilsonCompressor compressor(DaggerNo);
  Stencil.HaloExchange<vSpinColourVector,vHalfSpinColourVector,WilsonCompressor>(in,comm_buf,compressor);
  
  int skip = (disp==1) ? 0 : 1;

  int dirdisp = dir+skip*4;

PARALLEL_FOR_LOOP
  for(int ss=0;ss<Umu._grid->oSites();ss++){
    for(int s=0;s<Ls;s++){
      int sU=ss;
      int sF = s+Ls*sU; 
      DiracOpt::DhopDir(Stencil,Umu,comm_buf,sF,sU,in,out,dirdisp);
    }
  }

};

void WilsonFermion5D::DhopInternal(CartesianStencil & st, LebesgueOrder &lo,
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
PARALLEL_FOR_LOOP
      for(int ss=0;ss<U._grid->oSites();ss++){
	for(int s=0;s<Ls;s++){
	  //int sU=lo.Reorder(ss);
	  int sU=ss;
	  int sF = s+Ls*sU;
	  DiracOptHand::DhopSiteDag(st,U,comm_buf,sF,sU,in,out);
	}
      }
    } else { 
PARALLEL_FOR_LOOP
      for(int ss=0;ss<U._grid->oSites();ss++){
	for(int s=0;s<Ls;s++){
	  //	  int sU=lo.Reorder(ss);
	  int sU=ss;
	  int sF = s+Ls*sU;
	  DiracOpt::DhopSiteDag(st,U,comm_buf,sF,sU,in,out);
	}
      }
    }
  } else {
    if( HandOptDslash ) {
PARALLEL_FOR_LOOP
      for(int ss=0;ss<U._grid->oSites();ss++){
	for(int s=0;s<Ls;s++){
	  //	  int sU=lo.Reorder(ss);
	  int sU=ss;
	  int sF = s+Ls*sU;
	  DiracOptHand::DhopSite(st,U,comm_buf,sF,sU,in,out);
	}
      }

    } else { 
PARALLEL_FOR_LOOP
      for(int ss=0;ss<U._grid->oSites();ss++){
	for(int s=0;s<Ls;s++){
	  //	  int sU=lo.Reorder(ss);
	  int sU=ss;
	  int sF = s+Ls*sU; 
	  DiracOpt::DhopSite(st,U,comm_buf,sF,sU,in,out);
	}
      }
    }
  }
}
void WilsonFermion5D::DhopOE(const LatticeFermion &in, LatticeFermion &out,int dag)
{
  conformable(in._grid,FermionRedBlackGrid());    // verifies half grid
  conformable(in._grid,out._grid); // drops the cb check

  assert(in.checkerboard==Even);
  out.checkerboard = Odd;

  DhopInternal(StencilEven,LebesgueEvenOdd,UmuOdd,in,out,dag);
}
void WilsonFermion5D::DhopEO(const LatticeFermion &in, LatticeFermion &out,int dag)
{
  conformable(in._grid,FermionRedBlackGrid());    // verifies half grid
  conformable(in._grid,out._grid); // drops the cb check

  assert(in.checkerboard==Odd);
  out.checkerboard = Even;

  DhopInternal(StencilOdd,LebesgueEvenOdd,UmuEven,in,out,dag);
}
void WilsonFermion5D::Dhop(const LatticeFermion &in, LatticeFermion &out,int dag)
{
  conformable(in._grid,FermionGrid()); // verifies full grid
  conformable(in._grid,out._grid);

  out.checkerboard = in.checkerboard;

  DhopInternal(Stencil,Lebesgue,Umu,in,out,dag);
}
void WilsonFermion5D::DW(const LatticeFermion &in, LatticeFermion &out,int dag)
{
  out.checkerboard=in.checkerboard;
  Dhop(in,out,dag); // -0.5 is included
  axpy(out,4.0-M5,in,out);
}
}
}



