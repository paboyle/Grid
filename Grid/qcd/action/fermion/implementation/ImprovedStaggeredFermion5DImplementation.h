/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/action/fermion/ImprovedStaggeredFermion5D.cc

    Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

    See the full license in the file "LICENSE" in the top level distribution directory
*************************************************************************************/
/*  END LEGAL */
#include <Grid/qcd/action/fermion/FermionCore.h>
#include <Grid/qcd/action/fermion/ImprovedStaggeredFermion5D.h>
#include <Grid/perfmon/PerfCount.h>

#pragma once

NAMESPACE_BEGIN(Grid);

// 5d lattice for DWF.
template<class Impl>
ImprovedStaggeredFermion5D<Impl>::ImprovedStaggeredFermion5D(GridCartesian         &FiveDimGrid,
							     GridRedBlackCartesian &FiveDimRedBlackGrid,
							     GridCartesian         &FourDimGrid,
							     GridRedBlackCartesian &FourDimRedBlackGrid,
							     RealD _mass,
							     RealD _c1,RealD _c2, RealD _u0,
							     const ImplParams &p) :
  Kernels(p),
  _FiveDimGrid        (&FiveDimGrid),
  _FiveDimRedBlackGrid(&FiveDimRedBlackGrid),
  _FourDimGrid        (&FourDimGrid),
  _FourDimRedBlackGrid(&FourDimRedBlackGrid),
  Stencil    (&FiveDimGrid,npoint,Even,directions,displacements,p),
  StencilEven(&FiveDimRedBlackGrid,npoint,Even,directions,displacements,p), // source is Even
  StencilOdd (&FiveDimRedBlackGrid,npoint,Odd ,directions,displacements,p), // source is Odd
  mass(_mass),
  c1(_c1),
  c2(_c2),
  u0(_u0),
  Umu(&FourDimGrid),
  UmuEven(&FourDimRedBlackGrid),
  UmuOdd (&FourDimRedBlackGrid),
  UUUmu(&FourDimGrid),
  UUUmuEven(&FourDimRedBlackGrid),
  UUUmuOdd(&FourDimRedBlackGrid),
  _tmp(&FiveDimRedBlackGrid)
{

  // some assertions
  assert(FiveDimGrid._ndimension==5);
  assert(FourDimGrid._ndimension==4);
  assert(FourDimRedBlackGrid._ndimension==4);
  assert(FiveDimRedBlackGrid._ndimension==5);
  assert(FiveDimRedBlackGrid._checker_dim==1); // Don't checker the s direction

  // extent of fifth dim and not spread out
  Ls=FiveDimGrid._fdimensions[0];
  assert(FiveDimRedBlackGrid._fdimensions[0]==Ls);
  assert(FiveDimGrid._processors[0]         ==1);
  assert(FiveDimRedBlackGrid._processors[0] ==1);

  // Other dimensions must match the decomposition of the four-D fields 
  for(int d=0;d<4;d++){
    assert(FiveDimGrid._processors[d+1]         ==FourDimGrid._processors[d]);
    assert(FiveDimRedBlackGrid._processors[d+1] ==FourDimGrid._processors[d]);
    assert(FourDimRedBlackGrid._processors[d]   ==FourDimGrid._processors[d]);

    assert(FiveDimGrid._fdimensions[d+1]        ==FourDimGrid._fdimensions[d]);
    assert(FiveDimRedBlackGrid._fdimensions[d+1]==FourDimGrid._fdimensions[d]);
    assert(FourDimRedBlackGrid._fdimensions[d]  ==FourDimGrid._fdimensions[d]);

    assert(FiveDimGrid._simd_layout[d+1]        ==FourDimGrid._simd_layout[d]);
    assert(FiveDimRedBlackGrid._simd_layout[d+1]==FourDimGrid._simd_layout[d]);
    assert(FourDimRedBlackGrid._simd_layout[d]  ==FourDimGrid._simd_layout[d]);
  }

  if (Impl::LsVectorised) { 

    int nsimd = Simd::Nsimd();
    
    // Dimension zero of the five-d is the Ls direction
    assert(FiveDimGrid._simd_layout[0]        ==nsimd);
    assert(FiveDimRedBlackGrid._simd_layout[0]==nsimd);

    for(int d=0;d<4;d++){
      assert(FourDimGrid._simd_layout[d]==1);
      assert(FourDimRedBlackGrid._simd_layout[d]==1);
      assert(FiveDimRedBlackGrid._simd_layout[d+1]==1);
    }

  } else {
    
    // Dimension zero of the five-d is the Ls direction
    assert(FiveDimRedBlackGrid._simd_layout[0]==1);
    assert(FiveDimGrid._simd_layout[0]        ==1);

  }
  int LLs = FiveDimGrid._rdimensions[0];
  int vol4= FourDimGrid.oSites();
  Stencil.BuildSurfaceList(LLs,vol4);

  vol4=FourDimRedBlackGrid.oSites();
  StencilEven.BuildSurfaceList(LLs,vol4);
  StencilOdd.BuildSurfaceList(LLs,vol4);
}
template <class Impl>
void ImprovedStaggeredFermion5D<Impl>::CopyGaugeCheckerboards(void)
{
  pickCheckerboard(Even, UmuEven,  Umu);
  pickCheckerboard(Odd,  UmuOdd ,  Umu);
  pickCheckerboard(Even, UUUmuEven,UUUmu);
  pickCheckerboard(Odd,  UUUmuOdd, UUUmu);
}
template<class Impl>
ImprovedStaggeredFermion5D<Impl>::ImprovedStaggeredFermion5D(GaugeField &_Uthin,GaugeField &_Ufat,
							     GridCartesian         &FiveDimGrid,
							     GridRedBlackCartesian &FiveDimRedBlackGrid,
							     GridCartesian         &FourDimGrid,
							     GridRedBlackCartesian &FourDimRedBlackGrid,
							     RealD _mass,
							     RealD _c1,RealD _c2, RealD _u0,
							     const ImplParams &p) :
  ImprovedStaggeredFermion5D(FiveDimGrid,FiveDimRedBlackGrid,
			     FourDimGrid,FourDimRedBlackGrid,
			     _mass,_c1,_c2,_u0,p)
{
  ImportGauge(_Uthin,_Ufat);
}

///////////////////////////////////////////////////
// For MILC use; pass three link U's and 1 link U
///////////////////////////////////////////////////
template <class Impl>
void ImprovedStaggeredFermion5D<Impl>::ImportGaugeSimple(const GaugeField &_Utriple,const GaugeField &_Ufat) 
{
  /////////////////////////////////////////////////////////////////
  // Trivial import; phases and fattening and such like preapplied
  /////////////////////////////////////////////////////////////////
  for (int mu = 0; mu < Nd; mu++) {

    auto U = PeekIndex<LorentzIndex>(_Utriple, mu);
    Impl::InsertGaugeField(UUUmu,U,mu);

    U = adj( Cshift(U, mu, -3));
    Impl::InsertGaugeField(UUUmu,-U,mu+4);

    U = PeekIndex<LorentzIndex>(_Ufat, mu);
    Impl::InsertGaugeField(Umu,U,mu);

    U = adj( Cshift(U, mu, -1));
    Impl::InsertGaugeField(Umu,-U,mu+4);

  }
  CopyGaugeCheckerboards();
}
template <class Impl>
void ImprovedStaggeredFermion5D<Impl>::ImportGaugeSimple(const DoubledGaugeField &_UUU,const DoubledGaugeField &_U) 
{
  /////////////////////////////////////////////////////////////////
  // Trivial import; phases and fattening and such like preapplied
  /////////////////////////////////////////////////////////////////
  Umu   = _U;
  UUUmu = _UUU;
  CopyGaugeCheckerboards();
}
template<class Impl>
void ImprovedStaggeredFermion5D<Impl>::ImportGauge(const GaugeField &_Uthin,const GaugeField &_Ufat)
{
  ////////////////////////////////////////////////////////
  // Double Store should take two fields for Naik and one hop separately.
  ////////////////////////////////////////////////////////
  Impl::DoubleStore(GaugeGrid(), UUUmu, Umu, _Uthin, _Ufat );

  ////////////////////////////////////////////////////////
  // Apply scale factors to get the right fermion Kinetic term
  // Could pass coeffs into the double store to save work.
  // 0.5 ( U p(x+mu) - Udag(x-mu) p(x-mu) ) 
  ////////////////////////////////////////////////////////
  for (int mu = 0; mu < Nd; mu++) {

    auto U = PeekIndex<LorentzIndex>(Umu, mu);
    PokeIndex<LorentzIndex>(Umu, U*( 0.5*c1/u0), mu );
    
    U = PeekIndex<LorentzIndex>(Umu, mu+4);
    PokeIndex<LorentzIndex>(Umu, U*(-0.5*c1/u0), mu+4);

    U = PeekIndex<LorentzIndex>(UUUmu, mu);
    PokeIndex<LorentzIndex>(UUUmu, U*( 0.5*c2/u0/u0/u0), mu );
    
    U = PeekIndex<LorentzIndex>(UUUmu, mu+4);
    PokeIndex<LorentzIndex>(UUUmu, U*(-0.5*c2/u0/u0/u0), mu+4);
  }

  CopyGaugeCheckerboards();
}
template<class Impl>
void ImprovedStaggeredFermion5D<Impl>::DhopDir(const FermionField &in, FermionField &out,int dir5,int disp)
{
  int dir = dir5-1; // Maps to the ordering above in "directions" that is passed to stencil
                    // we drop off the innermost fifth dimension

  Compressor compressor;
  Stencil.HaloExchange(in,compressor);
  autoView( Umu_v   ,   Umu, CpuRead);
  autoView( UUUmu_v , UUUmu, CpuRead);
  autoView( in_v    ,  in, CpuRead);
  autoView( out_v   , out, CpuWrite);
  thread_for( ss,Umu.Grid()->oSites(),{
    for(int s=0;s<Ls;s++){
      int sU=ss;
      int sF = s+Ls*sU; 
      Kernels::DhopDirKernel(Stencil, Umu_v, UUUmu_v, Stencil.CommBuf(), sF, sU, in_v, out_v, dir, disp);
    }
  });
};

template<class Impl>
void ImprovedStaggeredFermion5D<Impl>::DerivInternal(StencilImpl & st,
						     DoubledGaugeField & U,
						     DoubledGaugeField & UUU,
						     GaugeField &mat,
						     const FermionField &A,
						     const FermionField &B,
						     int dag)
{
  // No force terms in multi-rhs solver staggered
  assert(0);
}

template<class Impl>
void ImprovedStaggeredFermion5D<Impl>::DhopDeriv(GaugeField &mat,
						 const FermionField &A,
						 const FermionField &B,
						 int dag)
{
  assert(0);
}

template<class Impl>
void ImprovedStaggeredFermion5D<Impl>::DhopDerivEO(GaugeField &mat,
						   const FermionField &A,
						   const FermionField &B,
						   int dag)
{
  assert(0);
}


template<class Impl>
void ImprovedStaggeredFermion5D<Impl>::DhopDerivOE(GaugeField &mat,
						   const FermionField &A,
						   const FermionField &B,
						   int dag)
{
  assert(0);
}

/*CHANGE */
template<class Impl>
void ImprovedStaggeredFermion5D<Impl>::DhopInternal(StencilImpl & st, 
						    DoubledGaugeField & U,DoubledGaugeField & UUU,
						    const FermionField &in, FermionField &out,int dag)
{
  if ( StaggeredKernelsStatic::Comms == StaggeredKernelsStatic::CommsAndCompute )
    DhopInternalOverlappedComms(st,U,UUU,in,out,dag);
  else
    DhopInternalSerialComms(st,U,UUU,in,out,dag);
}

template<class Impl>
void ImprovedStaggeredFermion5D<Impl>::DhopInternalOverlappedComms(StencilImpl & st, 
								   DoubledGaugeField & U,DoubledGaugeField & UUU,
								   const FermionField &in, FermionField &out,int dag)
{
  //  assert((dag==DaggerNo) ||(dag==DaggerYes));
  Compressor compressor; 

  int LLs = in.Grid()->_rdimensions[0];
  int len =  U.Grid()->oSites();

  st.Prepare();
  st.HaloGather(in,compressor);

  std::vector<std::vector<CommsRequest_t> > requests;
  st.CommunicateBegin(requests);

  //  st.HaloExchangeOptGather(in,compressor); // Wilson compressor
  st.CommsMergeSHM(compressor);// Could do this inside parallel region overlapped with comms

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  // Remove explicit thread mapping introduced for OPA reasons.
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  {
    int interior=1;
    int exterior=0;
    Kernels::DhopImproved(st,U,UUU,in,out,dag,interior,exterior);
  }

  st.CommsMerge(compressor);

  st.CommunicateComplete(requests);

  {
    int interior=0;
    int exterior=1;
    Kernels::DhopImproved(st,U,UUU,in,out,dag,interior,exterior);
  }
}

template<class Impl>
void ImprovedStaggeredFermion5D<Impl>::DhopInternalSerialComms(StencilImpl & st, 
						    DoubledGaugeField & U,DoubledGaugeField & UUU,
						    const FermionField &in, FermionField &out,int dag)
{
  Compressor compressor;
  int LLs = in.Grid()->_rdimensions[0];

  st.HaloExchange(in,compressor);
  
  // Dhop takes the 4d grid from U, and makes a 5d index for fermion
  {
    int interior=1;
    int exterior=1;
    Kernels::DhopImproved(st,U,UUU,in,out,dag,interior,exterior);
  }
}
/*CHANGE END*/



template<class Impl>
void ImprovedStaggeredFermion5D<Impl>::DhopOE(const FermionField &in, FermionField &out,int dag)
{
  conformable(in.Grid(),FermionRedBlackGrid());    // verifies half grid
  conformable(in.Grid(),out.Grid()); // drops the cb check

  assert(in.Checkerboard()==Even);
  out.Checkerboard() = Odd;

  DhopInternal(StencilEven,UmuOdd,UUUmuOdd,in,out,dag);
}
template<class Impl>
void ImprovedStaggeredFermion5D<Impl>::DhopEO(const FermionField &in, FermionField &out,int dag)
{
  conformable(in.Grid(),FermionRedBlackGrid());    // verifies half grid
  conformable(in.Grid(),out.Grid()); // drops the cb check

  assert(in.Checkerboard()==Odd);
  out.Checkerboard() = Even;

  DhopInternal(StencilOdd,UmuEven,UUUmuEven,in,out,dag);
}
template<class Impl>
void ImprovedStaggeredFermion5D<Impl>::Dhop(const FermionField &in, FermionField &out,int dag)
{
  conformable(in.Grid(),FermionGrid()); // verifies full grid
  conformable(in.Grid(),out.Grid());

  out.Checkerboard() = in.Checkerboard();

  DhopInternal(Stencil,Umu,UUUmu,in,out,dag);
}

/////////////////////////////////////////////////////////////////////////
// Implement the general interface. Here we use SAME mass on all slices
/////////////////////////////////////////////////////////////////////////
template <class Impl>
void ImprovedStaggeredFermion5D<Impl>::Mdir(const FermionField &in, FermionField &out, int dir, int disp) 
{
  DhopDir(in, out, dir, disp);
}
template <class Impl>
void ImprovedStaggeredFermion5D<Impl>::MdirAll(const FermionField &in, std::vector<FermionField> &out) 
{
  assert(0);
}
template <class Impl>
void ImprovedStaggeredFermion5D<Impl>::M(const FermionField &in, FermionField &out) 
{
  out.Checkerboard() = in.Checkerboard();
  Dhop(in, out, DaggerNo);
  axpy(out, mass, in, out);
}

template <class Impl>
void ImprovedStaggeredFermion5D<Impl>::Mdag(const FermionField &in, FermionField &out) 
{
  out.Checkerboard() = in.Checkerboard();
  Dhop(in, out, DaggerYes);
  axpy(out, mass, in, out);
}

template <class Impl>
void ImprovedStaggeredFermion5D<Impl>::Meooe(const FermionField &in, FermionField &out) 
{
  if (in.Checkerboard() == Odd) {
    DhopEO(in, out, DaggerNo);
  } else {
    DhopOE(in, out, DaggerNo);
  }
}
template <class Impl>
void ImprovedStaggeredFermion5D<Impl>::MeooeDag(const FermionField &in, FermionField &out) 
{
  if (in.Checkerboard() == Odd) {
    DhopEO(in, out, DaggerYes);
  } else {
    DhopOE(in, out, DaggerYes);
  }
}

template <class Impl>
void ImprovedStaggeredFermion5D<Impl>::Mooee(const FermionField &in, FermionField &out) 
{
  out.Checkerboard() = in.Checkerboard();
  typename FermionField::scalar_type scal(mass);
  out = scal * in;
}

template <class Impl>
void ImprovedStaggeredFermion5D<Impl>::MooeeDag(const FermionField &in, FermionField &out) 
{
  out.Checkerboard() = in.Checkerboard();
  Mooee(in, out);
}

template <class Impl>
void ImprovedStaggeredFermion5D<Impl>::MooeeInv(const FermionField &in, FermionField &out) 
{
  out.Checkerboard() = in.Checkerboard();
  out = (1.0 / (mass)) * in;
}

template <class Impl>
void ImprovedStaggeredFermion5D<Impl>::MooeeInvDag(const FermionField &in,FermionField &out) 
{
  out.Checkerboard() = in.Checkerboard();
  MooeeInv(in, out);
}

//////////////////////////////////////////////////////// 
// Conserved current - not yet implemented.
////////////////////////////////////////////////////////
template <class Impl>
void ImprovedStaggeredFermion5D<Impl>::ContractConservedCurrent(PropagatorField &q_in_1,
								PropagatorField &q_in_2,
								PropagatorField &q_out,
								PropagatorField &src,
								Current curr_type,
								unsigned int mu)
{
  assert(0);
}

template <class Impl>
void ImprovedStaggeredFermion5D<Impl>::SeqConservedCurrent(PropagatorField &q_in,
							   PropagatorField &q_out,
							   PropagatorField &src,
							   Current curr_type,
							   unsigned int mu, 
							   unsigned int tmin,
							   unsigned int tmax,
							   ComplexField &lattice_cmplx)
{
  assert(0);

}
  
NAMESPACE_END(Grid);



