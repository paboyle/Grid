/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/action/fermion/ImprovedStaggeredFermion.cc

Copyright (C) 2015

Author: Azusa Yamaguchi, Peter Boyle

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

See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
/*  END LEGAL */
#include <Grid/Grid.h>

#pragma once 

NAMESPACE_BEGIN(Grid);

/////////////////////////////////
// Constructor and gauge import
/////////////////////////////////

template <class Impl>
ImprovedStaggeredFermion<Impl>::ImprovedStaggeredFermion(GridCartesian &Fgrid, GridRedBlackCartesian &Hgrid, 
							 RealD _mass,
							 RealD _c1, RealD _c2,RealD _u0,
							 const ImplParams &p)
  : Kernels(p),
    _grid(&Fgrid),
    _cbgrid(&Hgrid),
    Stencil(&Fgrid, npoint, Even, directions, displacements,p),
    StencilEven(&Hgrid, npoint, Even, directions, displacements,p),  // source is Even
    StencilOdd(&Hgrid, npoint, Odd, directions, displacements,p),  // source is Odd
    mass(_mass),
    Umu(&Fgrid),
    UmuEven(&Hgrid),
    UmuOdd(&Hgrid),
    UUUmu(&Fgrid),
    UUUmuEven(&Hgrid),
    UUUmuOdd(&Hgrid) ,
    _tmp(&Hgrid)
{
  int vol4;
  int LLs=1;
  c1=_c1;
  c2=_c2;
  u0=_u0;
  vol4= _grid->oSites();
  Stencil.BuildSurfaceList(LLs,vol4);
  vol4= _cbgrid->oSites();
  StencilEven.BuildSurfaceList(LLs,vol4);
  StencilOdd.BuildSurfaceList(LLs,vol4);
}

template <class Impl>
ImprovedStaggeredFermion<Impl>::ImprovedStaggeredFermion(GaugeField &_Uthin, GaugeField &_Ufat, GridCartesian &Fgrid,
							 GridRedBlackCartesian &Hgrid, RealD _mass,
							 RealD _c1, RealD _c2,RealD _u0,
							 const ImplParams &p)
  : ImprovedStaggeredFermion(Fgrid,Hgrid,_mass,_c1,_c2,_u0,p)
{
  ImportGauge(_Uthin,_Ufat);
}

////////////////////////////////////////////////////////////
// Momentum space propagator should be 
// https://arxiv.org/pdf/hep-lat/9712010.pdf
//
// mom space action.
//   gamma_mu i ( c1 sin pmu + c2 sin 3 pmu ) + m
//
// must track through staggered flavour/spin reduction in literature to 
// turn to free propagator for the one component chi field, a la page 4/5
// of above link to implmement fourier based solver.
////////////////////////////////////////////////////////////
template <class Impl>
void ImprovedStaggeredFermion<Impl>::ImportGaugeSimple(const GaugeField &_Utriple,const GaugeField &_Ufat) 
{
  /////////////////////////////////////////////////////////////////
  // Trivial import; phases and fattening and such like preapplied
  /////////////////////////////////////////////////////////////////
  GaugeLinkField U(GaugeGrid());

  for (int mu = 0; mu < Nd; mu++) {

    U = PeekIndex<LorentzIndex>(_Utriple, mu);
    PokeIndex<LorentzIndex>(UUUmu, U, mu );

    U = adj( Cshift(U, mu, -3));
    PokeIndex<LorentzIndex>(UUUmu, -U, mu+4 );

    U = PeekIndex<LorentzIndex>(_Ufat, mu);
    PokeIndex<LorentzIndex>(Umu, U, mu);

    U = adj( Cshift(U, mu, -1));
    PokeIndex<LorentzIndex>(Umu, -U, mu+4);

  }
  CopyGaugeCheckerboards();
}
template <class Impl>
void ImprovedStaggeredFermion<Impl>::ImportGaugeSimple(const DoubledGaugeField &_UUU,const DoubledGaugeField &_U) 
{

  Umu   = _U;
  UUUmu = _UUU;
  CopyGaugeCheckerboards();
}

template <class Impl>
void ImprovedStaggeredFermion<Impl>::CopyGaugeCheckerboards(void)
{
  pickCheckerboard(Even, UmuEven,  Umu);
  pickCheckerboard(Odd,  UmuOdd ,  Umu);
  pickCheckerboard(Even, UUUmuEven,UUUmu);
  pickCheckerboard(Odd,  UUUmuOdd, UUUmu);
}
template <class Impl>
void ImprovedStaggeredFermion<Impl>::ImportGauge(const GaugeField &_Uthin,const GaugeField &_Ufat) 
{
  GaugeLinkField U(GaugeGrid());

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

    U = PeekIndex<LorentzIndex>(Umu, mu);
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

/////////////////////////////
// Implement the interface
/////////////////////////////

template <class Impl>
void ImprovedStaggeredFermion<Impl>::M(const FermionField &in, FermionField &out) 
{
  out.Checkerboard() = in.Checkerboard();
  Dhop(in, out, DaggerNo);
  axpy(out, mass, in, out);
}

template <class Impl>
void ImprovedStaggeredFermion<Impl>::Mdag(const FermionField &in, FermionField &out) 
{
  out.Checkerboard() = in.Checkerboard();
  Dhop(in, out, DaggerYes);
  axpy(out, mass, in, out);
}

template <class Impl>
void ImprovedStaggeredFermion<Impl>::Meooe(const FermionField &in, FermionField &out) 
{
  if (in.Checkerboard() == Odd) {
    DhopEO(in, out, DaggerNo);
  } else {
    DhopOE(in, out, DaggerNo);
  }
}
template <class Impl>
void ImprovedStaggeredFermion<Impl>::MeooeDag(const FermionField &in, FermionField &out) 
{
  if (in.Checkerboard() == Odd) {
    DhopEO(in, out, DaggerYes);
  } else {
    DhopOE(in, out, DaggerYes);
  }
}

template <class Impl>
void ImprovedStaggeredFermion<Impl>::Mooee(const FermionField &in, FermionField &out) 
{
  out.Checkerboard() = in.Checkerboard();
  typename FermionField::scalar_type scal(mass);
  out = scal * in;
}

template <class Impl>
void ImprovedStaggeredFermion<Impl>::MooeeDag(const FermionField &in, FermionField &out) 
{
  out.Checkerboard() = in.Checkerboard();
  Mooee(in, out);
}

template <class Impl>
void ImprovedStaggeredFermion<Impl>::MooeeInv(const FermionField &in, FermionField &out) 
{
  out.Checkerboard() = in.Checkerboard();
  out = (1.0 / (mass)) * in;
}

template <class Impl>
void ImprovedStaggeredFermion<Impl>::MooeeInvDag(const FermionField &in,FermionField &out) 
{
  out.Checkerboard() = in.Checkerboard();
  MooeeInv(in, out);
}

///////////////////////////////////
// Internal
///////////////////////////////////

template <class Impl>
void ImprovedStaggeredFermion<Impl>::DerivInternal(StencilImpl &st, DoubledGaugeField &U, DoubledGaugeField &UUU, 
						   GaugeField & mat,
						   const FermionField &A, const FermionField &B, int dag) 
{
  assert((dag == DaggerNo) || (dag == DaggerYes));

  Compressor compressor;

  FermionField Btilde(B.Grid());
  FermionField Atilde(B.Grid());
  Atilde = A;

  st.HaloExchange(B, compressor);

  for (int mu = 0; mu < Nd; mu++) {

    ////////////////////////
    // Call the single hop
    ////////////////////////
    autoView( U_v   , U, CpuRead);
    autoView( UUU_v , UUU, CpuRead);
    autoView( B_v      , B, CpuWrite);
    autoView( Btilde_v , Btilde, CpuWrite);
    thread_for(sss,B.Grid()->oSites(),{
      Kernels::DhopDirKernel(st, U_v, UUU_v, st.CommBuf(), sss, sss, B_v, Btilde_v, mu,1);
    });

    // Force in three link terms
    //
    //    Impl::InsertForce4D(mat, Btilde, Atilde, mu);
    //
    // dU_ac(x)/dt = i p_ab U_bc(x)
    //
    // => dS_f/dt = dS_f/dU_ac(x) . dU_ac(x)/dt =  i p_ab U_bc(x) dS_f/dU_ac(x) 
    //
    // One link: form fragments S_f = A U B 
    //
    //         write Btilde = U(x) B(x+mu)
    //
    // mat+= TraceIndex<SpinIndex>(outerProduct(Btilde,A)); 
    // 
    // Three link: form fragments S_f = A UUU B 
    //
    // mat+= outer ( A, UUUB) <-- Best take DhopDeriv with one linke or identity matrix
    // mat+= outer ( AU, UUB) <-- and then use covariant cshift?
    // mat+= outer ( AUU, UB) <-- Returned from call to DhopDir

    assert(0);// need to figure out the force interface with a blasted three link term.
    
  }
}

template <class Impl>
void ImprovedStaggeredFermion<Impl>::DhopDeriv(GaugeField &mat, const FermionField &U, const FermionField &V, int dag) 
{
  conformable(U.Grid(), _grid);
  conformable(U.Grid(), V.Grid());
  conformable(U.Grid(), mat.Grid());

  mat.Checkerboard() = U.Checkerboard();

  DerivInternal(Stencil, Umu, UUUmu, mat, U, V, dag);
}

template <class Impl>
void ImprovedStaggeredFermion<Impl>::DhopDerivOE(GaugeField &mat, const FermionField &U, const FermionField &V, int dag) 
{
  conformable(U.Grid(), _cbgrid);
  conformable(U.Grid(), V.Grid());
  conformable(U.Grid(), mat.Grid());

  assert(V.Checkerboard() == Even);
  assert(U.Checkerboard() == Odd);
  mat.Checkerboard() = Odd;

  DerivInternal(StencilEven, UmuOdd, UUUmuOdd, mat, U, V, dag);
}

template <class Impl>
void ImprovedStaggeredFermion<Impl>::DhopDerivEO(GaugeField &mat, const FermionField &U, const FermionField &V, int dag) 
{
  conformable(U.Grid(), _cbgrid);
  conformable(U.Grid(), V.Grid());
  conformable(U.Grid(), mat.Grid());

  assert(V.Checkerboard() == Odd);
  assert(U.Checkerboard() == Even);
  mat.Checkerboard() = Even;

  DerivInternal(StencilOdd, UmuEven, UUUmuEven, mat, U, V, dag);
}

template <class Impl>
void ImprovedStaggeredFermion<Impl>::Dhop(const FermionField &in, FermionField &out, int dag) 
{
  conformable(in.Grid(), _grid);  // verifies full grid
  conformable(in.Grid(), out.Grid());

  out.Checkerboard() = in.Checkerboard();

  DhopInternal(Stencil, Umu, UUUmu, in, out, dag);
}

template <class Impl>
void ImprovedStaggeredFermion<Impl>::DhopOE(const FermionField &in, FermionField &out, int dag) 
{
  conformable(in.Grid(), _cbgrid);    // verifies half grid
  conformable(in.Grid(), out.Grid());  // drops the cb check

  assert(in.Checkerboard() == Even);
  out.Checkerboard() = Odd;

  DhopInternal(StencilEven, UmuOdd, UUUmuOdd, in, out, dag);
}

template <class Impl>
void ImprovedStaggeredFermion<Impl>::DhopEO(const FermionField &in, FermionField &out, int dag) 
{
  conformable(in.Grid(), _cbgrid);    // verifies half grid
  conformable(in.Grid(), out.Grid());  // drops the cb check

  assert(in.Checkerboard() == Odd);
  out.Checkerboard() = Even;

  DhopInternal(StencilOdd, UmuEven, UUUmuEven, in, out, dag);
}

template <class Impl>
void ImprovedStaggeredFermion<Impl>::Mdir(const FermionField &in, FermionField &out, int dir, int disp) 
{
  DhopDir(in, out, dir, disp);
}
template <class Impl>
void ImprovedStaggeredFermion<Impl>::MdirAll(const FermionField &in, std::vector<FermionField> &out) 
{
  assert(0); // Not implemented yet
}

template <class Impl>
void ImprovedStaggeredFermion<Impl>::DhopDir(const FermionField &in, FermionField &out, int dir, int disp) 
{

  Compressor compressor;
  Stencil.HaloExchange(in, compressor);
  autoView( Umu_v   ,   Umu, CpuRead);
  autoView( UUUmu_v , UUUmu, CpuRead);
  autoView( in_v    ,  in, CpuRead);
  autoView( out_v   , out, CpuWrite);
  thread_for( sss, in.Grid()->oSites(),{
    Kernels::DhopDirKernel(Stencil, Umu_v, UUUmu_v, Stencil.CommBuf(), sss, sss, in_v, out_v, dir, disp);
  });
};


template <class Impl>
void ImprovedStaggeredFermion<Impl>::DhopInternal(StencilImpl &st, 
						  DoubledGaugeField &U,
						  DoubledGaugeField &UUU,
						  const FermionField &in,
						  FermionField &out, int dag) 
{
  if ( StaggeredKernelsStatic::Comms == StaggeredKernelsStatic::CommsAndCompute )
    DhopInternalOverlappedComms(st,U,UUU,in,out,dag);
  else
    DhopInternalSerialComms(st,U,UUU,in,out,dag);
}
template <class Impl>
void ImprovedStaggeredFermion<Impl>::DhopInternalOverlappedComms(StencilImpl &st, 
								 DoubledGaugeField &U,
								 DoubledGaugeField &UUU,
								 const FermionField &in,
								 FermionField &out, int dag) 
{
  Compressor compressor; 
  int len =  U.Grid()->oSites();

  st.Prepare();
  st.HaloGather(in,compressor);

  std::vector<std::vector<CommsRequest_t> > requests;
  st.CommunicateBegin(requests);

  st.CommsMergeSHM(compressor);

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  // Removed explicit thread comms
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  {
    int interior=1;
    int exterior=0;
    Kernels::DhopImproved(st,U,UUU,in,out,dag,interior,exterior);
  }

  st.CommunicateComplete(requests);

  // First to enter, last to leave timing
  st.CommsMerge(compressor);

  {
    int interior=0;
    int exterior=1;
    Kernels::DhopImproved(st,U,UUU,in,out,dag,interior,exterior);
  }
}


template <class Impl>
void ImprovedStaggeredFermion<Impl>::DhopInternalSerialComms(StencilImpl &st, 
							     DoubledGaugeField &U,
							     DoubledGaugeField &UUU,
							     const FermionField &in,
							     FermionField &out, int dag) 
{
  assert((dag == DaggerNo) || (dag == DaggerYes));

  Compressor compressor;
  st.HaloExchange(in, compressor);

  {
    int interior=1;
    int exterior=1;
    Kernels::DhopImproved(st,U,UUU,in,out,dag,interior,exterior);
  }
};

//////////////////////////////////////////////////////// 
// Conserved current - not yet implemented.
////////////////////////////////////////////////////////
template <class Impl>
void ImprovedStaggeredFermion<Impl>::ContractConservedCurrent(PropagatorField &q_in_1,
							      PropagatorField &q_in_2,
							      PropagatorField &q_out,
							      PropagatorField &src,
							      Current curr_type,
							      unsigned int mu)
{
  assert(0);
}

template <class Impl>
void ImprovedStaggeredFermion<Impl>::SeqConservedCurrent(PropagatorField &q_in,
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
