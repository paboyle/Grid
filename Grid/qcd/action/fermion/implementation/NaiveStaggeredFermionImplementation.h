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
NaiveStaggeredFermion<Impl>::NaiveStaggeredFermion(GridCartesian &Fgrid, GridRedBlackCartesian &Hgrid, 
						   RealD _mass,
						   RealD _c1, RealD _u0,
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
    _tmp(&Hgrid)
{
  int vol4;
  int LLs=1;
  c1=_c1;
  u0=_u0;
  vol4= _grid->oSites();
  Stencil.BuildSurfaceList(LLs,vol4);
  vol4= _cbgrid->oSites();
  StencilEven.BuildSurfaceList(LLs,vol4);
  StencilOdd.BuildSurfaceList(LLs,vol4);
}

template <class Impl>
NaiveStaggeredFermion<Impl>::NaiveStaggeredFermion(GaugeField &_U, GridCartesian &Fgrid,
						   GridRedBlackCartesian &Hgrid, RealD _mass,
						   RealD _c1, RealD _u0,
						   const ImplParams &p)
  : NaiveStaggeredFermion(Fgrid,Hgrid,_mass,_c1,_u0,p)
{
  ImportGauge(_U);
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
void NaiveStaggeredFermion<Impl>::CopyGaugeCheckerboards(void)
{
  pickCheckerboard(Even, UmuEven,  Umu);
  pickCheckerboard(Odd,  UmuOdd ,  Umu);
}
template <class Impl>
void NaiveStaggeredFermion<Impl>::ImportGauge(const GaugeField &_U) 
{
  GaugeLinkField U(GaugeGrid());
  DoubledGaugeField _UUU(GaugeGrid());
  ////////////////////////////////////////////////////////
  // Double Store should take two fields for Naik and one hop separately.
  // Discard teh Naik as Naive
  ////////////////////////////////////////////////////////
  Impl::DoubleStore(GaugeGrid(), _UUU, Umu, _U, _U );

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

  }

  CopyGaugeCheckerboards();
}

/////////////////////////////
// Implement the interface
/////////////////////////////

template <class Impl>
void NaiveStaggeredFermion<Impl>::M(const FermionField &in, FermionField &out) {
  out.Checkerboard() = in.Checkerboard();
  Dhop(in, out, DaggerNo);
  axpy(out, mass, in, out);
}

template <class Impl>
void NaiveStaggeredFermion<Impl>::Mdag(const FermionField &in, FermionField &out) {
  out.Checkerboard() = in.Checkerboard();
  Dhop(in, out, DaggerYes);
  axpy(out, mass, in, out);
}

template <class Impl>
void NaiveStaggeredFermion<Impl>::Meooe(const FermionField &in, FermionField &out) {
  if (in.Checkerboard() == Odd) {
    DhopEO(in, out, DaggerNo);
  } else {
    DhopOE(in, out, DaggerNo);
  }
}
template <class Impl>
void NaiveStaggeredFermion<Impl>::MeooeDag(const FermionField &in, FermionField &out) {
  if (in.Checkerboard() == Odd) {
    DhopEO(in, out, DaggerYes);
  } else {
    DhopOE(in, out, DaggerYes);
  }
}

template <class Impl>
void NaiveStaggeredFermion<Impl>::Mooee(const FermionField &in, FermionField &out) {
  out.Checkerboard() = in.Checkerboard();
  typename FermionField::scalar_type scal(mass);
  out = scal * in;
}

template <class Impl>
void NaiveStaggeredFermion<Impl>::MooeeDag(const FermionField &in, FermionField &out) {
  out.Checkerboard() = in.Checkerboard();
  Mooee(in, out);
}

template <class Impl>
void NaiveStaggeredFermion<Impl>::MooeeInv(const FermionField &in, FermionField &out) {
  out.Checkerboard() = in.Checkerboard();
  out = (1.0 / (mass)) * in;
}

template <class Impl>
void NaiveStaggeredFermion<Impl>::MooeeInvDag(const FermionField &in, FermionField &out) 
{
  out.Checkerboard() = in.Checkerboard();
  MooeeInv(in, out);
}

///////////////////////////////////
// Internal
///////////////////////////////////

template <class Impl>
void NaiveStaggeredFermion<Impl>::DerivInternal(StencilImpl &st, DoubledGaugeField &U,
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
    autoView( U_v      , U, CpuRead);
    autoView( B_v      , B, CpuWrite);
    autoView( Btilde_v , Btilde, CpuWrite);
    thread_for(sss,B.Grid()->oSites(),{
      Kernels::DhopDirKernel(st, U_v, U_v, st.CommBuf(), sss, sss, B_v, Btilde_v, mu,1);
    });

    assert(0);// need to figure out the force interface with a blasted three link term.
    
  }
}

template <class Impl>
void NaiveStaggeredFermion<Impl>::DhopDeriv(GaugeField &mat, const FermionField &U, const FermionField &V, int dag) {

  conformable(U.Grid(), _grid);
  conformable(U.Grid(), V.Grid());
  conformable(U.Grid(), mat.Grid());

  mat.Checkerboard() = U.Checkerboard();

  DerivInternal(Stencil, Umu, mat, U, V, dag);
}

template <class Impl>
void NaiveStaggeredFermion<Impl>::DhopDerivOE(GaugeField &mat, const FermionField &U, const FermionField &V, int dag) {

  conformable(U.Grid(), _cbgrid);
  conformable(U.Grid(), V.Grid());
  conformable(U.Grid(), mat.Grid());

  assert(V.Checkerboard() == Even);
  assert(U.Checkerboard() == Odd);
  mat.Checkerboard() = Odd;

  DerivInternal(StencilEven, UmuOdd, mat, U, V, dag);
}

template <class Impl>
void NaiveStaggeredFermion<Impl>::DhopDerivEO(GaugeField &mat, const FermionField &U, const FermionField &V, int dag) {

  conformable(U.Grid(), _cbgrid);
  conformable(U.Grid(), V.Grid());
  conformable(U.Grid(), mat.Grid());

  assert(V.Checkerboard() == Odd);
  assert(U.Checkerboard() == Even);
  mat.Checkerboard() = Even;

  DerivInternal(StencilOdd, UmuEven, mat, U, V, dag);
}

template <class Impl>
void NaiveStaggeredFermion<Impl>::Dhop(const FermionField &in, FermionField &out, int dag) 
{
  conformable(in.Grid(), _grid);  // verifies full grid
  conformable(in.Grid(), out.Grid());

  out.Checkerboard() = in.Checkerboard();

  DhopInternal(Stencil, Umu, in, out, dag);
}

template <class Impl>
void NaiveStaggeredFermion<Impl>::DhopOE(const FermionField &in, FermionField &out, int dag) 
{
  conformable(in.Grid(), _cbgrid);    // verifies half grid
  conformable(in.Grid(), out.Grid());  // drops the cb check

  assert(in.Checkerboard() == Even);
  out.Checkerboard() = Odd;

  DhopInternal(StencilEven, UmuOdd, in, out, dag);
}

template <class Impl>
void NaiveStaggeredFermion<Impl>::DhopEO(const FermionField &in, FermionField &out, int dag) 
{
  conformable(in.Grid(), _cbgrid);    // verifies half grid
  conformable(in.Grid(), out.Grid());  // drops the cb check

  assert(in.Checkerboard() == Odd);
  out.Checkerboard() = Even;

  DhopInternal(StencilOdd, UmuEven, in, out, dag);
}

template <class Impl>
void NaiveStaggeredFermion<Impl>::Mdir(const FermionField &in, FermionField &out, int dir, int disp) 
{
  DhopDir(in, out, dir, disp);
}
template <class Impl>
void NaiveStaggeredFermion<Impl>::MdirAll(const FermionField &in, std::vector<FermionField> &out) 
{
  assert(0); // Not implemented yet
}

template <class Impl>
void NaiveStaggeredFermion<Impl>::DhopDir(const FermionField &in, FermionField &out, int dir, int disp) 
{

  Compressor compressor;
  Stencil.HaloExchange(in, compressor);
  autoView( Umu_v   ,  Umu, CpuRead);
  autoView( in_v    ,  in, CpuRead);
  autoView( out_v   , out, CpuWrite);
  //  thread_for( sss, in.Grid()->oSites(),{
  //    Kernels::DhopDirKernel(Stencil, Umu_v, Stencil.CommBuf(), sss, sss, in_v, out_v, dir, disp);
  //  });
  assert(0);
};


template <class Impl>
void NaiveStaggeredFermion<Impl>::DhopInternal(StencilImpl &st,
					       DoubledGaugeField &U,
					       const FermionField &in,
					       FermionField &out, int dag) 
{
  if ( StaggeredKernelsStatic::Comms == StaggeredKernelsStatic::CommsAndCompute )
    DhopInternalOverlappedComms(st,U,in,out,dag);
  else
    DhopInternalSerialComms(st,U,in,out,dag);
}
template <class Impl>
void NaiveStaggeredFermion<Impl>::DhopInternalOverlappedComms(StencilImpl &st,
							      DoubledGaugeField &U,
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
    Kernels::DhopNaive(st,U,in,out,dag,interior,exterior);
  }

  st.CommunicateComplete(requests);

  // First to enter, last to leave timing
  st.CommsMerge(compressor);

  {
    int interior=0;
    int exterior=1;
    Kernels::DhopNaive(st,U,in,out,dag,interior,exterior);
  }
}

template <class Impl>
void NaiveStaggeredFermion<Impl>::DhopInternalSerialComms(StencilImpl &st,
							  DoubledGaugeField &U,
							  const FermionField &in,
							  FermionField &out, int dag) 
{
  assert((dag == DaggerNo) || (dag == DaggerYes));

  Compressor compressor;
  st.HaloExchange(in, compressor);

  {
    int interior=1;
    int exterior=1;
    Kernels::DhopNaive(st,U,in,out,dag,interior,exterior);
  }
};

//////////////////////////////////////////////////////// 
// Conserved current - not yet implemented.
////////////////////////////////////////////////////////
template <class Impl>
void NaiveStaggeredFermion<Impl>::ContractConservedCurrent(PropagatorField &q_in_1,
							      PropagatorField &q_in_2,
							      PropagatorField &q_out,
							      PropagatorField &src,
							      Current curr_type,
							      unsigned int mu)
{
  assert(0);
}

template <class Impl>
void NaiveStaggeredFermion<Impl>::SeqConservedCurrent(PropagatorField &q_in,
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
