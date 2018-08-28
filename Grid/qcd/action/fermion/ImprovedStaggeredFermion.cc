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
#include <Grid.h>

namespace Grid {
namespace QCD {

const std::vector<int> 
ImprovedStaggeredFermionStatic::directions({0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3});
const std::vector<int> 
ImprovedStaggeredFermionStatic::displacements({1, 1, 1, 1, -1, -1, -1, -1, 3, 3, 3, 3, -3, -3, -3, -3});

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
      Stencil(&Fgrid, npoint, Even, directions, displacements),
      StencilEven(&Hgrid, npoint, Even, directions, displacements),  // source is Even
      StencilOdd(&Hgrid, npoint, Odd, directions, displacements),  // source is Odd
      mass(_mass),
      Lebesgue(_grid),
      LebesgueEvenOdd(_cbgrid),
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
RealD ImprovedStaggeredFermion<Impl>::M(const FermionField &in, FermionField &out) {
  out.checkerboard = in.checkerboard;
  Dhop(in, out, DaggerNo);
  return axpy_norm(out, mass, in, out);
}

template <class Impl>
RealD ImprovedStaggeredFermion<Impl>::Mdag(const FermionField &in, FermionField &out) {
  out.checkerboard = in.checkerboard;
  Dhop(in, out, DaggerYes);
  return axpy_norm(out, mass, in, out);
}

template <class Impl>
void ImprovedStaggeredFermion<Impl>::Meooe(const FermionField &in, FermionField &out) {
  if (in.checkerboard == Odd) {
    DhopEO(in, out, DaggerNo);
  } else {
    DhopOE(in, out, DaggerNo);
  }
}
template <class Impl>
void ImprovedStaggeredFermion<Impl>::MeooeDag(const FermionField &in, FermionField &out) {
  if (in.checkerboard == Odd) {
    DhopEO(in, out, DaggerYes);
  } else {
    DhopOE(in, out, DaggerYes);
  }
}

template <class Impl>
void ImprovedStaggeredFermion<Impl>::Mooee(const FermionField &in, FermionField &out) {
  out.checkerboard = in.checkerboard;
  typename FermionField::scalar_type scal(mass);
  out = scal * in;
}

template <class Impl>
void ImprovedStaggeredFermion<Impl>::MooeeDag(const FermionField &in, FermionField &out) {
  out.checkerboard = in.checkerboard;
  Mooee(in, out);
}

template <class Impl>
void ImprovedStaggeredFermion<Impl>::MooeeInv(const FermionField &in, FermionField &out) {
  out.checkerboard = in.checkerboard;
  out = (1.0 / (mass)) * in;
}

template <class Impl>
void ImprovedStaggeredFermion<Impl>::MooeeInvDag(const FermionField &in,
                                      FermionField &out) {
  out.checkerboard = in.checkerboard;
  MooeeInv(in, out);
}

///////////////////////////////////
// Internal
///////////////////////////////////

template <class Impl>
void ImprovedStaggeredFermion<Impl>::DerivInternal(StencilImpl &st, DoubledGaugeField &U, DoubledGaugeField &UUU, 
						   GaugeField & mat,
						   const FermionField &A, const FermionField &B, int dag) {
  assert((dag == DaggerNo) || (dag == DaggerYes));

  Compressor compressor;

  FermionField Btilde(B._grid);
  FermionField Atilde(B._grid);
  Atilde = A;

  st.HaloExchange(B, compressor);

  for (int mu = 0; mu < Nd; mu++) {

    ////////////////////////
    // Call the single hop
    ////////////////////////
    PARALLEL_FOR_LOOP
    for (int sss = 0; sss < B._grid->oSites(); sss++) {
      Kernels::DhopDir(st, U, UUU, st.CommBuf(), sss, sss, B, Btilde, mu,1);
    }

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
void ImprovedStaggeredFermion<Impl>::DhopDeriv(GaugeField &mat, const FermionField &U, const FermionField &V, int dag) {

  conformable(U._grid, _grid);
  conformable(U._grid, V._grid);
  conformable(U._grid, mat._grid);

  mat.checkerboard = U.checkerboard;

  DerivInternal(Stencil, Umu, UUUmu, mat, U, V, dag);
}

template <class Impl>
void ImprovedStaggeredFermion<Impl>::DhopDerivOE(GaugeField &mat, const FermionField &U, const FermionField &V, int dag) {

  conformable(U._grid, _cbgrid);
  conformable(U._grid, V._grid);
  conformable(U._grid, mat._grid);

  assert(V.checkerboard == Even);
  assert(U.checkerboard == Odd);
  mat.checkerboard = Odd;

  DerivInternal(StencilEven, UmuOdd, UUUmuOdd, mat, U, V, dag);
}

template <class Impl>
void ImprovedStaggeredFermion<Impl>::DhopDerivEO(GaugeField &mat, const FermionField &U, const FermionField &V, int dag) {

  conformable(U._grid, _cbgrid);
  conformable(U._grid, V._grid);
  conformable(U._grid, mat._grid);

  assert(V.checkerboard == Odd);
  assert(U.checkerboard == Even);
  mat.checkerboard = Even;

  DerivInternal(StencilOdd, UmuEven, UUUmuEven, mat, U, V, dag);
}

template <class Impl>
void ImprovedStaggeredFermion<Impl>::Dhop(const FermionField &in, FermionField &out, int dag) {
  DhopCalls+=2;
  conformable(in._grid, _grid);  // verifies full grid
  conformable(in._grid, out._grid);

  out.checkerboard = in.checkerboard;

  DhopInternal(Stencil, Lebesgue, Umu, UUUmu, in, out, dag);
}

template <class Impl>
void ImprovedStaggeredFermion<Impl>::DhopOE(const FermionField &in, FermionField &out, int dag) {
  DhopCalls+=1;
  conformable(in._grid, _cbgrid);    // verifies half grid
  conformable(in._grid, out._grid);  // drops the cb check

  assert(in.checkerboard == Even);
  out.checkerboard = Odd;

  DhopInternal(StencilEven, LebesgueEvenOdd, UmuOdd, UUUmuOdd, in, out, dag);
}

template <class Impl>
void ImprovedStaggeredFermion<Impl>::DhopEO(const FermionField &in, FermionField &out, int dag) {
  DhopCalls+=1;
  conformable(in._grid, _cbgrid);    // verifies half grid
  conformable(in._grid, out._grid);  // drops the cb check

  assert(in.checkerboard == Odd);
  out.checkerboard = Even;

  DhopInternal(StencilOdd, LebesgueEvenOdd, UmuEven, UUUmuEven, in, out, dag);
}

template <class Impl>
void ImprovedStaggeredFermion<Impl>::Mdir(const FermionField &in, FermionField &out, int dir, int disp) {
  DhopDir(in, out, dir, disp);
}

template <class Impl>
void ImprovedStaggeredFermion<Impl>::DhopDir(const FermionField &in, FermionField &out, int dir, int disp) {

  Compressor compressor;
  Stencil.HaloExchange(in, compressor);

  PARALLEL_FOR_LOOP
  for (int sss = 0; sss < in._grid->oSites(); sss++) {
    Kernels::DhopDir(Stencil, Umu, UUUmu, Stencil.CommBuf(), sss, sss, in, out, dir, disp);
  }
};

template <class Impl>
void ImprovedStaggeredFermion<Impl>::DhopInternal(StencilImpl &st, LebesgueOrder &lo,
						  DoubledGaugeField &U,
						  DoubledGaugeField &UUU,
						  const FermionField &in,
						  FermionField &out, int dag) 
{
#ifdef GRID_OMP
  if ( StaggeredKernelsStatic::Comms == StaggeredKernelsStatic::CommsAndCompute )
    DhopInternalOverlappedComms(st,lo,U,UUU,in,out,dag);
  else
#endif
    DhopInternalSerialComms(st,lo,U,UUU,in,out,dag);
}
template <class Impl>
void ImprovedStaggeredFermion<Impl>::DhopInternalOverlappedComms(StencilImpl &st, LebesgueOrder &lo,
								 DoubledGaugeField &U,
								 DoubledGaugeField &UUU,
								 const FermionField &in,
								 FermionField &out, int dag) 
{
#ifdef GRID_OMP
  Compressor compressor; 
  int len =  U._grid->oSites();
  const int LLs =  1;

  DhopTotalTime   -= usecond();

  DhopFaceTime    -= usecond();
  st.Prepare();
  st.HaloGather(in,compressor);
  st.CommsMergeSHM(compressor);
  DhopFaceTime    += usecond();

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  // Ugly explicit thread mapping introduced for OPA reasons.
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  DhopComputeTime    -= usecond();
#pragma omp parallel 
  {
    int tid = omp_get_thread_num();
    int nthreads = omp_get_num_threads();
    int ncomms = CartesianCommunicator::nCommThreads;
    if (ncomms == -1) ncomms = 1;
    assert(nthreads > ncomms);

    if (tid >= ncomms) {
      nthreads -= ncomms;
      int ttid  = tid - ncomms;
      int n     = len;
      int chunk = n / nthreads;
      int rem   = n % nthreads;
      int myblock, myn;
      if (ttid < rem) {
        myblock = ttid * chunk + ttid;
        myn = chunk+1;
      } else {
        myblock = ttid*chunk + rem;
        myn = chunk;
      }

      // do the compute
      if (dag == DaggerYes) {
        for (int ss = myblock; ss < myblock+myn; ++ss) {
          int sU = ss;
	  // Interior = 1; Exterior = 0; must implement for staggered
          Kernels::DhopSiteDag(st,lo,U,UUU,st.CommBuf(),1,sU,in,out,1,0); 
        }
      } else {
        for (int ss = myblock; ss < myblock+myn; ++ss) {
	  // Interior = 1; Exterior = 0;
          int sU = ss;
          Kernels::DhopSite(st,lo,U,UUU,st.CommBuf(),1,sU,in,out,1,0);
        }
      }
    } else {
      st.CommunicateThreaded();
    }
  }
  DhopComputeTime    += usecond();

  // First to enter, last to leave timing
  DhopFaceTime    -= usecond();
  st.CommsMerge(compressor);
  DhopFaceTime    -= usecond();

  DhopComputeTime2    -= usecond();
  if (dag == DaggerYes) {
    int sz=st.surface_list.size();
    parallel_for (int ss = 0; ss < sz; ss++) {
      int sU = st.surface_list[ss];
      Kernels::DhopSiteDag(st,lo,U,UUU,st.CommBuf(),1,sU,in,out,0,1);
    }
  } else {
    int sz=st.surface_list.size();
    parallel_for (int ss = 0; ss < sz; ss++) {
      int sU = st.surface_list[ss];
      Kernels::DhopSite(st,lo,U,UUU,st.CommBuf(),1,sU,in,out,0,1);
    }
  }
  DhopComputeTime2    += usecond();
#else
  assert(0);
#endif
}


template <class Impl>
void ImprovedStaggeredFermion<Impl>::DhopInternalSerialComms(StencilImpl &st, LebesgueOrder &lo,
							     DoubledGaugeField &U,
							     DoubledGaugeField &UUU,
							     const FermionField &in,
							     FermionField &out, int dag) 
{
  assert((dag == DaggerNo) || (dag == DaggerYes));

  DhopTotalTime   -= usecond();

  DhopCommTime    -= usecond();
  Compressor compressor;
  st.HaloExchange(in, compressor);
  DhopCommTime    += usecond();

  DhopComputeTime -= usecond();
  if (dag == DaggerYes) {
    parallel_for (int sss = 0; sss < in._grid->oSites(); sss++) {
      Kernels::DhopSiteDag(st, lo, U, UUU, st.CommBuf(), 1, sss, in, out);
    }
  } else {
    parallel_for (int sss = 0; sss < in._grid->oSites(); sss++) {
      Kernels::DhopSite(st, lo, U, UUU, st.CommBuf(), 1, sss, in, out);
    }
  }
  DhopComputeTime += usecond();
  DhopTotalTime   += usecond();
};

  ////////////////////////////////////////////////////////////////
  // Reporting
  ////////////////////////////////////////////////////////////////
template<class Impl>
void ImprovedStaggeredFermion<Impl>::Report(void) 
{
  std::vector<int> latt = GridDefaultLatt();          
  RealD volume = 1;  for(int mu=0;mu<Nd;mu++) volume=volume*latt[mu];
  RealD NP = _grid->_Nprocessors;
  RealD NN = _grid->NodeCount();

  std::cout << GridLogMessage << "#### Dhop calls report " << std::endl;

  std::cout << GridLogMessage << "ImprovedStaggeredFermion Number of DhopEO Calls   : " 
	    << DhopCalls   << std::endl;
  std::cout << GridLogMessage << "ImprovedStaggeredFermion TotalTime   /Calls       : " 
	    << DhopTotalTime   / DhopCalls << " us" << std::endl;
  std::cout << GridLogMessage << "ImprovedStaggeredFermion CommTime    /Calls       : " 
	    << DhopCommTime    / DhopCalls << " us" << std::endl;
  std::cout << GridLogMessage << "ImprovedStaggeredFermion ComputeTime/Calls        : " 
	    << DhopComputeTime / DhopCalls << " us" << std::endl;

  // Average the compute time
  _grid->GlobalSum(DhopComputeTime);
  DhopComputeTime/=NP;

  RealD mflops = 1154*volume*DhopCalls/DhopComputeTime/2; // 2 for red black counting
  std::cout << GridLogMessage << "Average mflops/s per call                : " << mflops << std::endl;
  std::cout << GridLogMessage << "Average mflops/s per call per rank       : " << mflops/NP << std::endl;
  std::cout << GridLogMessage << "Average mflops/s per call per node       : " << mflops/NN << std::endl;
  
  RealD Fullmflops = 1154*volume*DhopCalls/(DhopTotalTime)/2; // 2 for red black counting
  std::cout << GridLogMessage << "Average mflops/s per call (full)         : " << Fullmflops << std::endl;
  std::cout << GridLogMessage << "Average mflops/s per call per rank (full): " << Fullmflops/NP << std::endl;
  std::cout << GridLogMessage << "Average mflops/s per call per node (full): " << Fullmflops/NN << std::endl;

  std::cout << GridLogMessage << "ImprovedStaggeredFermion Stencil"    <<std::endl;  Stencil.Report();
  std::cout << GridLogMessage << "ImprovedStaggeredFermion StencilEven"<<std::endl;  StencilEven.Report();
  std::cout << GridLogMessage << "ImprovedStaggeredFermion StencilOdd" <<std::endl;  StencilOdd.Report();
}
template<class Impl>
void ImprovedStaggeredFermion<Impl>::ZeroCounters(void) 
{
  DhopCalls       = 0;
  DhopTotalTime   = 0;
  DhopCommTime    = 0;
  DhopComputeTime = 0;
  DhopFaceTime    = 0;

  Stencil.ZeroCounters();
  StencilEven.ZeroCounters();
  StencilOdd.ZeroCounters();
}


//////////////////////////////////////////////////////// 
// Conserved current - not yet implemented.
////////////////////////////////////////////////////////
template <class Impl>
void ImprovedStaggeredFermion<Impl>::ContractConservedCurrent(PropagatorField &q_in_1,
                                                        PropagatorField &q_in_2,
                                                        PropagatorField &q_out,
                                                        Current curr_type,
                                                        unsigned int mu)
{
    assert(0);
}

template <class Impl>
void ImprovedStaggeredFermion<Impl>::SeqConservedCurrent(PropagatorField &q_in, 
                                              PropagatorField &q_out,
                                              Current curr_type,
                                              unsigned int mu,
                                              unsigned int tmin, 
                                              unsigned int tmax,
					      ComplexField &lattice_cmplx)
{
    assert(0);

}


FermOpStaggeredTemplateInstantiate(ImprovedStaggeredFermion);

  //AdjointFermOpTemplateInstantiate(ImprovedStaggeredFermion);
  //TwoIndexFermOpTemplateInstantiate(ImprovedStaggeredFermion);

}}
