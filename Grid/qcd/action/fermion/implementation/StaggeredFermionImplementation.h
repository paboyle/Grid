/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/action/fermion/StaggeredFermion.cc

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
StaggeredFermion<Impl>::StaggeredFermion(GridCartesian &Fgrid, GridRedBlackCartesian &Hgrid,
                                         RealD _mass,
                                         RealD _c1,
                                         RealD _u0,
                                         const ImplParams &p)
  : Kernels(p),
    _grid(&Fgrid),
    _cbgrid(&Hgrid),
    Stencil(&Fgrid, npoint, Even, directions, displacements,p),
    StencilEven(&Hgrid, npoint, Even, directions, displacements,p),  // source is Even
    StencilOdd(&Hgrid, npoint, Odd, directions, displacements,p),  // source is Odd
    mass(_mass),
    Lebesgue(_grid),
    LebesgueEvenOdd(_cbgrid),
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
StaggeredFermion<Impl>::StaggeredFermion(GaugeField &_Uthin, GridCartesian &Fgrid,
							 GridRedBlackCartesian &Hgrid, RealD _mass,
							 RealD _c1, RealD _u0,
							 const ImplParams &p)
  : StaggeredFermion(Fgrid,Hgrid,_mass,_c1,_u0,p)
{
  ImportGauge(_Uthin);
}

////////////////////////////////////////////////////////////
// Momentum space propagator should be 
// https://arxiv.org/pdf/hep-lat/9712010.pdf
//
// mom space action.
//   gamma_mu i ( c1 sin pmu) + m
//
// must track through staggered flavour/spin reduction in literature to 
// turn to free propagator for the one component chi field, a la page 4/5
// of above link to implmement fourier based solver.
////////////////////////////////////////////////////////////
template <class Impl>
void StaggeredFermion<Impl>::ImportGaugeSimple(const GaugeField &_Ufat)
{
  /////////////////////////////////////////////////////////////////
  // Trivial import; phases and fattening and such like preapplied
  /////////////////////////////////////////////////////////////////
  GaugeLinkField U(GaugeGrid());

  for (int mu = 0; mu < Nd; mu++) {

    U = PeekIndex<LorentzIndex>(_Ufat, mu);
    PokeIndex<LorentzIndex>(Umu, U, mu);

    U = adj( Cshift(U, mu, -1));
    PokeIndex<LorentzIndex>(Umu, -U, mu+4);

  }
  CopyGaugeCheckerboards();
}
template <class Impl>
void StaggeredFermion<Impl>::ImportGaugeSimple(const DoubledGaugeField &_U)
{

  Umu   = _U;
  CopyGaugeCheckerboards();
}

template <class Impl>
void StaggeredFermion<Impl>::CopyGaugeCheckerboards(void)
{
  pickCheckerboard(Even, UmuEven,  Umu);
  pickCheckerboard(Odd,  UmuOdd ,  Umu);
}
template <class Impl>
void StaggeredFermion<Impl>::ImportGauge(const GaugeField &_Uthin)
{
  GaugeLinkField U(GaugeGrid());

  ////////////////////////////////////////////////////////
  // Double Store should take two fields for Naik and one hop separately.
  ////////////////////////////////////////////////////////
  Impl::DoubleStore(GaugeGrid(), Umu, _Uthin);

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
RealD StaggeredFermion<Impl>::M(const FermionField &in, FermionField &out) {
  out.Checkerboard() = in.Checkerboard();
  Dhop(in, out, DaggerNo);
  return axpy_norm(out, mass, in, out);
}

template <class Impl>
RealD StaggeredFermion<Impl>::Mdag(const FermionField &in, FermionField &out) {
  out.Checkerboard() = in.Checkerboard();
  Dhop(in, out, DaggerYes);
  return axpy_norm(out, mass, in, out);
}

template <class Impl>
void StaggeredFermion<Impl>::Meooe(const FermionField &in, FermionField &out) {
  if (in.Checkerboard() == Odd) {
    DhopEO(in, out, DaggerNo);
  } else {
    DhopOE(in, out, DaggerNo);
  }
}
template <class Impl>
void StaggeredFermion<Impl>::MeooeDag(const FermionField &in, FermionField &out) {
  if (in.Checkerboard() == Odd) {
    DhopEO(in, out, DaggerYes);
  } else {
    DhopOE(in, out, DaggerYes);
  }
}

template <class Impl>
void StaggeredFermion<Impl>::Mooee(const FermionField &in, FermionField &out) {
  out.Checkerboard() = in.Checkerboard();
  typename FermionField::scalar_type scal(mass);
  out = scal * in;
}

template <class Impl>
void StaggeredFermion<Impl>::MooeeDag(const FermionField &in, FermionField &out) {
  out.Checkerboard() = in.Checkerboard();
  Mooee(in, out);
}

template <class Impl>
void StaggeredFermion<Impl>::MooeeInv(const FermionField &in, FermionField &out) {
  out.Checkerboard() = in.Checkerboard();
  out = (1.0 / (mass)) * in;
}

template <class Impl>
void StaggeredFermion<Impl>::MooeeInvDag(const FermionField &in,
						 FermionField &out) {
  out.Checkerboard() = in.Checkerboard();
  MooeeInv(in, out);
}

///////////////////////////////////
// Internal
///////////////////////////////////

template <class Impl>
void StaggeredFermion<Impl>::DerivInternal(StencilImpl &st,
                                           DoubledGaugeField &U,
                                           GaugeField & mat,
                                           const FermionField &A,
                                           const FermionField &B, int dag) {
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
    auto U_v   = U.View();
    auto B_v   = B.View();
    auto Btilde_v   = Btilde.View();
    thread_for(sss,B.Grid()->oSites(),{
      Kernels::DhopDirKernel(st, U_v, st.CommBuf(), sss, sss, B_v, Btilde_v, mu,1);
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
      
    assert(0);// need to figure out the force interface with a blasted three link term.
    
  }
}

template <class Impl>
void StaggeredFermion<Impl>::DhopDeriv(GaugeField &mat, const FermionField &U, const FermionField &V, int dag) {

  conformable(U.Grid(), _grid);
  conformable(U.Grid(), V.Grid());
  conformable(U.Grid(), mat.Grid());

  mat.Checkerboard() = U.Checkerboard();

  DerivInternal(Stencil, Umu, mat, U, V, dag);
}

template <class Impl>
void StaggeredFermion<Impl>::DhopDerivOE(GaugeField &mat, const FermionField &U, const FermionField &V, int dag) {

  conformable(U.Grid(), _cbgrid);
  conformable(U.Grid(), V.Grid());
  conformable(U.Grid(), mat.Grid());

  assert(V.Checkerboard() == Even);
  assert(U.Checkerboard() == Odd);
  mat.Checkerboard() = Odd;

  DerivInternal(StencilEven, UmuOdd, mat, U, V, dag);
}

template <class Impl>
void StaggeredFermion<Impl>::DhopDerivEO(GaugeField &mat, const FermionField &U, const FermionField &V, int dag) {

  conformable(U.Grid(), _cbgrid);
  conformable(U.Grid(), V.Grid());
  conformable(U.Grid(), mat.Grid());

  assert(V.Checkerboard() == Odd);
  assert(U.Checkerboard() == Even);
  mat.Checkerboard() = Even;

  DerivInternal(StencilOdd, UmuEven, mat, U, V, dag);
}

template <class Impl>
void StaggeredFermion<Impl>::Dhop(const FermionField &in, FermionField &out, int dag)
{
  DhopCalls+=2;
  conformable(in.Grid(), _grid);  // verifies full grid
  conformable(in.Grid(), out.Grid());

  out.Checkerboard() = in.Checkerboard();

  DhopInternal(Stencil, Lebesgue, Umu, in, out, dag);
}

template <class Impl>
void StaggeredFermion<Impl>::DhopOE(const FermionField &in, FermionField &out, int dag)
{
  DhopCalls+=1;
  conformable(in.Grid(), _cbgrid);    // verifies half grid
  conformable(in.Grid(), out.Grid());  // drops the cb check

  assert(in.Checkerboard() == Even);
  out.Checkerboard() = Odd;

  DhopInternal(StencilEven, LebesgueEvenOdd, UmuOdd, in, out, dag);
}

template <class Impl>
void StaggeredFermion<Impl>::DhopEO(const FermionField &in, FermionField &out, int dag)
{
  DhopCalls+=1;
  conformable(in.Grid(), _cbgrid);    // verifies half grid
  conformable(in.Grid(), out.Grid());  // drops the cb check

  assert(in.Checkerboard() == Odd);
  out.Checkerboard() = Even;

  DhopInternal(StencilOdd, LebesgueEvenOdd, UmuEven, in, out, dag);
}

template <class Impl>
void StaggeredFermion<Impl>::Mdir(const FermionField &in, FermionField &out, int dir, int disp) {
  DhopDir(in, out, dir, disp);
}

template <class Impl>
void StaggeredFermion<Impl>::DhopDir(const FermionField &in, FermionField &out, int dir, int disp) {

  Compressor compressor;
  Stencil.HaloExchange(in, compressor);
  auto Umu_v   =   Umu.View();
  auto in_v    =  in.View();
  auto out_v   = out.View();
  thread_for( sss, in.Grid()->oSites(),{
    Kernels::DhopDirKernel(Stencil, Umu_v, Stencil.CommBuf(), sss, sss, in_v, out_v, dir, disp);
  });
};

template <class Impl>
void StaggeredFermion<Impl>::DhopInternal(StencilImpl &st, LebesgueOrder &lo,
						  DoubledGaugeField &U,
						  const FermionField &in,
						  FermionField &out, int dag) 
{
#ifdef GRID_OMP
  if ( StaggeredKernelsStatic::Comms == StaggeredKernelsStatic::CommsAndCompute )
    DhopInternalOverlappedComms(st,lo,U,in,out,dag);
  else
#endif
    DhopInternalSerialComms(st,lo,U,in,out,dag);
}
template <class Impl>
void StaggeredFermion<Impl>::DhopInternalOverlappedComms(StencilImpl &st, LebesgueOrder &lo,
								 DoubledGaugeField &U,
								 const FermionField &in,
								 FermionField &out, int dag) 
{
#ifdef GRID_OMP
  Compressor compressor; 
  int len =  U.Grid()->oSites();
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
      auto U_v   = U.View();
      auto in_v  = in.View();
      auto out_v = out.View();
      if (dag == DaggerYes) {
        for (int ss = myblock; ss < myblock+myn; ++ss) {
          int sU = ss;
	  // Interior = 1; Exterior = 0; must implement for staggered
          Kernels::DhopSiteDag(st,lo,U_v,st.CommBuf(),1,sU,in_v,out_v,1,0);
        }
      } else {
        for (int ss = myblock; ss < myblock+myn; ++ss) {
	  // Interior = 1; Exterior = 0;
          int sU = ss;
          Kernels::DhopSite(st,lo,U_v,st.CommBuf(),1,sU,in_v,out_v,1,0);
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
  {
    auto U_v   = U.View();
    auto in_v  = in.View();
    auto out_v = out.View();
    if (dag == DaggerYes) {
      int sz=st.surface_list.size();
      thread_for(ss,sz,{
	int sU = st.surface_list[ss];
	Kernels::DhopSiteDag(st,lo,U_v,st.CommBuf(),1,sU,in_v,out_v,0,1);
      });
    } else {
      int sz=st.surface_list.size();
      thread_for(ss,sz,{
	int sU = st.surface_list[ss];
	Kernels::DhopSite(st,lo,U_v,st.CommBuf(),1,sU,in_v,out_v,0,1);
      });
    }
  }
  DhopComputeTime2    += usecond();
#else
  assert(0);
#endif
}


template <class Impl>
void StaggeredFermion<Impl>::DhopInternalSerialComms(StencilImpl &st, LebesgueOrder &lo,
							     DoubledGaugeField &U,
							     const FermionField &in,
							     FermionField &out, int dag) 
{
  assert((dag == DaggerNo) || (dag == DaggerYes));

  DhopTotalTime   -= usecond();

  DhopCommTime    -= usecond();
  Compressor compressor;
  st.HaloExchange(in, compressor);
  DhopCommTime    += usecond();

  auto U_v   =   U.View();
  auto in_v  =  in.View();
  auto out_v = out.View();
  DhopComputeTime -= usecond();
  if (dag == DaggerYes) {
    thread_for(sss, in.Grid()->oSites(),{
      Kernels::DhopSiteDag(st, lo, U_v, st.CommBuf(), 1, sss, in_v, out_v);
    });
  } else {
    thread_for(sss, in.Grid()->oSites(),{
      Kernels::DhopSite(st, lo, U_v, st.CommBuf(), 1, sss, in_v, out_v);
    });
  }
  DhopComputeTime += usecond();
  DhopTotalTime   += usecond();
};

  ////////////////////////////////////////////////////////////////
  // Reporting
  ////////////////////////////////////////////////////////////////
template<class Impl>
void StaggeredFermion<Impl>::Report(void)
{
  Coordinate latt = _grid->GlobalDimensions();
  RealD volume = 1;  for(int mu=0;mu<Nd;mu++) volume=volume*latt[mu];
  RealD NP = _grid->_Nprocessors;
  RealD NN = _grid->NodeCount();

  std::cout << GridLogMessage << "#### Dhop calls report " << std::endl;

  std::cout << GridLogMessage << "StaggeredFermion Number of DhopEO Calls   : "
	    << DhopCalls   << std::endl;
  std::cout << GridLogMessage << "StaggeredFermion TotalTime   /Calls       : "
	    << DhopTotalTime   / DhopCalls << " us" << std::endl;
  std::cout << GridLogMessage << "StaggeredFermion CommTime    /Calls       : "
	    << DhopCommTime    / DhopCalls << " us" << std::endl;
  std::cout << GridLogMessage << "StaggeredFermion ComputeTime/Calls        : "
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

  std::cout << GridLogMessage << "StaggeredFermion Stencil"    <<std::endl;  Stencil.Report();
  std::cout << GridLogMessage << "StaggeredFermion StencilEven"<<std::endl;  StencilEven.Report();
  std::cout << GridLogMessage << "StaggeredFermion StencilOdd" <<std::endl;  StencilOdd.Report();
}
template<class Impl>
void StaggeredFermion<Impl>::ZeroCounters(void)
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
void StaggeredFermion<Impl>::ContractConservedCurrent(PropagatorField &q_in_1,
							      PropagatorField &q_in_2,
							      PropagatorField &q_out,
							      Current curr_type,
							      unsigned int mu)
{
  assert(0);
}

template <class Impl>
void StaggeredFermion<Impl>::SeqConservedCurrent(PropagatorField &q_in,
                                                         PropagatorField &q_out,
                                                         Current curr_type,
                                                         unsigned int mu, 
                                                         unsigned int tmin,
                                              unsigned int tmax,
					      ComplexField &lattice_cmplx)
{
  assert(0);

}

NAMESPACE_END(Grid);
