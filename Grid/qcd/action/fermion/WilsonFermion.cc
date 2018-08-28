
/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/action/fermion/WilsonFermion.cc

Copyright (C) 2015

Author: Peter Boyle <pabobyle@ph.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: Peter Boyle <peterboyle@Peters-MacBook-Pro-2.local>
Author: paboyle <paboyle@ph.ed.ac.uk>

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
#include <Grid/qcd/action/fermion/FermionCore.h>
#include <Grid/qcd/action/fermion/WilsonFermion.h>

namespace Grid {
namespace QCD {

const std::vector<int> WilsonFermionStatic::directions({0, 1, 2, 3, 0, 1, 2, 3});
const std::vector<int> WilsonFermionStatic::displacements({1, 1, 1, 1, -1, -1, -1, -1});
int WilsonFermionStatic::HandOptDslash;

/////////////////////////////////
// Constructor and gauge import
/////////////////////////////////

template <class Impl>
WilsonFermion<Impl>::WilsonFermion(GaugeField &_Umu, GridCartesian &Fgrid,
                                   GridRedBlackCartesian &Hgrid, RealD _mass,
                                   const ImplParams &p,
                                   const WilsonAnisotropyCoefficients &anis)
    : Kernels(p),
      _grid(&Fgrid),
      _cbgrid(&Hgrid),
      Stencil(&Fgrid, npoint, Even, directions, displacements),
      StencilEven(&Hgrid, npoint, Even, directions,displacements),  // source is Even
      StencilOdd(&Hgrid, npoint, Odd, directions,displacements),  // source is Odd
      mass(_mass),
      Lebesgue(_grid),
      LebesgueEvenOdd(_cbgrid),
      Umu(&Fgrid),
      UmuEven(&Hgrid),
      UmuOdd(&Hgrid),
      _tmp(&Hgrid),
      anisotropyCoeff(anis)
{
  // Allocate the required comms buffer
  ImportGauge(_Umu);
  if  (anisotropyCoeff.isAnisotropic){
    diag_mass = mass + 1.0 + (Nd-1)*(anisotropyCoeff.nu / anisotropyCoeff.xi_0);
  } else {
    diag_mass = 4.0 + mass;
  }


}

template <class Impl>
void WilsonFermion<Impl>::ImportGauge(const GaugeField &_Umu) {
  GaugeField HUmu(_Umu._grid);

  //Here multiply the anisotropy coefficients
  if (anisotropyCoeff.isAnisotropic)
  {

    for (int mu = 0; mu < Nd; mu++)
    {
      GaugeLinkField U_dir = (-0.5)*PeekIndex<LorentzIndex>(_Umu, mu);
      if (mu != anisotropyCoeff.t_direction)
        U_dir *= (anisotropyCoeff.nu / anisotropyCoeff.xi_0);

      PokeIndex<LorentzIndex>(HUmu, U_dir, mu);
    }
  }
  else
  {
    HUmu = _Umu * (-0.5);
  }
  Impl::DoubleStore(GaugeGrid(), Umu, HUmu);
  pickCheckerboard(Even, UmuEven, Umu);
  pickCheckerboard(Odd, UmuOdd, Umu);
}

/////////////////////////////
// Implement the interface
/////////////////////////////

template <class Impl>
RealD WilsonFermion<Impl>::M(const FermionField &in, FermionField &out) {
  out.checkerboard = in.checkerboard;
  Dhop(in, out, DaggerNo);
  return axpy_norm(out, diag_mass, in, out);
}

template <class Impl>
RealD WilsonFermion<Impl>::Mdag(const FermionField &in, FermionField &out) {
  out.checkerboard = in.checkerboard;
  Dhop(in, out, DaggerYes);
  return axpy_norm(out, diag_mass, in, out);
}

template <class Impl>
void WilsonFermion<Impl>::Meooe(const FermionField &in, FermionField &out) {
  if (in.checkerboard == Odd) {
    DhopEO(in, out, DaggerNo);
  } else {
    DhopOE(in, out, DaggerNo);
  }
}

template <class Impl>
void WilsonFermion<Impl>::MeooeDag(const FermionField &in, FermionField &out) {
  if (in.checkerboard == Odd) {
    DhopEO(in, out, DaggerYes);
  } else {
    DhopOE(in, out, DaggerYes);
  }
}
  
template <class Impl>
void WilsonFermion<Impl>::Mooee(const FermionField &in, FermionField &out) {
  out.checkerboard = in.checkerboard;
  typename FermionField::scalar_type scal(diag_mass);
  out = scal * in;
}

template <class Impl>
void WilsonFermion<Impl>::MooeeDag(const FermionField &in, FermionField &out) {
  out.checkerboard = in.checkerboard;
  Mooee(in, out);
}

template<class Impl>
void WilsonFermion<Impl>::MooeeInv(const FermionField &in, FermionField &out) {
  out.checkerboard = in.checkerboard;
  out = (1.0/(diag_mass))*in;
}
  
template<class Impl>
void WilsonFermion<Impl>::MooeeInvDag(const FermionField &in, FermionField &out) {
  out.checkerboard = in.checkerboard;
  MooeeInv(in,out);
}
template<class Impl>
void WilsonFermion<Impl>::MomentumSpacePropagator(FermionField &out, const FermionField &in,RealD _m,std::vector<double> twist)
{  
  typedef typename FermionField::vector_type vector_type;
  typedef typename FermionField::scalar_type ScalComplex;
  typedef Lattice<iSinglet<vector_type> > LatComplex;
  
  // what type LatticeComplex 
  conformable(_grid,out._grid);
  
  Gamma::Algebra Gmu [] = {
    Gamma::Algebra::GammaX,
    Gamma::Algebra::GammaY,
    Gamma::Algebra::GammaZ,
    Gamma::Algebra::GammaT
  };
  
  std::vector<int> latt_size   = _grid->_fdimensions;
  
  FermionField   num  (_grid); num  = zero;
  LatComplex    wilson(_grid); wilson= zero;
  LatComplex     one  (_grid); one = ScalComplex(1.0,0.0);
  
  LatComplex denom(_grid); denom= zero;
  LatComplex kmu(_grid); 
  ScalComplex ci(0.0,1.0);
  // momphase = n * 2pi / L
  for(int mu=0;mu<Nd;mu++) {
    
    LatticeCoordinate(kmu,mu);
    
    RealD TwoPiL =  M_PI * 2.0/ latt_size[mu];
    
    kmu = TwoPiL * kmu;
    kmu = kmu + TwoPiL * one * twist[mu];//momentum for twisted boundary conditions
    
    wilson = wilson + 2.0*sin(kmu*0.5)*sin(kmu*0.5); // Wilson term
    
    num = num - sin(kmu)*ci*(Gamma(Gmu[mu])*in);    // derivative term
    
    denom=denom + sin(kmu)*sin(kmu);
  }
  
  wilson = wilson + _m;     // 2 sin^2 k/2 + m
  
  num   = num + wilson*in;     // -i gmu sin k + 2 sin^2 k/2 + m
  
  denom= denom+wilson*wilson; // sin^2 k + (2 sin^2 k/2 + m)^2
  
  denom= one/denom;
  
  out = num*denom; // [ -i gmu sin k + 2 sin^2 k/2 + m] / [ sin^2 k + (2 sin^2 k/2 + m)^2 ]
  
}
  

///////////////////////////////////
// Internal
///////////////////////////////////

template <class Impl>
void WilsonFermion<Impl>::DerivInternal(StencilImpl &st, DoubledGaugeField &U,
                                        GaugeField &mat, const FermionField &A,
                                        const FermionField &B, int dag) {
  assert((dag == DaggerNo) || (dag == DaggerYes));

  Compressor compressor(dag);

  FermionField Btilde(B._grid);
  FermionField Atilde(B._grid);
  Atilde = A;//redundant

  st.HaloExchange(B, compressor);

  for (int mu = 0; mu < Nd; mu++) {
    ////////////////////////////////////////////////////////////////////////
    // Flip gamma (1+g)<->(1-g) if dag
    ////////////////////////////////////////////////////////////////////////
    int gamma = mu;
    if (!dag) gamma += Nd;

    ////////////////////////
    // Call the single hop
    ////////////////////////
    parallel_for (int sss = 0; sss < B._grid->oSites(); sss++) {
      Kernels::DhopDir(st, U, st.CommBuf(), sss, sss, B, Btilde, mu, gamma);
    }

    //////////////////////////////////////////////////
    // spin trace outer product
    //////////////////////////////////////////////////
    Impl::InsertForce4D(mat, Btilde, Atilde, mu);
  }
}

template <class Impl>
void WilsonFermion<Impl>::DhopDeriv(GaugeField &mat, const FermionField &U, const FermionField &V, int dag) {
  conformable(U._grid, _grid);
  conformable(U._grid, V._grid);
  conformable(U._grid, mat._grid);

  mat.checkerboard = U.checkerboard;

  DerivInternal(Stencil, Umu, mat, U, V, dag);
}

template <class Impl>
void WilsonFermion<Impl>::DhopDerivOE(GaugeField &mat, const FermionField &U, const FermionField &V, int dag) {
  conformable(U._grid, _cbgrid);
  conformable(U._grid, V._grid);
  //conformable(U._grid, mat._grid); not general, leaving as a comment (Guido)
  // Motivation: look at the SchurDiff operator
  
  assert(V.checkerboard == Even);
  assert(U.checkerboard == Odd);
  mat.checkerboard = Odd;

  DerivInternal(StencilEven, UmuOdd, mat, U, V, dag);
}

template <class Impl>
void WilsonFermion<Impl>::DhopDerivEO(GaugeField &mat, const FermionField &U, const FermionField &V, int dag) {
  conformable(U._grid, _cbgrid);
  conformable(U._grid, V._grid);
  //conformable(U._grid, mat._grid);

  assert(V.checkerboard == Odd);
  assert(U.checkerboard == Even);
  mat.checkerboard = Even;

  DerivInternal(StencilOdd, UmuEven, mat, U, V, dag);
}

template <class Impl>
void WilsonFermion<Impl>::Dhop(const FermionField &in, FermionField &out, int dag) {
  conformable(in._grid, _grid);  // verifies full grid
  conformable(in._grid, out._grid);

  out.checkerboard = in.checkerboard;

  DhopInternal(Stencil, Lebesgue, Umu, in, out, dag);
}

template <class Impl>
void WilsonFermion<Impl>::DhopOE(const FermionField &in, FermionField &out, int dag) {
  conformable(in._grid, _cbgrid);    // verifies half grid
  conformable(in._grid, out._grid);  // drops the cb check

  assert(in.checkerboard == Even);
  out.checkerboard = Odd;

  DhopInternal(StencilEven, LebesgueEvenOdd, UmuOdd, in, out, dag);
}

template <class Impl>
void WilsonFermion<Impl>::DhopEO(const FermionField &in, FermionField &out,int dag) {
  conformable(in._grid, _cbgrid);    // verifies half grid
  conformable(in._grid, out._grid);  // drops the cb check

  assert(in.checkerboard == Odd);
  out.checkerboard = Even;

  DhopInternal(StencilOdd, LebesgueEvenOdd, UmuEven, in, out, dag);
}

template <class Impl>
void WilsonFermion<Impl>::Mdir(const FermionField &in, FermionField &out, int dir, int disp) {
  DhopDir(in, out, dir, disp);
}

template <class Impl>
void WilsonFermion<Impl>::DhopDir(const FermionField &in, FermionField &out, int dir, int disp) {
  int skip = (disp == 1) ? 0 : 1;
  int dirdisp = dir + skip * 4;
  int gamma = dir + (1 - skip) * 4;

  DhopDirDisp(in, out, dirdisp, gamma, DaggerNo);
};

template <class Impl>
void WilsonFermion<Impl>::DhopDirDisp(const FermionField &in, FermionField &out,int dirdisp, int gamma, int dag) {
  Compressor compressor(dag);

  Stencil.HaloExchange(in, compressor);

  parallel_for (int sss = 0; sss < in._grid->oSites(); sss++) {
    Kernels::DhopDir(Stencil, Umu, Stencil.CommBuf(), sss, sss, in, out, dirdisp, gamma);
  }
} 
/*Change starts*/
template <class Impl>
void WilsonFermion<Impl>::DhopInternal(StencilImpl &st, LebesgueOrder &lo,
                                       DoubledGaugeField &U,
                                       const FermionField &in,
                                       FermionField &out, int dag) {
#ifdef GRID_OMP
  if ( WilsonKernelsStatic::Comms == WilsonKernelsStatic::CommsAndCompute )
    DhopInternalOverlappedComms(st,lo,U,in,out,dag);
  else
#endif 
    DhopInternalSerial(st,lo,U,in,out,dag);

}

template <class Impl>
void WilsonFermion<Impl>::DhopInternalOverlappedComms(StencilImpl &st, LebesgueOrder &lo,
                                       DoubledGaugeField &U,
                                       const FermionField &in,
                                       FermionField &out, int dag) {
  assert((dag == DaggerNo) || (dag == DaggerYes));
#ifdef GRID_OMP
  Compressor compressor;
  int len =  U._grid->oSites();
  const int LLs =  1;

  st.Prepare();
  st.HaloGather(in,compressor);
  st.CommsMergeSHM(compressor);
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

        for (int sss = myblock; sss < myblock+myn; ++sss) {
         Kernels::DhopSiteDag(st, lo, U, st.CommBuf(), sss, sss, 1, 1, in, out);
       }
     } else {
        for (int sss = myblock; sss < myblock+myn; ++sss) {
         Kernels::DhopSite(st, lo, U, st.CommBuf(), sss, sss, 1, 1, in, out);
       }
    } //else

    } else {
      st.CommunicateThreaded();
    }

  Compressor compressor(dag);

  if (dag == DaggerYes) {
    parallel_for (int sss = 0; sss < in._grid->oSites(); sss++) {
      Kernels::DhopSiteDag(st, lo, U, st.CommBuf(), sss, sss, 1, 1, in, out);
    }
  } else {
    parallel_for (int sss = 0; sss < in._grid->oSites(); sss++) {
      Kernels::DhopSite(st, lo, U, st.CommBuf(), sss, sss, 1, 1, in, out);
    }
  }

  }  //pragma
#else
  assert(0);
#endif
};


template <class Impl>
void WilsonFermion<Impl>::DhopInternalSerial(StencilImpl &st, LebesgueOrder &lo,
                                       DoubledGaugeField &U,
                                       const FermionField &in,
                                       FermionField &out, int dag) {
  assert((dag == DaggerNo) || (dag == DaggerYes));
  Compressor compressor(dag);
  st.HaloExchange(in, compressor);

  if (dag == DaggerYes) {
    parallel_for (int sss = 0; sss < in._grid->oSites(); sss++) {
      Kernels::DhopSiteDag(st, lo, U, st.CommBuf(), sss, sss, 1, 1, in, out);
    }
  } else {
    parallel_for (int sss = 0; sss < in._grid->oSites(); sss++) {
      Kernels::DhopSite(st, lo, U, st.CommBuf(), sss, sss, 1, 1, in, out);
    }
  }
};
/*Change ends */

/*******************************************************************************
 * Conserved current utilities for Wilson fermions, for contracting propagators
 * to make a conserved current sink or inserting the conserved current 
 * sequentially.
 ******************************************************************************/
template <class Impl>
void WilsonFermion<Impl>::ContractConservedCurrent(PropagatorField &q_in_1,
                                                   PropagatorField &q_in_2,
                                                   PropagatorField &q_out,
                                                   Current curr_type,
                                                   unsigned int mu)
{
    Gamma g5(Gamma::Algebra::Gamma5);
    conformable(_grid, q_in_1._grid);
    conformable(_grid, q_in_2._grid);
    conformable(_grid, q_out._grid);
    PropagatorField tmp1(_grid), tmp2(_grid);
    q_out = zero;

    // Forward, need q1(x + mu), q2(x). Backward, need q1(x), q2(x + mu).
    // Inefficient comms method but not performance critical.
    tmp1 = Cshift(q_in_1, mu, 1);
    tmp2 = Cshift(q_in_2, mu, 1);
    parallel_for (unsigned int sU = 0; sU < Umu._grid->oSites(); ++sU)
    {
        Kernels::ContractConservedCurrentSiteFwd(tmp1._odata[sU],
                                                 q_in_2._odata[sU],
                                                 q_out._odata[sU],
                                                 Umu, sU, mu);
        Kernels::ContractConservedCurrentSiteBwd(q_in_1._odata[sU],
                                                 tmp2._odata[sU],
                                                 q_out._odata[sU],
                                                 Umu, sU, mu);
    }
}


template <class Impl>
void WilsonFermion<Impl>::SeqConservedCurrent(PropagatorField &q_in, 
                                              PropagatorField &q_out,
                                              Current curr_type,
                                              unsigned int mu,
                                              unsigned int tmin, 
                                              unsigned int tmax,
					      ComplexField &lattice_cmplx)
{
    conformable(_grid, q_in._grid);
    conformable(_grid, q_out._grid);
    PropagatorField tmpFwd(_grid), tmpBwd(_grid), tmp(_grid);
    unsigned int tshift = (mu == Tp) ? 1 : 0;
    unsigned int LLt    = GridDefaultLatt()[Tp];

    q_out = zero;
    LatticeInteger coords(_grid);
    LatticeCoordinate(coords, Tp);

    // Need q(x + mu) and q(x - mu).
    tmp = Cshift(q_in, mu, 1);
    tmpFwd = tmp*lattice_cmplx;
    tmp = lattice_cmplx*q_in;
    tmpBwd = Cshift(tmp, mu, -1);

    parallel_for (unsigned int sU = 0; sU < Umu._grid->oSites(); ++sU)
    {
        // Compute the sequential conserved current insertion only if our simd
        // object contains a timeslice we need.
        vInteger t_mask   = ((coords._odata[sU] >= tmin) &&
                             (coords._odata[sU] <= tmax));
        Integer timeSlices = Reduce(t_mask);

        if (timeSlices > 0)
        {
            Kernels::SeqConservedCurrentSiteFwd(tmpFwd._odata[sU], 
                                                q_out._odata[sU], 
                                                Umu, sU, mu, t_mask);
        }

        // Repeat for backward direction.
        t_mask     = ((coords._odata[sU] >= (tmin + tshift)) && 
                      (coords._odata[sU] <= (tmax + tshift)));

	//if tmax = LLt-1 (last timeslice) include timeslice 0 if the time is shifted (mu=3)	
	unsigned int t0 = 0;
	if((tmax==LLt-1) && (tshift==1)) t_mask = (t_mask || (coords._odata[sU] == t0 ));

        timeSlices = Reduce(t_mask);

        if (timeSlices > 0)
        {
            Kernels::SeqConservedCurrentSiteBwd(tmpBwd._odata[sU], 
                                                q_out._odata[sU], 
                                                Umu, sU, mu, t_mask);
        }
    }


}

FermOpTemplateInstantiate(WilsonFermion);
AdjointFermOpTemplateInstantiate(WilsonFermion);
TwoIndexFermOpTemplateInstantiate(WilsonFermion);
GparityFermOpTemplateInstantiate(WilsonFermion);
}
}
