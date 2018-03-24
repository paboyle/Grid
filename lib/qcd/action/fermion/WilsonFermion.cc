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

NAMESPACE_BEGIN(Grid);

const std::vector<int> WilsonFermionStatic::directions({0, 1, 2, 3, 0, 1, 2, 3});
const std::vector<int> WilsonFermionStatic::displacements({1, 1, 1, 1, -1, -1, -1, -1});
int WilsonFermionStatic::HandOptDslash;

/////////////////////////////////
// Constructor and gauge import
/////////////////////////////////

template <class Impl>
WilsonFermion<Impl>::WilsonFermion(GaugeField &_Umu, GridCartesian &Fgrid,
                                   GridRedBlackCartesian &Hgrid, RealD _mass,
                                   const ImplParams &p)
  : Impl(p),
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
    _tmp(&Hgrid)
{
  // Allocate the required comms buffer
  ImportGauge(_Umu);
}

template <class Impl>
void WilsonFermion<Impl>::ImportGauge(const GaugeField &_Umu) {
  GaugeField HUmu(_Umu.Grid());
  HUmu = _Umu * (-0.5);
  Impl::DoubleStore(GaugeGrid(), Umu, HUmu);
  pickCheckerboard(Even, UmuEven, Umu);
  pickCheckerboard(Odd, UmuOdd, Umu);
}

/////////////////////////////
// Implement the interface
/////////////////////////////

template <class Impl>
RealD WilsonFermion<Impl>::M(const FermionField &in, FermionField &out) {
  out.Checkerboard() = in.Checkerboard();
  Dhop(in, out, DaggerNo);
  return axpy_norm(out, 4 + mass, in, out);
}

template <class Impl>
RealD WilsonFermion<Impl>::Mdag(const FermionField &in, FermionField &out) {
  out.Checkerboard() = in.Checkerboard();
  Dhop(in, out, DaggerYes);
  return axpy_norm(out, 4 + mass, in, out);
}

template <class Impl>
void WilsonFermion<Impl>::Meooe(const FermionField &in, FermionField &out) {
  if (in.Checkerboard() == Odd) {
    DhopEO(in, out, DaggerNo);
  } else {
    DhopOE(in, out, DaggerNo);
  }
}

template <class Impl>
void WilsonFermion<Impl>::MeooeDag(const FermionField &in, FermionField &out) {
  if (in.Checkerboard() == Odd) {
    DhopEO(in, out, DaggerYes);
  } else {
    DhopOE(in, out, DaggerYes);
  }
}
  
template <class Impl>
void WilsonFermion<Impl>::Mooee(const FermionField &in, FermionField &out) {
  out.Checkerboard() = in.Checkerboard();
  typename FermionField::scalar_type scal(4.0 + mass);
  out = scal * in;
}

template <class Impl>
void WilsonFermion<Impl>::MooeeDag(const FermionField &in, FermionField &out) {
  out.Checkerboard() = in.Checkerboard();
  Mooee(in, out);
}

template<class Impl>
void WilsonFermion<Impl>::MooeeInv(const FermionField &in, FermionField &out) {
  out.Checkerboard() = in.Checkerboard();
  out = (1.0/(4.0+mass))*in;
}
  
template<class Impl>
void WilsonFermion<Impl>::MooeeInvDag(const FermionField &in, FermionField &out) {
  out.Checkerboard() = in.Checkerboard();
  MooeeInv(in,out);
}
template<class Impl>
void WilsonFermion<Impl>::MomentumSpacePropagator(FermionField &out, const FermionField &in,RealD _m) 
{  
  typedef typename FermionField::vector_type vector_type;
  typedef typename FermionField::scalar_type ScalComplex;
  typedef Lattice<iSinglet<vector_type> > LatComplex;
  
  // what type LatticeComplex 
  conformable(_grid,out.Grid());
  
  Gamma::Algebra Gmu [] = {
    Gamma::Algebra::GammaX,
    Gamma::Algebra::GammaY,
    Gamma::Algebra::GammaZ,
    Gamma::Algebra::GammaT
  };
  
  Coordinate latt_size   = _grid->_fdimensions;
  
  FermionField   num  (_grid); num  = Zero();
  LatComplex    wilson(_grid); wilson= Zero();
  LatComplex     one  (_grid); one = ScalComplex(1.0,0.0);
  
  LatComplex denom(_grid); denom= Zero();
  LatComplex kmu(_grid); 
  ScalComplex ci(0.0,1.0);
  // momphase = n * 2pi / L
  for(int mu=0;mu<Nd;mu++) {
    
    LatticeCoordinate(kmu,mu);
    
    RealD TwoPiL =  M_PI * 2.0/ latt_size[mu];
    
    kmu = TwoPiL * kmu;
    
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

  FermionField Btilde(B.Grid());
  FermionField Atilde(B.Grid());
  Atilde = A;

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
    auto U_v = U.View();
    auto B_v = B.View();
    auto Btilde_v = Btilde.View();
    auto st_v = st.View();
    thread_loop( (int sss = 0; sss < B.Grid()->oSites(); sss++) ,{
      Kernels::DhopDirK(st_v, U_v, st.CommBuf(), sss, sss, B_v, Btilde_v, mu, gamma);
    });

    //////////////////////////////////////////////////
    // spin trace outer product
    //////////////////////////////////////////////////
    Impl::InsertForce4D(mat, Btilde, Atilde, mu);
  }
}

template <class Impl>
void WilsonFermion<Impl>::DhopDeriv(GaugeField &mat, const FermionField &U, const FermionField &V, int dag) {
  conformable(U.Grid(), _grid);
  conformable(U.Grid(), V.Grid());
  conformable(U.Grid(), mat.Grid());

  mat.Checkerboard() = U.Checkerboard();

  DerivInternal(Stencil, Umu, mat, U, V, dag);
}

template <class Impl>
void WilsonFermion<Impl>::DhopDerivOE(GaugeField &mat, const FermionField &U, const FermionField &V, int dag) {
  conformable(U.Grid(), _cbgrid);
  conformable(U.Grid(), V.Grid());
  //conformable(U.Grid(), mat.Grid()); not general, leaving as a comment (Guido)
  // Motivation: look at the SchurDiff operator
  
  assert(V.Checkerboard() == Even);
  assert(U.Checkerboard() == Odd);
  mat.Checkerboard() = Odd;

  DerivInternal(StencilEven, UmuOdd, mat, U, V, dag);
}

template <class Impl>
void WilsonFermion<Impl>::DhopDerivEO(GaugeField &mat, const FermionField &U, const FermionField &V, int dag) {
  conformable(U.Grid(), _cbgrid);
  conformable(U.Grid(), V.Grid());
  //conformable(U.Grid(), mat.Grid());

  assert(V.Checkerboard() == Odd);
  assert(U.Checkerboard() == Even);
  mat.Checkerboard() = Even;

  DerivInternal(StencilOdd, UmuEven, mat, U, V, dag);
}

template <class Impl>
void WilsonFermion<Impl>::Dhop(const FermionField &in, FermionField &out, int dag) {
  conformable(in.Grid(), _grid);  // verifies full grid
  conformable(in.Grid(), out.Grid());

  out.Checkerboard() = in.Checkerboard();

  DhopInternal(Stencil, Lebesgue, Umu, in, out, dag);
}

template <class Impl>
void WilsonFermion<Impl>::DhopOE(const FermionField &in, FermionField &out, int dag) {
  conformable(in.Grid(), _cbgrid);    // verifies half grid
  conformable(in.Grid(), out.Grid());  // drops the cb check

  assert(in.Checkerboard() == Even);
  out.Checkerboard() = Odd;

  DhopInternal(StencilEven, LebesgueEvenOdd, UmuOdd, in, out, dag);
}

template <class Impl>
void WilsonFermion<Impl>::DhopEO(const FermionField &in, FermionField &out,int dag) {
  conformable(in.Grid(), _cbgrid);    // verifies half grid
  conformable(in.Grid(), out.Grid());  // drops the cb check

  assert(in.Checkerboard() == Odd);
  out.Checkerboard() = Even;

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
  auto in_v = in.View();
  auto out_v = in.View();
  auto Umu_v = Umu.View();
  auto Stencil_v = Stencil.View();
  thread_loop( (int sss = 0; sss < in.Grid()->oSites(); sss++) ,{
    Kernels::DhopDirK(Stencil_v, Umu_v, Stencil.CommBuf(), sss, sss, in_v, out_v, dirdisp, gamma);
  });
};

template <class Impl>
void WilsonFermion<Impl>::DhopInternal(StencilImpl &st, LebesgueOrder &lo,
                                       DoubledGaugeField &U,
                                       const FermionField &in,
                                       FermionField &out, int dag) {
  assert((dag == DaggerNo) || (dag == DaggerYes));

  Compressor compressor(dag);
  st.HaloExchange(in, compressor);

  int Opt = WilsonKernelsStatic::Opt;
  auto U_v  = U.View();
  auto in_v = in.View();
  auto out_v= out.View();
  auto st_v = st.View();
  if (dag == DaggerYes) {
    accelerator_loop( sss,in_v, {
      Kernels::DhopSiteDag(Opt,st_v, U_v, st.CommBuf(), sss, sss, 1, 1, in_v, out_v);
    });
  } else {
    accelerator_loop( sss,in_v, {
      Kernels::DhopSite(Opt,st_v, U_v, st.CommBuf(), sss, sss, 1, 1, in_v, out_v);
    });
  }
};

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
  conformable(_grid, q_in_1.Grid());
  conformable(_grid, q_in_2.Grid());
  conformable(_grid, q_out.Grid());
  PropagatorField tmp1(_grid), tmp2(_grid);
  q_out = Zero();

  // Forward, need q1(x + mu), q2(x). Backward, need q1(x), q2(x + mu).
  // Inefficient comms method but not performance critical.
  tmp1 = Cshift(q_in_1, mu, 1);
  tmp2 = Cshift(q_in_2, mu, 1);
  auto tmp1_v  =  tmp1.View();
  auto tmp2_v  =  tmp2.View();
  auto q_in_1_v=q_in_1.View();
  auto q_in_2_v=q_in_2.View();
  auto q_out_v = q_out.View();
  auto Umu_v   =   Umu.View();
  thread_loop( (unsigned int sU = 0; sU < Umu.Grid()->oSites(); ++sU), {
      Kernels::ContractConservedCurrentSiteFwd(tmp1_v[sU],
					       q_in_2_v[sU],
					       q_out_v[sU],
					       Umu_v, sU, mu);
      Kernels::ContractConservedCurrentSiteBwd(q_in_1_v[sU],
					       tmp2_v[sU],
					       q_out_v[sU],
					       Umu_v, sU, mu);
  });
}

template <class Impl>
void WilsonFermion<Impl>::SeqConservedCurrent(PropagatorField &q_in, 
                                              PropagatorField &q_out,
                                              Current curr_type,
                                              unsigned int mu,
                                              std::vector<Real> mom,
                                              unsigned int tmin, 
                                              unsigned int tmax)
{
  conformable(_grid, q_in.Grid());
  conformable(_grid, q_out.Grid());
  Lattice<iSinglet<Simd>> ph(_grid), coor(_grid);
  Complex i(0.0,1.0);
  PropagatorField tmpFwd(_grid), tmpBwd(_grid), tmp(_grid);
  unsigned int tshift = (mu == Tp) ? 1 : 0;
  unsigned int LLt    = GridDefaultLatt()[Tp];

  // Momentum projection
  ph = Zero();
  for(unsigned int mu = 0; mu < Nd - 1; mu++)
    {
      LatticeCoordinate(coor, mu);
      ph = ph + mom[mu]*coor*((1./(_grid->_fdimensions[mu])));
    }
  ph = exp((Real)(2*M_PI)*i*ph);

  q_out = Zero();
  LatticeInteger coords(_grid);
  LatticeCoordinate(coords, Tp);

  // Need q(x + mu) and q(x - mu).
  tmp = Cshift(q_in, mu, 1);
  tmpFwd = tmp*ph;
  tmp = ph*q_in;
  tmpBwd = Cshift(tmp, mu, -1);

  auto coords_v = coords.View();
  auto tmpFwd_v = tmpFwd.View();
  auto tmpBwd_v = tmpBwd.View();
  auto Umu_v    = Umu.View();
  auto q_out_v  = q_out.View();

  thread_loop( (unsigned int sU = 0; sU < Umu.Grid()->oSites(); ++sU), {

    // Compute the sequential conserved current insertion only if our simd
    // object contains a timeslice we need.
    vInteger t_mask   = ((coords_v[sU] >= tmin) &&
			 (coords_v[sU] <= tmax));
    Integer timeSlices = Reduce(t_mask);

    if (timeSlices > 0) {
      Kernels::SeqConservedCurrentSiteFwd(tmpFwd_v[sU], 
					  q_out_v[sU], 
					  Umu_v, sU, mu, t_mask);
    }

    // Repeat for backward direction.
    t_mask     = ((coords_v[sU] >= (tmin + tshift)) && 
		  (coords_v[sU] <= (tmax + tshift)));
    
    //if tmax = LLt-1 (last timeslice) include timeslice 0 if the time is shifted (mu=3)	
    unsigned int t0 = 0;
    if((tmax==LLt-1) && (tshift==1)) t_mask = (t_mask || (coords_v[sU] == t0 ));
    
    timeSlices = Reduce(t_mask);

    if (timeSlices > 0) {
      Kernels::SeqConservedCurrentSiteBwd(tmpBwd_v[sU], 
					  q_out_v[sU], 
					  Umu_v, sU, mu, t_mask);
    }
  });
}

FermOpTemplateInstantiate(WilsonFermion);
AdjointFermOpTemplateInstantiate(WilsonFermion);
TwoIndexFermOpTemplateInstantiate(WilsonFermion);
GparityFermOpTemplateInstantiate(WilsonFermion);

NAMESPACE_END(Grid);
