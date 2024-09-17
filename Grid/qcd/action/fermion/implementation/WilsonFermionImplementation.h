/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/action/fermion/WilsonFermion.cc

Copyright (C) 2022

Author: Peter Boyle <pabobyle@ph.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: Peter Boyle <peterboyle@Peters-MacBook-Pro-2.local>
Author: paboyle <paboyle@ph.ed.ac.uk>
Author: Fabian Joswig <fabian.joswig@ed.ac.uk>

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

/////////////////////////////////
// Constructor and gauge import
/////////////////////////////////

template <class Impl>
WilsonFermion<Impl>::WilsonFermion(GaugeField &_Umu, GridCartesian &Fgrid,
                                   GridRedBlackCartesian &Hgrid, RealD _mass,
                                   const ImplParams &p,
                                   const WilsonAnisotropyCoefficients &anis)
  :
    Kernels(p),
    _grid(&Fgrid),
    _cbgrid(&Hgrid),
    Stencil(&Fgrid, npoint, Even, directions, displacements,p),
    StencilEven(&Hgrid, npoint, Even, directions,displacements,p),  // source is Even
    StencilOdd(&Hgrid, npoint, Odd, directions,displacements,p),  // source is Odd
    mass(_mass),
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

  int vol4;
  vol4=Fgrid.oSites();
  Stencil.BuildSurfaceList(1,vol4);
  vol4=Hgrid.oSites();
  StencilEven.BuildSurfaceList(1,vol4);
  StencilOdd.BuildSurfaceList(1,vol4);
}

template <class Impl>
void WilsonFermion<Impl>::ImportGauge(const GaugeField &_Umu)
{
  GaugeField HUmu(_Umu.Grid());

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
void WilsonFermion<Impl>::M(const FermionField &in, FermionField &out)
{
  out.Checkerboard() = in.Checkerboard();
  Dhop(in, out, DaggerNo);
  axpy(out, diag_mass, in, out);
}

template <class Impl>
void WilsonFermion<Impl>::Mdag(const FermionField &in, FermionField &out)
{
  out.Checkerboard() = in.Checkerboard();
  Dhop(in, out, DaggerYes);
  axpy(out, diag_mass, in, out);
}

template <class Impl>
void WilsonFermion<Impl>::Meooe(const FermionField &in, FermionField &out)
{
  if (in.Checkerboard() == Odd) {
    DhopEO(in, out, DaggerNo);
  } else {
    DhopOE(in, out, DaggerNo);
  }
}

template <class Impl>
void WilsonFermion<Impl>::MeooeDag(const FermionField &in, FermionField &out)
{
  if (in.Checkerboard() == Odd) {
    DhopEO(in, out, DaggerYes);
  } else {
    DhopOE(in, out, DaggerYes);
  }
}

template <class Impl>
void WilsonFermion<Impl>::Mooee(const FermionField &in, FermionField &out)
{
  out.Checkerboard() = in.Checkerboard();
  typename FermionField::scalar_type scal(diag_mass);
  out = scal * in;
}

template <class Impl>
void WilsonFermion<Impl>::MooeeDag(const FermionField &in, FermionField &out)
{
  out.Checkerboard() = in.Checkerboard();
  Mooee(in, out);
}

template<class Impl>
void WilsonFermion<Impl>::MooeeInv(const FermionField &in, FermionField &out)
{
  out.Checkerboard() = in.Checkerboard();
  out = (1.0/(diag_mass))*in;
}

template<class Impl>
void WilsonFermion<Impl>::MooeeInvDag(const FermionField &in, FermionField &out)
{
  out.Checkerboard() = in.Checkerboard();
  MooeeInv(in,out);
}
template<class Impl>
void WilsonFermion<Impl>::MomentumSpacePropagator(FermionField &out, const FermionField &in,RealD _m,std::vector<double> twist)
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

    int Ls=1;
    Kernels::DhopDirKernel(st, U, st.CommBuf(), Ls, B.Grid()->oSites(), B, Btilde, mu, gamma);

    //////////////////////////////////////////////////
    // spin trace outer product
    //////////////////////////////////////////////////
    Impl::InsertForce4D(mat, Btilde, Atilde, mu);
  }
}

template <class Impl>
void WilsonFermion<Impl>::DhopDeriv(GaugeField &mat, const FermionField &U, const FermionField &V, int dag)
{
  conformable(U.Grid(), _grid);
  conformable(U.Grid(), V.Grid());
  conformable(U.Grid(), mat.Grid());

  mat.Checkerboard() = U.Checkerboard();

  DerivInternal(Stencil, Umu, mat, U, V, dag);
}

template <class Impl>
void WilsonFermion<Impl>::DhopDerivOE(GaugeField &mat, const FermionField &U, const FermionField &V, int dag)
{
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
void WilsonFermion<Impl>::DhopDerivEO(GaugeField &mat, const FermionField &U, const FermionField &V, int dag)
{
  conformable(U.Grid(), _cbgrid);
  conformable(U.Grid(), V.Grid());
  //conformable(U.Grid(), mat.Grid());

  assert(V.Checkerboard() == Odd);
  assert(U.Checkerboard() == Even);
  mat.Checkerboard() = Even;

  DerivInternal(StencilOdd, UmuEven, mat, U, V, dag);
}

template <class Impl>
void WilsonFermion<Impl>::Dhop(const FermionField &in, FermionField &out, int dag)
{
  conformable(in.Grid(), _grid);  // verifies full grid
  conformable(in.Grid(), out.Grid());

  out.Checkerboard() = in.Checkerboard();

  DhopInternal(Stencil, Umu, in, out, dag);
}

template <class Impl>
void WilsonFermion<Impl>::DhopOE(const FermionField &in, FermionField &out, int dag)
{
  conformable(in.Grid(), _cbgrid);    // verifies half grid
  conformable(in.Grid(), out.Grid());  // drops the cb check

  assert(in.Checkerboard() == Even);
  out.Checkerboard() = Odd;

  DhopInternal(StencilEven, UmuOdd, in, out, dag);
}

template <class Impl>
void WilsonFermion<Impl>::DhopEO(const FermionField &in, FermionField &out,int dag)
{
  conformable(in.Grid(), _cbgrid);    // verifies half grid
  conformable(in.Grid(), out.Grid());  // drops the cb check

  assert(in.Checkerboard() == Odd);
  out.Checkerboard() = Even;

  DhopInternal(StencilOdd, UmuEven, in, out, dag);
}

template <class Impl>
void WilsonFermion<Impl>::Mdir(const FermionField &in, FermionField &out, int dir, int disp)
{
  DhopDir(in, out, dir, disp);
}
template <class Impl>
void WilsonFermion<Impl>::MdirAll(const FermionField &in, std::vector<FermionField> &out)
{
  DhopDirAll(in, out);
}

template <class Impl>
void WilsonFermion<Impl>::DhopDir(const FermionField &in, FermionField &out, int dir, int disp)
{
  Compressor compressor(DaggerNo);
  Stencil.HaloExchange(in, compressor);

  int skip = (disp == 1) ? 0 : 1;
  int dirdisp = dir + skip * 4;
  int gamma = dir + (1 - skip) * 4;

  DhopDirCalc(in, out, dirdisp, gamma, DaggerNo);
};
template <class Impl>
void WilsonFermion<Impl>::DhopDirAll(const FermionField &in, std::vector<FermionField> &out)
{
  Compressor compressor(DaggerNo);
  Stencil.HaloExchange(in, compressor);

  assert((out.size()==8)||(out.size()==9));
  for(int dir=0;dir<Nd;dir++){
    for(int disp=-1;disp<=1;disp+=2){

      int skip = (disp == 1) ? 0 : 1;
      int dirdisp = dir + skip * 4;
      int gamma = dir + (1 - skip) * 4;

      DhopDirCalc(in, out[dirdisp], dirdisp, gamma, DaggerNo);
    }
  }
}
template <class Impl>
void WilsonFermion<Impl>::DhopDirCalc(const FermionField &in, FermionField &out,int dirdisp, int gamma, int dag)
{
  int Ls=1;
  uint64_t Nsite=in.oSites();
  Kernels::DhopDirKernel(Stencil, Umu, Stencil.CommBuf(), Ls, Nsite, in, out, dirdisp, gamma);
};

template <class Impl>
void WilsonFermion<Impl>::DhopInternal(StencilImpl &st, 
                                       DoubledGaugeField &U,
                                       const FermionField &in,
                                       FermionField &out, int dag)
{
#ifdef GRID_OMP
  if ( WilsonKernelsStatic::Comms == WilsonKernelsStatic::CommsAndCompute )
    DhopInternalOverlappedComms(st,U,in,out,dag);
  else
#endif
    DhopInternalSerial(st,U,in,out,dag);
}

template <class Impl>
void WilsonFermion<Impl>::DhopInternalOverlappedComms(StencilImpl &st, 
						      DoubledGaugeField &U,
						      const FermionField &in,
						      FermionField &out, int dag)
{
  GRID_TRACE("DhopOverlapped");
  assert((dag == DaggerNo) || (dag == DaggerYes));

  Compressor compressor(dag);
  int len =  U.Grid()->oSites();

  /////////////////////////////
  // Start comms  // Gather intranode and extra node differentiated??
  /////////////////////////////
  std::vector<std::vector<CommsRequest_t> > requests;
  st.Prepare();
  {
    GRID_TRACE("Gather");
    st.HaloGather(in,compressor);
  }

  tracePush("Communication");
  st.CommunicateBegin(requests);

  /////////////////////////////
  // Overlap with comms
  /////////////////////////////
  {
    GRID_TRACE("MergeSHM");
    st.CommsMergeSHM(compressor);
  }

  /////////////////////////////
  // do the compute interior
  /////////////////////////////
  int Opt = WilsonKernelsStatic::Opt;
  if (dag == DaggerYes) {
    GRID_TRACE("DhopDagInterior");
    Kernels::DhopDagKernel(Opt,st,U,st.CommBuf(),1,U.oSites(),in,out,1,0);
  } else {
    GRID_TRACE("DhopInterior");
    Kernels::DhopKernel(Opt,st,U,st.CommBuf(),1,U.oSites(),in,out,1,0);
  }

  /////////////////////////////
  // Complete comms
  /////////////////////////////
  st.CommunicateComplete(requests);
  tracePop("Communication");

  {
    GRID_TRACE("Merge");
    st.CommsMerge(compressor);
  }
  /////////////////////////////
  // do the compute exterior
  /////////////////////////////

  if (dag == DaggerYes) {
    GRID_TRACE("DhopDagExterior");
    Kernels::DhopDagKernel(Opt,st,U,st.CommBuf(),1,U.oSites(),in,out,0,1);
  } else {
    GRID_TRACE("DhopExterior");
    Kernels::DhopKernel(Opt,st,U,st.CommBuf(),1,U.oSites(),in,out,0,1);
  }
};


template <class Impl>
void WilsonFermion<Impl>::DhopInternalSerial(StencilImpl &st, 
					     DoubledGaugeField &U,
					     const FermionField &in,
					     FermionField &out, int dag)
{
  GRID_TRACE("DhopSerial");
  assert((dag == DaggerNo) || (dag == DaggerYes));
  Compressor compressor(dag);
  {
    GRID_TRACE("HaloExchange");
    st.HaloExchange(in, compressor);
  }

  int Opt = WilsonKernelsStatic::Opt;
  if (dag == DaggerYes) {
    GRID_TRACE("DhopDag");
    Kernels::DhopDagKernel(Opt,st,U,st.CommBuf(),1,U.oSites(),in,out);
  } else {
    GRID_TRACE("Dhop");
    Kernels::DhopKernel(Opt,st,U,st.CommBuf(),1,U.oSites(),in,out);
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
                                                   PropagatorField &src,
                                                   Current curr_type,
                                                   unsigned int mu)
{
  if(curr_type != Current::Vector)
  {
    std::cout << GridLogError << "Only the conserved vector current is implemented so far." << std::endl;
    exit(1);
  }

  Gamma g5(Gamma::Algebra::Gamma5);
  conformable(_grid, q_in_1.Grid());
  conformable(_grid, q_in_2.Grid());
  conformable(_grid, q_out.Grid());
  auto UGrid= this->GaugeGrid();

  PropagatorField tmp_shifted(UGrid);
  PropagatorField g5Lg5(UGrid);
  PropagatorField R(UGrid);
  PropagatorField gmuR(UGrid);

    Gamma::Algebra Gmu [] = {
    Gamma::Algebra::GammaX,
    Gamma::Algebra::GammaY,
    Gamma::Algebra::GammaZ,
    Gamma::Algebra::GammaT,
  };
  Gamma gmu=Gamma(Gmu[mu]);

  g5Lg5=g5*q_in_1*g5;
  tmp_shifted=Cshift(q_in_2,mu,1);
  Impl::multLinkField(R,this->Umu,tmp_shifted,mu);
  gmuR=gmu*R;

  q_out=adj(g5Lg5)*R;
  q_out-=adj(g5Lg5)*gmuR;

  tmp_shifted=Cshift(q_in_1,mu,1);
  Impl::multLinkField(g5Lg5,this->Umu,tmp_shifted,mu);
  g5Lg5=g5*g5Lg5*g5;
  R=q_in_2;
  gmuR=gmu*R;

  q_out-=adj(g5Lg5)*R;
  q_out-=adj(g5Lg5)*gmuR;
}


template <class Impl>
void WilsonFermion<Impl>::SeqConservedCurrent(PropagatorField &q_in,
                                              PropagatorField &q_out,
                                              PropagatorField &src,
                                              Current curr_type,
                                              unsigned int mu,
                                              unsigned int tmin,
                                              unsigned int tmax,
					      ComplexField &lattice_cmplx)
{
  if(curr_type != Current::Vector)
  {
    std::cout << GridLogError << "Only the conserved vector current is implemented so far." << std::endl;
    exit(1);
  }

  int tshift = (mu == Nd-1) ? 1 : 0;
  unsigned int LLt    = GridDefaultLatt()[Tp];
  conformable(_grid, q_in.Grid());
  conformable(_grid, q_out.Grid());
  auto UGrid= this->GaugeGrid();

  PropagatorField tmp(UGrid);
  PropagatorField Utmp(UGrid);
  PropagatorField L(UGrid);
  PropagatorField zz (UGrid);
  zz=Zero();
  LatticeInteger lcoor(UGrid); LatticeCoordinate(lcoor,Nd-1);

    Gamma::Algebra Gmu [] = {
    Gamma::Algebra::GammaX,
    Gamma::Algebra::GammaY,
    Gamma::Algebra::GammaZ,
    Gamma::Algebra::GammaT,
  };
  Gamma gmu=Gamma(Gmu[mu]);

  tmp = Cshift(q_in,mu,1);
  Impl::multLinkField(Utmp,this->Umu,tmp,mu);
  tmp = ( Utmp*lattice_cmplx - gmu*Utmp*lattice_cmplx ); // Forward hop
  tmp = where((lcoor>=tmin),tmp,zz); // Mask the time
  q_out = where((lcoor<=tmax),tmp,zz); // Position of current complicated

  tmp = q_in *lattice_cmplx;
  tmp = Cshift(tmp,mu,-1);
  Impl::multLinkField(Utmp,this->Umu,tmp,mu+Nd); // Adjoint link
  tmp = -( Utmp + gmu*Utmp );
  // Mask the time
  if (tmax == LLt - 1 && tshift == 1){ // quick fix to include timeslice 0 if tmax + tshift is over the last timeslice
    unsigned int t0 = 0;
    tmp = where(((lcoor==t0) || (lcoor>=tmin+tshift)),tmp,zz);
  } else {
    tmp = where((lcoor>=tmin+tshift),tmp,zz);
  }
  q_out+= where((lcoor<=tmax+tshift),tmp,zz); // Position of current complicated
}

NAMESPACE_END(Grid);
