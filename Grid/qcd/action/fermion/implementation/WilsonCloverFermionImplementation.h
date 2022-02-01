/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./lib/qcd/action/fermion/WilsonCloverFermionImplementation.h

    Copyright (C) 2017 - 2022

    Author: paboyle <paboyle@ph.ed.ac.uk>
    Author: Guido Cossu <guido.cossu@ed.ac.uk>
    Author: Daniel Richtmann <daniel.richtmann@gmail.com>

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

#include <Grid/Grid.h>
#include <Grid/qcd/spin/Dirac.h>
#include <Grid/qcd/action/fermion/WilsonCloverFermion.h>

NAMESPACE_BEGIN(Grid);

template<class Impl>
WilsonCloverFermion<Impl>::WilsonCloverFermion(GaugeField&                         _Umu,
                                               GridCartesian&                      Fgrid,
                                               GridRedBlackCartesian&              Hgrid,
                                               const RealD                         _mass,
                                               const RealD                         _csw_r,
                                               const RealD                         _csw_t,
                                               const WilsonAnisotropyCoefficients& clover_anisotropy,
                                               const ImplParams&                   impl_p)
  : WilsonFermion<Impl>(_Umu, Fgrid, Hgrid, _mass, impl_p, clover_anisotropy)
  , CloverTerm(&Fgrid)
  , CloverTermInv(&Fgrid)
  , CloverTermEven(&Hgrid)
  , CloverTermOdd(&Hgrid)
  , CloverTermInvEven(&Hgrid)
  , CloverTermInvOdd(&Hgrid)
  , CloverTermDagEven(&Hgrid)
  , CloverTermDagOdd(&Hgrid)
  , CloverTermInvDagEven(&Hgrid)
  , CloverTermInvDagOdd(&Hgrid) {
  assert(Nd == 4); // require 4 dimensions

  if(clover_anisotropy.isAnisotropic) {
    csw_r     = _csw_r * 0.5 / clover_anisotropy.xi_0;
    diag_mass = _mass + 1.0 + (Nd - 1) * (clover_anisotropy.nu / clover_anisotropy.xi_0);
  } else {
    csw_r     = _csw_r * 0.5;
    diag_mass = 4.0 + _mass;
  }
  csw_t = _csw_t * 0.5;

  if(csw_r == 0)
    std::cout << GridLogWarning << "Initializing WilsonCloverFermion with csw_r = 0" << std::endl;
  if(csw_t == 0)
    std::cout << GridLogWarning << "Initializing WilsonCloverFermion with csw_t = 0" << std::endl;

  ImportGauge(_Umu);
}

// *NOT* EO
template <class Impl>
void WilsonCloverFermion<Impl>::M(const FermionField &in, FermionField &out)
{
  FermionField temp(out.Grid());

  // Wilson term
  out.Checkerboard() = in.Checkerboard();
  this->Dhop(in, out, DaggerNo);

  // Clover term
  Mooee(in, temp);

  out += temp;
}

template <class Impl>
void WilsonCloverFermion<Impl>::Mdag(const FermionField &in, FermionField &out)
{
  FermionField temp(out.Grid());

  // Wilson term
  out.Checkerboard() = in.Checkerboard();
  this->Dhop(in, out, DaggerYes);

  // Clover term
  MooeeDag(in, temp);

  out += temp;
}

template <class Impl>
void WilsonCloverFermion<Impl>::ImportGauge(const GaugeField &_Umu)
{
  double t0 = usecond();
  WilsonFermion<Impl>::ImportGauge(_Umu);
  double t1 = usecond();
  GridBase *grid = _Umu.Grid();
  typename Impl::GaugeLinkField Bx(grid), By(grid), Bz(grid), Ex(grid), Ey(grid), Ez(grid);

  double t2 = usecond();
  // Compute the field strength terms mu>nu
  WilsonLoops<Impl>::FieldStrength(Bx, _Umu, Zdir, Ydir);
  WilsonLoops<Impl>::FieldStrength(By, _Umu, Zdir, Xdir);
  WilsonLoops<Impl>::FieldStrength(Bz, _Umu, Ydir, Xdir);
  WilsonLoops<Impl>::FieldStrength(Ex, _Umu, Tdir, Xdir);
  WilsonLoops<Impl>::FieldStrength(Ey, _Umu, Tdir, Ydir);
  WilsonLoops<Impl>::FieldStrength(Ez, _Umu, Tdir, Zdir);

  double t3 = usecond();
  // Compute the Clover Operator acting on Colour and Spin
  // multiply here by the clover coefficients for the anisotropy
  CloverTerm  = Helpers::fillCloverYZ(Bx) * csw_r;
  CloverTerm += Helpers::fillCloverXZ(By) * csw_r;
  CloverTerm += Helpers::fillCloverXY(Bz) * csw_r;
  CloverTerm += Helpers::fillCloverXT(Ex) * csw_t;
  CloverTerm += Helpers::fillCloverYT(Ey) * csw_t;
  CloverTerm += Helpers::fillCloverZT(Ez) * csw_t;
  CloverTerm += diag_mass;

  double t4 = usecond();
  int lvol = _Umu.Grid()->lSites();
  int DimRep = Impl::Dimension;

  double t5 = usecond();
  {
    autoView(CTv,CloverTerm,CpuRead);
    autoView(CTIv,CloverTermInv,CpuWrite);
    thread_for(site, lvol, {
      Coordinate lcoor;
      grid->LocalIndexToLocalCoor(site, lcoor);
      Eigen::MatrixXcd EigenCloverOp = Eigen::MatrixXcd::Zero(Ns * DimRep, Ns * DimRep);
      Eigen::MatrixXcd EigenInvCloverOp = Eigen::MatrixXcd::Zero(Ns * DimRep, Ns * DimRep);
      typename SiteClover::scalar_object Qx = Zero(), Qxinv = Zero();
      peekLocalSite(Qx, CTv, lcoor);
      //if (csw!=0){
      for (int j = 0; j < Ns; j++)
	for (int k = 0; k < Ns; k++)
	  for (int a = 0; a < DimRep; a++)
	    for (int b = 0; b < DimRep; b++){
	      auto zz =  Qx()(j, k)(a, b);
	      EigenCloverOp(a + j * DimRep, b + k * DimRep) = std::complex<double>(zz);
	    }
      //   if (site==0) std::cout << "site =" << site << "\n" << EigenCloverOp << std::endl;
      
      EigenInvCloverOp = EigenCloverOp.inverse();
      //std::cout << EigenInvCloverOp << std::endl;
      for (int j = 0; j < Ns; j++)
	for (int k = 0; k < Ns; k++)
	  for (int a = 0; a < DimRep; a++)
	    for (int b = 0; b < DimRep; b++)
	      Qxinv()(j, k)(a, b) = EigenInvCloverOp(a + j * DimRep, b + k * DimRep);
      //    if (site==0) std::cout << "site =" << site << "\n" << EigenInvCloverOp << std::endl;
      //  }
      pokeLocalSite(Qxinv, CTIv, lcoor);
    });
  }

  double t6 = usecond();
  // Separate the even and odd parts
  pickCheckerboard(Even, CloverTermEven, CloverTerm);
  pickCheckerboard(Odd, CloverTermOdd, CloverTerm);

  pickCheckerboard(Even, CloverTermDagEven, adj(CloverTerm));
  pickCheckerboard(Odd, CloverTermDagOdd, adj(CloverTerm));

  pickCheckerboard(Even, CloverTermInvEven, CloverTermInv);
  pickCheckerboard(Odd, CloverTermInvOdd, CloverTermInv);

  pickCheckerboard(Even, CloverTermInvDagEven, adj(CloverTermInv));
  pickCheckerboard(Odd, CloverTermInvDagOdd, adj(CloverTermInv));
  double t7 = usecond();

#if 0
  std::cout << GridLogMessage << "WilsonCloverFermion::ImportGauge timings:"
            << " WilsonFermion::Importgauge = " << (t1 - t0) / 1e6
            << ", allocations = "               << (t2 - t1) / 1e6
            << ", field strength = "            << (t3 - t2) / 1e6
            << ", fill clover = "               << (t4 - t3) / 1e6
            << ", misc = "                      << (t5 - t4) / 1e6
            << ", inversions = "                << (t6 - t5) / 1e6
            << ", pick cbs = "                  << (t7 - t6) / 1e6
            << ", total = "                     << (t7 - t0) / 1e6
            << std::endl;
#endif
}

template <class Impl>
void WilsonCloverFermion<Impl>::Mooee(const FermionField &in, FermionField &out)
{
  this->MooeeInternal(in, out, DaggerNo, InverseNo);
}

template <class Impl>
void WilsonCloverFermion<Impl>::MooeeDag(const FermionField &in, FermionField &out)
{
  this->MooeeInternal(in, out, DaggerYes, InverseNo);
}

template <class Impl>
void WilsonCloverFermion<Impl>::MooeeInv(const FermionField &in, FermionField &out)
{
  this->MooeeInternal(in, out, DaggerNo, InverseYes);
}

template <class Impl>
void WilsonCloverFermion<Impl>::MooeeInvDag(const FermionField &in, FermionField &out)
{
  this->MooeeInternal(in, out, DaggerYes, InverseYes);
}

template <class Impl>
void WilsonCloverFermion<Impl>::MooeeInternal(const FermionField &in, FermionField &out, int dag, int inv)
{
  out.Checkerboard() = in.Checkerboard();
  CloverField *Clover;
  assert(in.Checkerboard() == Odd || in.Checkerboard() == Even);

  if (dag)
  {
    if (in.Grid()->_isCheckerBoarded)
    {
      if (in.Checkerboard() == Odd)
      {
        Clover = (inv) ? &CloverTermInvDagOdd : &CloverTermDagOdd;
      }
      else
      {
        Clover = (inv) ? &CloverTermInvDagEven : &CloverTermDagEven;
      }
      Helpers::multCloverField(out, *Clover, in);
    }
    else
    {
      Clover = (inv) ? &CloverTermInv : &CloverTerm;
      Helpers::multCloverField(out, *Clover, in); // don't bother with adj, hermitian anyway
    }
  }
  else
  {
    if (in.Grid()->_isCheckerBoarded)
    {

      if (in.Checkerboard() == Odd)
      {
        //  std::cout << "Calling clover term Odd" << std::endl;
        Clover = (inv) ? &CloverTermInvOdd : &CloverTermOdd;
      }
      else
      {
        //  std::cout << "Calling clover term Even" << std::endl;
        Clover = (inv) ? &CloverTermInvEven : &CloverTermEven;
      }
      Helpers::multCloverField(out, *Clover, in);
      //  std::cout << GridLogMessage << "*Clover.Checkerboard() "  << (*Clover).Checkerboard() << std::endl;
    }
    else
    {
      Clover = (inv) ? &CloverTermInv : &CloverTerm;
      Helpers::multCloverField(out, *Clover, in);
    }
  }
} // MooeeInternal

// Derivative parts unpreconditioned pseudofermions
template <class Impl>
void WilsonCloverFermion<Impl>::MDeriv(GaugeField &force, const FermionField &X, const FermionField &Y, int dag)
{
  conformable(X.Grid(), Y.Grid());
  conformable(X.Grid(), force.Grid());
  GaugeLinkField force_mu(force.Grid()), lambda(force.Grid());
  GaugeField clover_force(force.Grid());
  PropagatorField Lambda(force.Grid());

  // Guido: Here we are hitting some performance issues:
  // need to extract the components of the DoubledGaugeField
  // for each call
  // Possible solution
  // Create a vector object to store them? (cons: wasting space)
  std::vector<GaugeLinkField> U(Nd, this->Umu.Grid());

  Impl::extractLinkField(U, this->Umu);

  force = Zero();
  // Derivative of the Wilson hopping term
  this->DhopDeriv(force, X, Y, dag);

  ///////////////////////////////////////////////////////////
  // Clover term derivative
  ///////////////////////////////////////////////////////////
  Impl::outerProductImpl(Lambda, X, Y);
  //std::cout << "Lambda:" << Lambda << std::endl;

  Gamma::Algebra sigma[] = {
      Gamma::Algebra::SigmaXY,
      Gamma::Algebra::SigmaXZ,
      Gamma::Algebra::SigmaXT,
      Gamma::Algebra::MinusSigmaXY,
      Gamma::Algebra::SigmaYZ,
      Gamma::Algebra::SigmaYT,
      Gamma::Algebra::MinusSigmaXZ,
      Gamma::Algebra::MinusSigmaYZ,
      Gamma::Algebra::SigmaZT,
      Gamma::Algebra::MinusSigmaXT,
      Gamma::Algebra::MinusSigmaYT,
      Gamma::Algebra::MinusSigmaZT};

  /*
    sigma_{\mu \nu}=
    | 0         sigma[0]  sigma[1]  sigma[2] |
    | sigma[3]    0       sigma[4]  sigma[5] |
    | sigma[6]  sigma[7]     0      sigma[8] |
    | sigma[9]  sigma[10] sigma[11]   0      |
  */

  int count = 0;
  clover_force = Zero();
  for (int mu = 0; mu < 4; mu++)
  {
    force_mu = Zero();
    for (int nu = 0; nu < 4; nu++)
    {
      if (mu == nu)
      continue;

      RealD factor;
      if (nu == 4 || mu == 4)
      {
        factor = 2.0 * csw_t;
      }
      else
      {
        factor = 2.0 * csw_r;
      }
      PropagatorField Slambda = Gamma(sigma[count]) * Lambda; // sigma checked
      Impl::TraceSpinImpl(lambda, Slambda);                   // traceSpin ok
      force_mu -= factor*Helpers::Cmunu(U, lambda, mu, nu);                   // checked
      count++;
    }

    pokeLorentz(clover_force, U[mu] * force_mu, mu);
  }
  //clover_force *= csw;
  force += clover_force;
}

// Derivative parts
template <class Impl>
void WilsonCloverFermion<Impl>::MooDeriv(GaugeField &mat, const FermionField &X, const FermionField &Y, int dag)
{
  assert(0);
}

// Derivative parts
template <class Impl>
void WilsonCloverFermion<Impl>::MeeDeriv(GaugeField &mat, const FermionField &U, const FermionField &V, int dag)
{
  assert(0); // not implemented yet
}

NAMESPACE_END(Grid);
