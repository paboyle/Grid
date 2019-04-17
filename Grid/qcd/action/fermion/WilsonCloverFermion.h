/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./lib/qcd/action/fermion/WilsonCloverFermion.h

    Copyright (C) 2017

    Author: Guido Cossu <guido.cossu@ed.ac.uk>
    Author: David Preti <>

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

#ifndef GRID_QCD_WILSON_CLOVER_FERMION_H
#define GRID_QCD_WILSON_CLOVER_FERMION_H

#include <Grid/Grid.h>

namespace Grid
{
namespace QCD
{

///////////////////////////////////////////////////////////////////
// Wilson Clover
//
// Operator ( with anisotropy coefficients):
//
// Q =   1 + (Nd-1)/xi_0 + m
//     + W_t + (nu/xi_0) * W_s
//     - 1/2*[ csw_t * sum_s (sigma_ts F_ts) + (csw_s/xi_0) * sum_ss (sigma_ss F_ss)  ]
//
// s spatial, t temporal directions.
// where W_t and W_s are the temporal and spatial components of the
// Wilson Dirac operator
//
// csw_r = csw_t to recover the isotropic version
//////////////////////////////////////////////////////////////////

template <class Impl>
class WilsonCloverFermion : public WilsonFermion<Impl>
{
public:
  // Types definitions
  INHERIT_IMPL_TYPES(Impl);
  template <typename vtype>
  using iImplClover = iScalar<iMatrix<iMatrix<vtype, Impl::Dimension>, Ns>>;
  typedef iImplClover<Simd> SiteCloverType;
  typedef Lattice<SiteCloverType> CloverFieldType;

public:
  typedef WilsonFermion<Impl> WilsonBase;

  virtual int    ConstEE(void)     { return 0; };
  virtual void Instantiatable(void){};
  // Constructors
  WilsonCloverFermion(GaugeField &_Umu, GridCartesian &Fgrid,
                      GridRedBlackCartesian &Hgrid,
                      const RealD _mass,
                      const RealD _csw_r = 0.0,
                      const RealD _csw_t = 0.0,
                      const WilsonAnisotropyCoefficients &clover_anisotropy = WilsonAnisotropyCoefficients(),
                      const ImplParams &impl_p = ImplParams()) : WilsonFermion<Impl>(_Umu,
                                                                                     Fgrid,
                                                                                     Hgrid,
                                                                                     _mass, impl_p, clover_anisotropy),
                                                                 CloverTerm(&Fgrid),
                                                                 CloverTermInv(&Fgrid),
                                                                 CloverTermEven(&Hgrid),
                                                                 CloverTermOdd(&Hgrid),
                                                                 CloverTermInvEven(&Hgrid),
                                                                 CloverTermInvOdd(&Hgrid),
                                                                 CloverTermDagEven(&Hgrid),
                                                                 CloverTermDagOdd(&Hgrid),
                                                                 CloverTermInvDagEven(&Hgrid),
                                                                 CloverTermInvDagOdd(&Hgrid)
  {
    assert(Nd == 4); // require 4 dimensions

    if (clover_anisotropy.isAnisotropic)
    {
      csw_r = _csw_r * 0.5 / clover_anisotropy.xi_0;
      diag_mass = _mass + 1.0 + (Nd - 1) * (clover_anisotropy.nu / clover_anisotropy.xi_0);
    }
    else
    {
      csw_r = _csw_r * 0.5;
      diag_mass = 4.0 + _mass;
    }
    csw_t = _csw_t * 0.5;

    if (csw_r == 0)
      std::cout << GridLogWarning << "Initializing WilsonCloverFermion with csw_r = 0" << std::endl;
    if (csw_t == 0)
      std::cout << GridLogWarning << "Initializing WilsonCloverFermion with csw_t = 0" << std::endl;

    ImportGauge(_Umu);
  }

  virtual RealD M(const FermionField &in, FermionField &out);
  virtual RealD Mdag(const FermionField &in, FermionField &out);

  virtual void Mooee(const FermionField &in, FermionField &out);
  virtual void MooeeDag(const FermionField &in, FermionField &out);
  virtual void MooeeInv(const FermionField &in, FermionField &out);
  virtual void MooeeInvDag(const FermionField &in, FermionField &out);
  virtual void MooeeInternal(const FermionField &in, FermionField &out, int dag, int inv);

  //virtual void MDeriv(GaugeField &mat, const FermionField &U, const FermionField &V, int dag);
  virtual void MooDeriv(GaugeField &mat, const FermionField &U, const FermionField &V, int dag);
  virtual void MeeDeriv(GaugeField &mat, const FermionField &U, const FermionField &V, int dag);

  void ImportGauge(const GaugeField &_Umu);

  // Derivative parts unpreconditioned pseudofermions
  void MDeriv(GaugeField &force, const FermionField &X, const FermionField &Y, int dag)
  {
    conformable(X._grid, Y._grid);
    conformable(X._grid, force._grid);
    GaugeLinkField force_mu(force._grid), lambda(force._grid);
    GaugeField clover_force(force._grid);
    PropagatorField Lambda(force._grid);

    // Guido: Here we are hitting some performance issues:
    // need to extract the components of the DoubledGaugeField
    // for each call
    // Possible solution
    // Create a vector object to store them? (cons: wasting space)
    std::vector<GaugeLinkField> U(Nd, this->Umu._grid);

    Impl::extractLinkField(U, this->Umu);

    force = zero;
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
    clover_force = zero;
    for (int mu = 0; mu < 4; mu++)
    {
      force_mu = zero;
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
        force_mu -= factor*Cmunu(U, lambda, mu, nu);                   // checked
        count++;
      }

      pokeLorentz(clover_force, U[mu] * force_mu, mu);
    }
    //clover_force *= csw;
    force += clover_force;
  }

  // Computing C_{\mu \nu}(x) as in Eq.(B.39) in Zbigniew Sroczynski's PhD thesis
  GaugeLinkField Cmunu(std::vector<GaugeLinkField> &U, GaugeLinkField &lambda, int mu, int nu)
  {
    conformable(lambda._grid, U[0]._grid);
    GaugeLinkField out(lambda._grid), tmp(lambda._grid);
    // insertion in upper staple
    // please check redundancy of shift operations

    // C1+
    tmp = lambda * U[nu];
    out = Impl::ShiftStaple(Impl::CovShiftForward(tmp, nu, Impl::CovShiftBackward(U[mu], mu, Impl::CovShiftIdentityBackward(U[nu], nu))), mu);

    // C2+
    tmp = U[mu] * Impl::ShiftStaple(adj(lambda), mu);
    out += Impl::ShiftStaple(Impl::CovShiftForward(U[nu], nu, Impl::CovShiftBackward(tmp, mu, Impl::CovShiftIdentityBackward(U[nu], nu))), mu);

    // C3+
    tmp = U[nu] * Impl::ShiftStaple(adj(lambda), nu);
    out += Impl::ShiftStaple(Impl::CovShiftForward(U[nu], nu, Impl::CovShiftBackward(U[mu], mu, Impl::CovShiftIdentityBackward(tmp, nu))), mu);

    // C4+
    out += Impl::ShiftStaple(Impl::CovShiftForward(U[nu], nu, Impl::CovShiftBackward(U[mu], mu, Impl::CovShiftIdentityBackward(U[nu], nu))), mu) * lambda;

    // insertion in lower staple
    // C1-
    out -= Impl::ShiftStaple(lambda, mu) * Impl::ShiftStaple(Impl::CovShiftBackward(U[nu], nu, Impl::CovShiftBackward(U[mu], mu, U[nu])), mu);

    // C2-
    tmp = adj(lambda) * U[nu];
    out -= Impl::ShiftStaple(Impl::CovShiftBackward(tmp, nu, Impl::CovShiftBackward(U[mu], mu, U[nu])), mu);

    // C3-
    tmp = lambda * U[nu];
    out -= Impl::ShiftStaple(Impl::CovShiftBackward(U[nu], nu, Impl::CovShiftBackward(U[mu], mu, tmp)), mu);

    // C4-
    out -= Impl::ShiftStaple(Impl::CovShiftBackward(U[nu], nu, Impl::CovShiftBackward(U[mu], mu, U[nu])), mu) * lambda;

    return out;
  }

private:
  // here fixing the 4 dimensions, make it more general?

  RealD csw_r;                                               // Clover coefficient - spatial
  RealD csw_t;                                               // Clover coefficient - temporal
  RealD diag_mass;                                           // Mass term
  CloverFieldType CloverTerm, CloverTermInv;                 // Clover term
  CloverFieldType CloverTermEven, CloverTermOdd;             // Clover term EO
  CloverFieldType CloverTermInvEven, CloverTermInvOdd;       // Clover term Inv EO
  CloverFieldType CloverTermDagEven, CloverTermDagOdd;       // Clover term Dag EO
  CloverFieldType CloverTermInvDagEven, CloverTermInvDagOdd; // Clover term Inv Dag EO

  // eventually these can be compressed into 6x6 blocks instead of the 12x12
  // using the DeGrand-Rossi basis for the gamma matrices
  CloverFieldType fillCloverYZ(const GaugeLinkField &F)
  {
    CloverFieldType T(F._grid);
    T = zero;
    PARALLEL_FOR_LOOP
    for (int i = 0; i < CloverTerm._grid->oSites(); i++)
    {
      T._odata[i]()(0, 1) = timesMinusI(F._odata[i]()());
      T._odata[i]()(1, 0) = timesMinusI(F._odata[i]()());
      T._odata[i]()(2, 3) = timesMinusI(F._odata[i]()());
      T._odata[i]()(3, 2) = timesMinusI(F._odata[i]()());
    }

    return T;
  }

  CloverFieldType fillCloverXZ(const GaugeLinkField &F)
  {
    CloverFieldType T(F._grid);
    T = zero;
    PARALLEL_FOR_LOOP
    for (int i = 0; i < CloverTerm._grid->oSites(); i++)
    {
      T._odata[i]()(0, 1) = -F._odata[i]()();
      T._odata[i]()(1, 0) = F._odata[i]()();
      T._odata[i]()(2, 3) = -F._odata[i]()();
      T._odata[i]()(3, 2) = F._odata[i]()();
    }

    return T;
  }

  CloverFieldType fillCloverXY(const GaugeLinkField &F)
  {
    CloverFieldType T(F._grid);
    T = zero;
    PARALLEL_FOR_LOOP
    for (int i = 0; i < CloverTerm._grid->oSites(); i++)
    {

      T._odata[i]()(0, 0) = timesMinusI(F._odata[i]()());
      T._odata[i]()(1, 1) = timesI(F._odata[i]()());
      T._odata[i]()(2, 2) = timesMinusI(F._odata[i]()());
      T._odata[i]()(3, 3) = timesI(F._odata[i]()());
    }

    return T;
  }

  CloverFieldType fillCloverXT(const GaugeLinkField &F)
  {
    CloverFieldType T(F._grid);
    T = zero;
    PARALLEL_FOR_LOOP
    for (int i = 0; i < CloverTerm._grid->oSites(); i++)
    {
      T._odata[i]()(0, 1) = timesI(F._odata[i]()());
      T._odata[i]()(1, 0) = timesI(F._odata[i]()());
      T._odata[i]()(2, 3) = timesMinusI(F._odata[i]()());
      T._odata[i]()(3, 2) = timesMinusI(F._odata[i]()());
    }

    return T;
  }

  CloverFieldType fillCloverYT(const GaugeLinkField &F)
  {
    CloverFieldType T(F._grid);
    T = zero;
    PARALLEL_FOR_LOOP
    for (int i = 0; i < CloverTerm._grid->oSites(); i++)
    {
      T._odata[i]()(0, 1) = -(F._odata[i]()());
      T._odata[i]()(1, 0) = (F._odata[i]()());
      T._odata[i]()(2, 3) = (F._odata[i]()());
      T._odata[i]()(3, 2) = -(F._odata[i]()());
    }

    return T;
  }

  CloverFieldType fillCloverZT(const GaugeLinkField &F)
  {
    CloverFieldType T(F._grid);
    T = zero;
    PARALLEL_FOR_LOOP
    for (int i = 0; i < CloverTerm._grid->oSites(); i++)
    {
      T._odata[i]()(0, 0) = timesI(F._odata[i]()());
      T._odata[i]()(1, 1) = timesMinusI(F._odata[i]()());
      T._odata[i]()(2, 2) = timesMinusI(F._odata[i]()());
      T._odata[i]()(3, 3) = timesI(F._odata[i]()());
    }

    return T;
  }
};
}
}

#endif // GRID_QCD_WILSON_CLOVER_FERMION_H
