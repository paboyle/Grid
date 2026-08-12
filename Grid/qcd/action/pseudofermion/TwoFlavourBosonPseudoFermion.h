/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./lib/qcd/action/pseudofermion/TwoFlavourBosonPseudoFermion.h

    Copyright (C) 2026

Author: Peter Boyle <pboyle@bnl.gov>

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
#pragma once

NAMESPACE_BEGIN(Grid);

///////////////////////////////////////////////////////////////////////////////
// Two flavour BOSON (wrong-sign) pseudofermion for any FermionOperator B:
//
//    S2 = chi^dag Bdag B chi = |B chi|^2
//
//    integral  ==>  det( Bdag B )^-1 = |det B|^-2
//
// A compensator monomial: supplies an INVERSE determinant with NO solve in
// the force or the action -- both are matrix multiplies.  The only solve is
// the heatbath  chi = B^-1 eta,  once per trajectory (for B = the
// Pauli-Villars operator this is a mass-one solve, trivially cheap).
//
// Primary use: two instances with B = PV cancel the |det PV|^2 excess of
// TwoFlavourPVdagMPseudoFermionAction down to the DWF quotient
// |det M|^2/|det PV|^2 (two unsquared instances rather than one squared
// kernel: first powers of PV in the force, milder).  Being generic in B it
// also serves Hasenbusch-chain compensation at intermediate masses, or any
// future inverse-det bookkeeping.  (Sibling of the domain-decomposed boson
// in DomainDecomposedBoundaryTwoFlavourBosonPseudoFermion.h, without the
// boundary machinery.)
//
// Heatbath exact by construction: S2 after refresh = |B B^-1 eta|^2 = |eta|^2.
///////////////////////////////////////////////////////////////////////////////
template<class Impl>
class TwoFlavourBosonPseudoFermionAction : public Action<typename Impl::GaugeField> {
public:
  INHERIT_IMPL_TYPES(Impl);

private:
  FermionOperator<Impl> & BOp;  // the operator whose |det|^-2 is supplied

  LinearFunction<FermionField> &HeatbathSolver; // b -> B^-1 b (heatbath only)

  FermionField Chi; // the pseudo fermion field for this trajectory

public:
  TwoFlavourBosonPseudoFermionAction(FermionOperator<Impl>  &_BOp,
				     LinearFunction<FermionField> & HS
				     ) : BOp(_BOp),
					 HeatbathSolver(HS),
					 Chi(_BOp.FermionGrid())
  {};

  virtual std::string action_name(){return "TwoFlavourBosonPseudoFermionAction";}

  virtual std::string LogParameters(){
    std::stringstream sstream;
    sstream << GridLogMessage << "["<<action_name()<<"] has no parameters" << std::endl;
    return sstream.str();
  }

  virtual void refresh(const GaugeField &U, GridSerialRNG &sRNG, GridParallelRNG& pRNG) {
    // P(chi) = e^{- chi^dag BdagB chi} ; chi = B^-1 eta ; P(eta) = e^{-eta^dag eta}
    // e^{-x^2/2 sig^2} => sig^2 = 0.5 ; eta enters with width 1/sqrt(2).
    RealD scale = std::sqrt(0.5);
    FermionField eta(BOp.FermionGrid());
    gaussian(pRNG,eta);
    eta = eta * scale;
    refresh(U,eta);
  }

  // Deterministic-noise variant (test hook):
  // after this,  S(U) == norm2(eta)  exactly (to solver tolerance).
  void refresh(const GaugeField &U, const FermionField &eta) {
    BOp.ImportGauge(U);
    Chi = Zero();
    HeatbathSolver(eta,Chi);       // Chi = B^-1 eta : the ONLY solve
    std::cout << GridLogMessage << action_name() << " refresh |Chi|^2 = "<< norm2(Chi)<<std::endl;
  }

  //////////////////////////////////////////////////////
  // S2 = |B chi|^2  -- matrix multiply only
  //////////////////////////////////////////////////////
  virtual RealD S(const GaugeField &U) {
    BOp.ImportGauge(U);

    FermionField w(BOp.FermionGrid());
    BOp.M(Chi,w);                  // w = B chi
    RealD action = norm2(w);
    return action;
  }

  //////////////////////////////////////////////////////
  // dS2 = chi^dag dBdag w + w^dag dB chi ,   w = B chi
  // NO solves.
  //////////////////////////////////////////////////////
  virtual void deriv(const GaugeField &U,GaugeField & dSdU) {
    BOp.ImportGauge(U);

    FermionField w(BOp.FermionGrid());
    GaugeField force(BOp.GaugeGrid());

    BOp.M(Chi,w);                  // w = B chi

    BOp.MDeriv(force, Chi, w, DaggerYes); dSdU =      force;
    BOp.MDeriv(force, w, Chi, DaggerNo ); dSdU = dSdU+force;

    dSdU *= -1.0;  // Grid action sign convention (cf TwoFlavourRatio.h)
  };
};

NAMESPACE_END(Grid);
