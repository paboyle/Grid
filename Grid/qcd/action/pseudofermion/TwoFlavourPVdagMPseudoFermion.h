/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./lib/qcd/action/pseudofermion/TwoFlavourPVdagMPseudoFermion.h

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
// Two flavour pseudofermion on the COMPOSITE operator  F = PVdag M :
//
//    S1 = phi^dag (Fdag F)^-1 phi
//
//    integral  ==>  det( Mdag PV PVdag M ) = |det M|^2 |det PV|^2
//
// i.e. the target two-flavour |det M|^2 TIMES an excess |det PV|^2, to be
// cancelled by compensator monomials (two TwoFlavourBosonPseudoFermionAction
// instances on PV, each contributing |det PV|^-2, net |det PV|^-4; together
// with this action's |det PV|^2 the ensemble carries |det M|^2/|det PV|^2 --
// the standard DWF quotient).
//
// Why this shape: F = PVdag M is exactly the operator the non-Hermitian
// multigrid coarsens, so its cycles precondition (Fdag F) natively.  The
// outer solve is CG on a Hermitian positive definite system -- the
// non-normality is quarantined inside the preconditioner.  One solve per
// force evaluation; all other force ingredients are matrix multiplies.
//
// Heatbath is exact by OPERATOR ALGEBRA (no wall/projection identity):
//    refresh:  phi = Fdag eta   ==>   S1 = eta^dag F (Fdag F)^-1 Fdag eta
//                                        = |eta|^2       (to solver tol)
//
// Hasenbusch: nothing here requires PVOp to have mass one.  Any
// (heavier,lighter) pair F(m1,m2) = D^dag(m1) D(m2) works, each rung
// coarsenable by the same machinery; compensate intermediate-mass excess
// dets with boson monomials on the heavier operator.
//
// Solver slots map  b -> (Fdag F)^-1 b  (full grid, zero guess imposed
// internally).  The class is agnostic to the implementation: plain CG on
// the normal equations for testing; sequential MG solves of Fdag and F, or
// preconditioned CG with a frozen-cycle G Gdag preconditioner in production.
///////////////////////////////////////////////////////////////////////////////
template<class Impl>
class TwoFlavourPVdagMPseudoFermionAction : public Action<typename Impl::GaugeField> {
public:
  INHERIT_IMPL_TYPES(Impl);

private:
  FermionOperator<Impl> & PVOp; // the heavier / Pauli-Villars operator
  FermionOperator<Impl> & MOp;  // the lighter operator

  LinearFunction<FermionField> &DerivSolver;  // b -> (FdagF)^-1 b, MD tolerance
  LinearFunction<FermionField> &ActionSolver; // b -> (FdagF)^-1 b, accept/reject tolerance

  FermionField Phi; // the pseudo fermion field for this trajectory

  ////////////////////////////////////////////////////////////////////
  // F = PVdag M  and  Fdag = Mdag PV
  ////////////////////////////////////////////////////////////////////
  void Fapply(const FermionField &in, FermionField &out) {
    FermionField tmp(MOp.FermionGrid());
    MOp.M(in,tmp);
    PVOp.Mdag(tmp,out);
  }
  void FdagApply(const FermionField &in, FermionField &out) {
    FermionField tmp(MOp.FermionGrid());
    PVOp.M(in,tmp);
    MOp.Mdag(tmp,out);
  }

public:
  TwoFlavourPVdagMPseudoFermionAction(FermionOperator<Impl>  &_PVOp,
				      FermionOperator<Impl>  &_MOp,
				      LinearFunction<FermionField> & DS,
				      LinearFunction<FermionField> & AS
				      ) : PVOp(_PVOp),
					  MOp(_MOp),
					  DerivSolver(DS),
					  ActionSolver(AS),
					  Phi(_MOp.FermionGrid())
  {};

  virtual std::string action_name(){return "TwoFlavourPVdagMPseudoFermionAction";}

  virtual std::string LogParameters(){
    std::stringstream sstream;
    sstream << GridLogMessage << "["<<action_name()<<"] has no parameters" << std::endl;
    return sstream.str();
  }

  virtual void refresh(const GaugeField &U, GridSerialRNG &sRNG, GridParallelRNG& pRNG) {
    // P(phi) = e^{- phi^dag (FdagF)^-1 phi} ; phi = Fdag eta ; P(eta) = e^{-eta^dag eta}
    // e^{-x^2/2 sig^2} => sig^2 = 0.5 ; eta enters with width 1/sqrt(2).
    RealD scale = std::sqrt(0.5);
    FermionField eta(MOp.FermionGrid());
    gaussian(pRNG,eta);
    eta = eta * scale;
    refresh(U,eta);
  }

  // Deterministic-noise variant (test hook, TwoFlavourEvenOddRatio idiom):
  // after this,  S(U) == norm2(eta)  exactly (to solver tolerance).
  void refresh(const GaugeField &U, const FermionField &eta) {
    PVOp.ImportGauge(U);
    MOp.ImportGauge(U);
    FdagApply(eta,Phi);            // NO solve: heatbath is two matmuls
    std::cout << GridLogMessage << action_name() << " refresh |Phi|^2 = "<< norm2(Phi)<<std::endl;
  }

  //////////////////////////////////////////////////////
  // S1 = phi^dag (FdagF)^-1 phi
  //////////////////////////////////////////////////////
  virtual RealD S(const GaugeField &U) {
    PVOp.ImportGauge(U);
    MOp.ImportGauge(U);

    FermionField X(MOp.FermionGrid());
    X = Zero();
    ActionSolver(Phi,X);                       // X = (FdagF)^-1 phi

    RealD action = real(innerProduct(Phi,X));  // Hermitian positive kernel
    return action;
  }

  //////////////////////////////////////////////////////
  // dS1 = - X^dag [ dFdag F + Fdag dF ] X ,   X = (FdagF)^-1 phi,  Y = F X
  //
  // dF    = dPVdag M + PVdag dM   ==>  with A = M X,  B = PV Y :
  //
  // dS1 = - X^dag dMdag B - B^dag dM X - A^dag dPV Y - Y^dag dPVdag A
  //
  // ONE solve; A,Y,B by matrix multiply (Y = PVdag A reuses A).
  //////////////////////////////////////////////////////
  virtual void deriv(const GaugeField &U,GaugeField & dSdU) {
    PVOp.ImportGauge(U);
    MOp.ImportGauge(U);

    FermionField X(MOp.FermionGrid());
    FermionField Y(MOp.FermionGrid());
    FermionField A(MOp.FermionGrid());
    FermionField B(MOp.FermionGrid());

    GaugeField force(MOp.GaugeGrid());

    X = Zero();
    DerivSolver(Phi,X);            // X = (FdagF)^-1 phi
    MOp.M(X,A);                    // A = M X
    PVOp.Mdag(A,Y);                // Y = PVdag A = F X
    PVOp.M(Y,B);                   // B = PV Y

    // dS1 = -( X^dag dMdag B + B^dag dM X + A^dag dPV Y + Y^dag dPVdag A )
    MOp.MDeriv (force, X, B, DaggerYes); dSdU =      -force;
    MOp.MDeriv (force, B, X, DaggerNo ); dSdU = dSdU -force;
    PVOp.MDeriv(force, A, Y, DaggerNo ); dSdU = dSdU -force;
    PVOp.MDeriv(force, Y, A, DaggerYes); dSdU = dSdU -force;

    dSdU *= -1.0;  // Grid action sign convention (cf TwoFlavourRatio.h)
  };
};

NAMESPACE_END(Grid);
