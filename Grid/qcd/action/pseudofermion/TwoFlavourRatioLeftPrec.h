/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./lib/qcd/action/pseudofermion/TwoFlavourRatioLeftPrec.h

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
// Two flavour ratio with LEFT-PRECONDITIONED solves.
//
// Same action content as TwoFlavourRatio.h:
//
//    S = phi^dag V (Mdag M)^-1 Vdag phi     ==>   det[ Mdag M / Vdag V ]
//
// (V = NumOp the heavier / Pauli-Villars operator, M = DenOp the lighter),
// but organised around the composite
//
//    F = Vdag M
//
// which is the 2-hop-coarsenable operator the non-Hermitian multigrid
// serves.  Solving M X = b as F X = Vdag b is LEFT PRECONDITIONING by
// Vdag; the determinant/action layer is the standard quotient, and all
// novelty is confined to the solver contract.
//
// TwoFlavourRatio.h is tied to a normal-equations solver: one (MdagM)^-1
// solve, then Y = M X gives Mdag^-1 Vdag phi almost free.  The left-
// preconditioned idiom is DIFFERENT: the chain
//
//    b = Vdag phi
//    z : Fdag z = b          (adjoint  F solve)
//    Y = V z                 (= Mdag^-1 Vdag phi -- harvested from solve 1)
//    s = Vdag Y              (= Vdag V z)
//    X : F X = s             (forward F solve;  X = (MdagM)^-1 Vdag phi)
//
// yields Y BEFORE X (so S(U) needs only the adjoint solve), with Y's
// accuracy independent of the second solve.  Force terms are then the
// standard four MDeriv insertions of TwoFlavourRatio.
//
// Solver slots are LinearFunctions with the F-SOLVE contract (solution
// overwritten, zero guess imposed internally):
//    ForwardSolver(b,x)  :  F x    = b
//    AdjointSolver(b,z)  :  Fdag z = b
// implemented in production by the multigrid-GCR stack (forward cycle and
// adjoint cycle); in tests by CG on the composite normal equations.
// HeatbathSolver(b,x) : x = (Vdag V)^-1 b -- heavy operator, plain CG.
//
// Heatbath is exact by operator algebra:  phi = V (VdagV)^-1 Mdag eta
// ==>  S = | Mdag^-1 Vdag phi |^2 = |eta|^2  (to solver tolerance); the
// deterministic refresh(U,eta) hook below is the test point.
//
// Hasenbusch: nothing requires V to have mass one; any (heavier,lighter)
// pair works, F(V,M) = Vdag M coarsenable by the same machinery, rungs'
// solves are F-family (mrhs-batchable, mass-shared coarse space).
///////////////////////////////////////////////////////////////////////////////
template<class Impl>
class TwoFlavourRatioLeftPrecPseudoFermionAction : public Action<typename Impl::GaugeField> {
public:
  INHERIT_IMPL_TYPES(Impl);

private:
  FermionOperator<Impl> & NumOp;// V
  FermionOperator<Impl> & DenOp;// M

  LinearFunction<FermionField> &DerivForwardSolver;  // F x = b, MD tolerance
  LinearFunction<FermionField> &DerivAdjointSolver;  // Fdag z = b, MD tolerance
  LinearFunction<FermionField> &ActionAdjointSolver; // Fdag z = b, accept/reject tolerance
  LinearFunction<FermionField> &HeatbathSolver;      // (VdagV)^-1 b, heavy op

  FermionField Phi; // the pseudo fermion field for this trajectory

public:
  TwoFlavourRatioLeftPrecPseudoFermionAction(FermionOperator<Impl>  &_NumOp,
					     FermionOperator<Impl>  &_DenOp,
					     LinearFunction<FermionField> & DFS,
					     LinearFunction<FermionField> & DAS,
					     LinearFunction<FermionField> & AAS,
					     LinearFunction<FermionField> & HS
					     ) : NumOp(_NumOp),
						 DenOp(_DenOp),
						 DerivForwardSolver(DFS),
						 DerivAdjointSolver(DAS),
						 ActionAdjointSolver(AAS),
						 HeatbathSolver(HS),
						 Phi(_NumOp.FermionGrid())
  {};

  virtual std::string action_name(){return "TwoFlavourRatioLeftPrecPseudoFermionAction";}

  virtual std::string LogParameters(){
    std::stringstream sstream;
    sstream << GridLogMessage << "["<<action_name()<<"] has no parameters" << std::endl;
    return sstream.str();
  }

  virtual void refresh(const GaugeField &U, GridSerialRNG &sRNG, GridParallelRNG& pRNG) {
    // P(phi) = e^{- phi^dag V (MdagM)^-1 Vdag phi} ; phi = Vdag^-1 Mdag eta
    // e^{-x^2/2 sig^2} => sig^2 = 0.5 ; eta enters with width 1/sqrt(2).
    RealD scale = std::sqrt(0.5);
    FermionField eta(NumOp.FermionGrid());
    gaussian(pRNG,eta);
    eta = eta * scale;
    refresh(U,eta);
  }

  // Deterministic-noise variant (test hook):
  // after this,  S(U) == norm2(eta)  exactly (to solver tolerance).
  void refresh(const GaugeField &U, const FermionField &eta) {
    NumOp.ImportGauge(U);
    DenOp.ImportGauge(U);

    FermionField tmp(NumOp.FermionGrid());
    FermionField w  (NumOp.FermionGrid());

    DenOp.Mdag(eta,tmp);           // tmp = Mdag eta
    w = Zero();
    HeatbathSolver(tmp,w);         // w   = (VdagV)^-1 Mdag eta
    NumOp.M(w,Phi);                // Phi = V (VdagV)^-1 Mdag eta = Vdag^-1 Mdag eta
    std::cout << GridLogMessage << action_name() << " refresh |Phi|^2 = "<< norm2(Phi)<<std::endl;
  }

  //////////////////////////////////////////////////////
  // S = phi^dag V (MdagM)^-1 Vdag phi = | Mdag^-1 Vdag phi |^2
  // ONE adjoint F solve:  Y = V Fdag^-1 Vdag phi = Mdag^-1 Vdag phi
  //////////////////////////////////////////////////////
  virtual RealD S(const GaugeField &U) {
    NumOp.ImportGauge(U);
    DenOp.ImportGauge(U);

    FermionField b(NumOp.FermionGrid());
    FermionField z(NumOp.FermionGrid());
    FermionField Y(NumOp.FermionGrid());

    NumOp.Mdag(Phi,b);             // b = Vdag phi
    z = Zero();
    ActionAdjointSolver(b,z);      // Fdag z = b
    NumOp.M(z,Y);                  // Y = V z = Mdag^-1 Vdag phi

    RealD action = norm2(Y);
    return action;
  }

  //////////////////////////////////////////////////////
  // dS/du = phi^dag dV (MdagM)^-1 Vdag phi
  //       - phi^dag V (MdagM)^-1 [ Mdag dM + dMdag M ] (MdagM)^-1 Vdag phi
  //       + phi^dag V (MdagM)^-1 dVdag phi
  // Identical force insertions to TwoFlavourRatio.h; X and Y from the
  // left-preconditioned chain (Y harvested from the adjoint solve).
  //////////////////////////////////////////////////////
  virtual void deriv(const GaugeField &U,GaugeField & dSdU) {
    NumOp.ImportGauge(U);
    DenOp.ImportGauge(U);

    FermionField b(NumOp.FermionGrid());
    FermionField z(NumOp.FermionGrid());
    FermionField Y(NumOp.FermionGrid());
    FermionField s(NumOp.FermionGrid());
    FermionField X(NumOp.FermionGrid());

    GaugeField force(NumOp.GaugeGrid());

    NumOp.Mdag(Phi,b);             // b = Vdag phi
    z = Zero();
    DerivAdjointSolver(b,z);       // Fdag z = b
    NumOp.M(z,Y);                  // Y = V z = Mdag^-1 Vdag phi   (solve-1 harvest)
    NumOp.Mdag(Y,s);               // s = Vdag V z
    X = Zero();
    DerivForwardSolver(s,X);       // F X = s  ==>  X = (MdagM)^-1 Vdag phi

    // phi^dag V (MdagM)^-1 dVdag phi
    NumOp.MDeriv(force , X, Phi, DaggerYes);  dSdU  =      force;
    // phi^dag dV (MdagM)^-1 Vdag phi
    NumOp.MDeriv(force , Phi, X, DaggerNo );  dSdU  = dSdU+force;
    //  - phi^dag V (MdagM)^-1 Mdag dM (MdagM)^-1 Vdag phi
    //  - phi^dag V (MdagM)^-1 dMdag M (MdagM)^-1 Vdag phi
    DenOp.MDeriv(force, Y, X, DaggerNo );     dSdU  = dSdU-force;
    DenOp.MDeriv(force, X, Y, DaggerYes);     dSdU  = dSdU-force;

    dSdU *= -1.0;
  };
};

NAMESPACE_END(Grid);
