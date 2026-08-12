/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./lib/qcd/action/pseudofermion/TwoFlavourRatio4DPseudoFermion.h

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
// Two flavour ratio with FOUR dimensional pseudofermion, UNpreconditioned
// (full grid) solves.
//
// Companion to TwoFlavourRatioEO4DPseudoFermion.h but with the solver
// plumbing exposed as LinearFunction<FermionField> objects that already
// know their operator -- the natural interface for the non-Hermitian
// multigrid GCR stack (PVdagM), which solves M and Mdag DIRECTLY rather
// than through SchurRedBlack normal equations.
//
// Why: with 5D pseudofermions the squared-operator formulation hands
// normal-equation solvers (MdagM)^-1 phi AND Mdag^-1 phi from ONE Krylov
// space; a direct solver must solve twice, halving its per-solve gain.
// The 4D pseudofermion action needs one M^-1 and one M^-dag solve per
// force evaluation FOR BOTH solver families, so the direct-solver gain
// carries through undiluted.  In addition phi4 is Ls-agnostic, so the
// force can be evaluated with a reduced-Ls operator pair while the
// accept/reject uses full Ls (inexact force, exact action).
//
// Solver slots (all full-grid 5D LinearFunctions, solution overwritten,
// zero guess imposed internally):
//   DerivMinvSolver     : x = M^-1 b        (DenOp)
//   DerivMdagInvSolver  : x = M^-dag b      (DenOp).  For G5R5-hermitian
//                         actions this may be implemented by the caller as
//                         G5R5 . DerivMinvSolver . G5R5 -- no adjoint
//                         multigrid needed.
//   ActionMinvSolver    : x = M^-1 b        (DenOp, accept/reject tolerance)
//   HeatbathVinvSolver  : x = V^-1 b        (NumOp)
//
// 4D <-> 5D wall maps: the action is S = | P (M^-1 V) Pdag phi4 |^2 where
// (P,Pdag) MUST be a mutually adjoint pair for S and deriv to be
// consistent.  Two candidate conventions, selected by solution_walls:
//   true  : P = P_- psi(0) + P_+ psi(Ls-1)   (solution walls, matches
//           ExportPhysicalFermionSolution) and Pdag its literal adjoint.
//   false : P = P_+ psi(0) + P_- psi(Ls-1)   (source walls, Pdag matches
//           ImportUnphysicalFermion).
// The heatbath is exact iff [P M^-1 V Pdag][P V^-1 M Pdag] = 1 (the 4D
// effective-operator composition identity); which convention satisfies it
// is settled numerically by the refresh test  S == 0.5*|eta4|^2  exactly.
///////////////////////////////////////////////////////////////////////////////
template<class Impl>
class TwoFlavourRatio4DPseudoFermionAction : public Action<typename Impl::GaugeField> {
public:
  INHERIT_IMPL_TYPES(Impl);

private:
  typedef FermionOperator<Impl> FermOp;
  FermionOperator<Impl> & NumOp;// the basic operator (V)
  FermionOperator<Impl> & DenOp;// the basic operator (M)

  LinearFunction<FermionField> &DerivMinvSolver;
  LinearFunction<FermionField> &DerivMdagInvSolver;
  LinearFunction<FermionField> &ActionMinvSolver;
  LinearFunction<FermionField> &HeatbathVinvSolver;

  FermionField phi4; // the pseudo fermion field for this trajectory

  int solution_walls; // wall convention for the (P,Pdag) pair; see header

  ////////////////////////////////////////////////////////////////////
  // The mutually adjoint 4D <-> 5D pair.
  //   Wall4D    : q4  = P psi5      (extract)
  //   Wall4DAdj : psi5 = Pdag q4    (insert; literal adjoint of Wall4D)
  ////////////////////////////////////////////////////////////////////
  void Wall4D(const FermionField &psi5, FermionField &q4)
  {
    int Ls = NumOp.FermionGrid()->_fdimensions[0];
    FermionField tmp(NumOp.FermionGrid());
    if ( solution_walls ) {
      // q4 = P_- psi(0) + P_+ psi(Ls-1)
      axpby_ssp_pminus(tmp, 0., psi5, 1., psi5, 0, 0);
      axpby_ssp_pplus (tmp, 1., tmp , 1., psi5, 0, Ls-1);
    } else {
      // q4 = P_+ psi(0) + P_- psi(Ls-1)
      axpby_ssp_pplus (tmp, 0., psi5, 1., psi5, 0, 0);
      axpby_ssp_pminus(tmp, 1., tmp , 1., psi5, 0, Ls-1);
    }
    ExtractSlice(q4, tmp, 0, 0);
  }
  void Wall4DAdj(const FermionField &q4, FermionField &psi5)
  {
    int Ls = NumOp.FermionGrid()->_fdimensions[0];
    FermionField tmp(NumOp.FermionGrid());
    tmp = Zero();
    InsertSlice(q4, tmp, 0   , 0);
    InsertSlice(q4, tmp, Ls-1, 0);
    if ( solution_walls ) {
      // psi(0) = P_- q4 ; psi(Ls-1) = P_+ q4
      axpby_ssp_pminus(tmp, 0., tmp, 1., tmp, 0   , 0);
      axpby_ssp_pplus (tmp, 0., tmp, 1., tmp, Ls-1, Ls-1);
    } else {
      // psi(0) = P_+ q4 ; psi(Ls-1) = P_- q4
      axpby_ssp_pplus (tmp, 0., tmp, 1., tmp, 0   , 0);
      axpby_ssp_pminus(tmp, 0., tmp, 1., tmp, Ls-1, Ls-1);
    }
    psi5 = tmp;
  }

public:
  TwoFlavourRatio4DPseudoFermionAction(FermionOperator<Impl>  &_NumOp,
				       FermionOperator<Impl>  &_DenOp,
				       LinearFunction<FermionField> & DMS,
				       LinearFunction<FermionField> & DMDS,
				       LinearFunction<FermionField> & AMS,
				       LinearFunction<FermionField> & HVS,
				       int _solution_walls = 1
				       ) : NumOp(_NumOp),
					   DenOp(_DenOp),
					   DerivMinvSolver(DMS),
					   DerivMdagInvSolver(DMDS),
					   ActionMinvSolver(AMS),
					   HeatbathVinvSolver(HVS),
					   phi4(_NumOp.GaugeGrid()),
					   solution_walls(_solution_walls)
  {};

  virtual std::string action_name(){return "TwoFlavourRatio4DPseudoFermionAction";}

  virtual std::string LogParameters(){
    std::stringstream sstream;
    sstream << GridLogMessage << "["<<action_name()<<"] solution_walls " << solution_walls << std::endl;
    return sstream.str();
  }

  virtual void refresh(const GaugeField &U, GridSerialRNG &sRNG, GridParallelRNG& pRNG) {

    // P(phi4) = e^{- phi4^dag Beff^dag Beff phi4}   ;  Beff = P M^-1 V Pdag
    //
    // NumOp == V
    // DenOp == M
    //
    // Take phi4 = P V^-1 M Pdag eta4   ( = Beff^-1 eta4 by the composition
    // identity; verified numerically by S == 0.5 |eta4|^2 after refresh )
    //
    // P(eta) = e^{- eta^dag eta} ;  e^{-x^2/2 sig^2} => sig^2 = 0.5
    // so eta enters with width 1/sqrt(2).
    //
    RealD scale = std::sqrt(0.5);

    FermionField eta4(NumOp.GaugeGrid());
    FermionField eta5(NumOp.FermionGrid());
    FermionField tmp (NumOp.FermionGrid());
    FermionField phi5(NumOp.FermionGrid());

    gaussian(pRNG,eta4);

    NumOp.ImportGauge(U);
    DenOp.ImportGauge(U);

    Wall4DAdj(eta4,eta5);            // eta5 = Pdag eta4
    DenOp.M(eta5,tmp);               // tmp  = M eta5
    phi5 = Zero();
    HeatbathVinvSolver(tmp,phi5);    // phi5 = V^-1 M eta5
    Wall4D(phi5,phi4);               // phi4 = P phi5
    phi4 = phi4*scale;

    std::cout << GridLogMessage << "4d pf (non-EO) refresh "<< norm2(phi4)<<"\n";
  };

  //////////////////////////////////////////////////////
  // S = phi4^dag (Pdag^dag V^dag M^-dag P^dag) (P M^-1 V Pdag) phi4
  //   = | P M^-1 V Pdag phi4 |^2
  //////////////////////////////////////////////////////
  virtual RealD S(const GaugeField &U) {

    NumOp.ImportGauge(U);
    DenOp.ImportGauge(U);

    FermionField Y4  (NumOp.GaugeGrid());
    FermionField phi5(NumOp.FermionGrid());
    FermionField X   (NumOp.FermionGrid());
    FermionField Y   (NumOp.FermionGrid());

    Wall4DAdj(phi4,phi5);            // phi5 = Pdag phi4
    NumOp.M(phi5,X);                 // X    = V phi5
    Y = Zero();
    ActionMinvSolver(X,Y);           // Y    = M^-1 V phi5
    Wall4D(Y,Y4);                    // Y4   = P Y

    RealD action = norm2(Y4);

    return action;
  };

  //////////////////////////////////////////////////////
  // dS/du = 2 Re [ (M^-dag Pdag w4)^dag dV Pdag phi4 ]
  //       - 2 Re [ (M^-dag Pdag w4)^dag dM (M^-1 V Pdag phi4) ]
  // with w4 = P M^-1 V Pdag phi4.
  // Two first-power solves: one M^-1, one M^-dag.
  //////////////////////////////////////////////////////
  virtual void deriv(const GaugeField &U,GaugeField & dSdU) {

    NumOp.ImportGauge(U);
    DenOp.ImportGauge(U);

    FermionField  phi5           (NumOp.FermionGrid());
    FermionField  Vphi           (NumOp.FermionGrid());
    FermionField  MinvVphi       (NumOp.FermionGrid());
    FermionField  w4             (NumOp.GaugeGrid());
    FermionField  Y              (NumOp.FermionGrid());
    FermionField  MdagInvPdagW   (NumOp.FermionGrid());

    GaugeField   force(NumOp.GaugeGrid());

    Wall4DAdj(phi4,phi5);              // phi5     = Pdag phi4
    NumOp.M(phi5,Vphi);                // Vphi     = V phi5
    MinvVphi = Zero();
    DerivMinvSolver(Vphi,MinvVphi);    // MinvVphi = M^-1 V phi5
    std::cout << GridLogMessage << "4d pf (non-EO) deriv solve "<< norm2(MinvVphi)<<"\n";

    // Project onto the physical 4D subspace and back: Y = Pdag P MinvVphi.
    // Pdag here MUST be the literal adjoint of the P used in S, else the
    // force is inconsistent with the action.
    Wall4D(MinvVphi,w4);               // w4 = P MinvVphi
    Wall4DAdj(w4,Y);                   // Y  = Pdag w4

    MdagInvPdagW = Zero();
    DerivMdagInvSolver(Y,MdagInvPdagW);   // = M^-dag Pdag w4  (adjoint solve)
    std::cout << GridLogMessage << "4d pf (non-EO) deriv solve dag "<< norm2(MdagInvPdagW)<<"\n";

    // phi^dag (Pdag' Vdag Mdag^-1 P') (dV) Pdag phi   + h.c.
    NumOp.MDeriv(force, MdagInvPdagW, phi5, DaggerNo );  dSdU=force;
    NumOp.MDeriv(force, phi5, MdagInvPdagW, DaggerYes);  dSdU=dSdU+force;

    // - phi^dag ( ... Mdag^-1 ) dM ( M^-1 V ... ) phi  + h.c.
    DenOp.MDeriv(force, MdagInvPdagW, MinvVphi, DaggerNo );  dSdU=dSdU-force;
    DenOp.MDeriv(force, MinvVphi, MdagInvPdagW, DaggerYes);  dSdU=dSdU-force;

    dSdU *= -1.0;
  };
};

NAMESPACE_END(Grid);
