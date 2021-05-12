/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/action/pseudofermion/DomainDecomposedTwoFlavourBoundary.h

    Copyright (C) 2021

Author: Peter Boyle <paboyle@ph.ed.ac.uk>

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

///////////////////////////////////////
// Two flavour ratio
///////////////////////////////////////
template<class Impl>
class DomainBoundaryPseudoFermionAction : public Action<typename Impl::GaugeField> {
public:
  INHERIT_IMPL_TYPES(Impl);

private:
  SchurFactoredFermionOperator<Impl> & NumOp;// the basic operator
  SchurFactoredFermionOperator<Impl> & DenOp;// the basic operator

  OperatorFunction<FermionField> &DerivativeSolver;
  OperatorFunction<FermionField> &ActionSolver;

  FermionField Phi; // the pseudo fermion field for this trajectory

public:
  DomainBoundaryPseudoFermionAction(SchurFactoredFermionOperator<Impl>  &_NumOp, 
				    SchurFactoredFermionOperator<Impl>  &_DenOp,
				    OperatorFunction<FermionField> & DS,
				    OperatorFunction<FermionField> & AS
				    ) : NumOp(_NumOp), DenOp(_DenOp),
					DerivativeSolver(DS), ActionSolver(AS),
					Phi(_NumOp.FermionGrid()) {};
      
  virtual std::string action_name(){return "DomainBoundaryPseudoFermionRatioAction";}

 
  virtual std::string LogParameters(){
    std::stringstream sstream;
    return sstream.str();
  }  
  
  virtual void refresh(const GaugeField &U, GridParallelRNG& pRNG)
  {
    // P(phi) = e^{- phi^dag P^dag Rdag^-1 R^-1 P phi}
    //
    // NumOp == P
    // DenOp == R
    //
    // Take phi = P^{-1} R eta  ; eta = R^-1 P Phi
    //
    // P(eta) = e^{- eta^dag eta}
    //
    // e^{x^2/2 sig^2} => sig^2 = 0.5.
    // 
    // So eta should be of width sig = 1/sqrt(2) and must multiply by 0.707....
    //
    RealD scale = std::sqrt(0.5);

    NumOp.ImportGauge(U);
    DenOp.ImportGauge(U);

    FermionField eta(NumOp.FermOp.FermionGrid());
    FermionField tmp(NumOp.FermOp.FermionGrid());

    gaussian(pRNG,eta);

    NumOp.ProjectBoundaryBar(eta);

    NumOp.R(eta,tmp);
    DenOp.RInv(tmp,Phi);

    Phi=Phi*scale;
	
    NumOp.ProjectBoundaryBar(Phi);
  };

  //////////////////////////////////////////////////////
  // S = phi^dag V (Mdag M)^-1 Vdag phi
  //////////////////////////////////////////////////////
  virtual RealD S(const GaugeField &U) {

    NumOp.ImportGauge(U);
    DenOp.ImportGauge(U);

    FermionField X(NumOp.FermOp.FermionGrid());
    FermionField Y(NumOp.FermOp.FermionGrid());

    NumOp.R(Phi,Y);
    DenOp.RInv(Y,X);

    RealD action = norm2(X);

    return action;
  };

  //////////////////////////////////////////////////////
  // dS/du = phi^dag dV (Mdag M)^-1 V^dag  phi
  //       - phi^dag V (Mdag M)^-1 [ Mdag dM + dMdag M ]  (Mdag M)^-1 V^dag  phi
  //       + phi^dag V (Mdag M)^-1 dV^dag  phi
  //////////////////////////////////////////////////////
  virtual void deriv(const GaugeField &U,GaugeField & dSdU) {

    NumOp.ImportGauge(U);
    DenOp.ImportGauge(U);

    GridBase *fgrid = NumOp.FermOp.FermionGrid();
    GridBase *ugrid = NumOp.FermOp.GaugeGrid();

    FermionField  X(fgrid);
    FermionField  Y(fgrid);
    FermionField  tmp(fgrid);

    GaugeField   force(ugrid);	

    FermionField DobiDdbPhi(fgrid);      // Vector A in my notes
    FermionField DoiDdDobiDdbPhi(fgrid); // Vector B in my notes
    FermionField DiDdbP_Phi(fgrid);      // Vector C in my notes
    FermionField DidRinvP_Phi(fgrid);    // Vector D in my notes
    FermionField DoidRinvDagRinvP_Phi(fgrid);    // Vector E in my notes
    FermionField DobidDddDoidRinvDagRinvP_Phi(fgrid);    // Vector F in my notes
    
    FermionField P_Phi(fgrid);
    FermionField RinvP_Phi(fgrid);
    FermionField RinvDagRinvP_Phi(fgrid);

    // P term
    NumOp.dBoundaryBar(Phi,tmp);
    NumOp.dOmegaBarInv(tmp,DobiDdbPhi);        // Vector A
    NumOp.dBoundary(DobiDdbPhi,tmp);
    NumOp.dOmegaInv(tmp,DoiDdDobiDdbPhi);      // Vector B
    P_Phi  = Phi - DoiDdDobiDdbPhi;
    NumOp.ProjectBoundaryBar(P_Phi);

    // R^-1 P term
    DenOp.dBoundaryBar(P_Phi,tmp);
    DenOp.Dinverse(tmp,DiDdbP_Phi);            // Vector C
    RinvP_Phi = P_Phi - DiDdbP_Phi;
    DenOp.ProjectBoundaryBar(RinvP_Phi);

    // R^-dagger R^-1 P term
    DenOp.DinverseDag(RinvP_Phi,DidRinvP_Phi); // Vector D
    RinvDagRinvP_Phi = RinvP_Phi - DidRinvP_Phi;
    DenOp.ProjectBoundaryBar(RinvDagRinvP_Phi);

    // P^dag R^-dagger R^-1 P term
    NumOp.dOmegaDagInv(RinvDagRinvP_Phi,DoidRinvDagRinvP_Phi); // Vector E
    NumOp.dBoundaryDag(DoidRinvDagRinvP_Phi,tmp);
    NumOp.dOmegaBarDagInv(tmp,DobidDddDoidRinvDagRinvP_Phi);   // Vector F


    dSdU=Zero();

    // phi^dag V (Mdag M)^-1 dV^dag  phi
    X = DobiDdbPhi;
    Y = DobidDddDoidRinvDagRinvP_Phi;
    NumOp.DirichletFermOp.MDeriv(force,X,Y,DaggerYes); dSdU=dsdU+force;
    NumOp.DirichletFermOp.MDeriv(force,Y,X,DaggerNo);  dSdU=dSdU+force;

    X = DoiDdDobiDdbPhi;
    Y = DoidRinvDagRinvP_Phi;
    NumOp.DirichletFermOp.MDeriv(force,X,Y,DaggerYes); dSdU=dsdU+force;
    NumOp.DirichletFermOp.MDeriv(force,Y,X,DaggerNo);  dSdU=dSdU+force;

    X = DiDdbP_Phi;
    Y = DidRinvP_Phi;
    NumOp.DirichletFermOp.MDeriv(force,X,Y,DaggerYes); dSdU=dsdU+force;
    NumOp.DirichletFermOp.MDeriv(force,Y,X,DaggerNo);  dSdU=dSdU+force;

    dSdU *= -1.0;
    //dSdU = - Ta(dSdU);

  };
};

NAMESPACE_END(Grid);

#endif
