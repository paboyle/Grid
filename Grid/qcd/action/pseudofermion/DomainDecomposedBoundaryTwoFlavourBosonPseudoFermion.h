/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/action/pseudofermion/DomainDecomposedTwoFlavourBoundaryBoson.h

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
template<class ImplD,class ImplF>
class DomainDecomposedBoundaryTwoFlavourBosonPseudoFermion : public Action<typename ImplD::GaugeField> {
public:
  INHERIT_IMPL_TYPES(ImplD);

private:
  SchurFactoredFermionOperator<ImplD,ImplF> & NumOp;// the basic operator
  RealD InnerStoppingCondition;
  RealD ActionStoppingCondition;
  RealD DerivativeStoppingCondition;
  FermionField Phi; // the pseudo fermion field for this trajectory
public:
  DomainDecomposedBoundaryTwoFlavourBosonPseudoFermion(SchurFactoredFermionOperator<ImplD,ImplF>  &_NumOp,RealD _DerivativeTol, RealD _ActionTol, RealD _InnerTol=1.0e-6)
    : NumOp(_NumOp), 
      DerivativeStoppingCondition(_DerivativeTol),
      ActionStoppingCondition(_ActionTol),
      InnerStoppingCondition(_InnerTol),
      Phi(_NumOp.FermionGrid()) {};

  virtual std::string action_name(){return "DomainDecomposedBoundaryTwoFlavourBosonPseudoFermion";}

  virtual std::string LogParameters(){
    std::stringstream sstream;
    return sstream.str();
  }  
  
  virtual void refresh(const GaugeField &U, GridSerialRNG& sRNG, GridParallelRNG& pRNG)
  {
    // P(phi) = e^{- phi^dag P^dag P phi}
    //
    // NumOp == P
    //
    // Take phi = P^{-1} eta  ; eta = P Phi
    //
    // P(eta) = e^{- eta^dag eta}
    //
    // e^{x^2/2 sig^2} => sig^2 = 0.5.
    // 
    // So eta should be of width sig = 1/sqrt(2) and must multiply by 0.707....
    //
    RealD scale = std::sqrt(0.5);

    NumOp.tolinner=InnerStoppingCondition;
    NumOp.tol=ActionStoppingCondition;
    NumOp.ImportGauge(U);

    FermionField eta(NumOp.FermionGrid());

    gaussian(pRNG,eta);    eta=eta*scale;
    
    NumOp.ProjectBoundaryBar(eta);
    //DumpSliceNorm("eta",eta);
    NumOp.RInv(eta,Phi);

    //DumpSliceNorm("Phi",Phi);

  };

  //////////////////////////////////////////////////////
  // S = phi^dag Pdag P phi
  //////////////////////////////////////////////////////
  virtual RealD S(const GaugeField &U) {

    NumOp.tolinner=InnerStoppingCondition;
    NumOp.tol=ActionStoppingCondition;
    NumOp.ImportGauge(U);

    FermionField Y(NumOp.FermionGrid());

    NumOp.R(Phi,Y);

    RealD action = norm2(Y);

    return action;
  };

  virtual void deriv(const GaugeField &U,GaugeField & dSdU)
  {
    NumOp.tolinner=InnerStoppingCondition;
    NumOp.tol=DerivativeStoppingCondition;
    NumOp.ImportGauge(U);

    GridBase *fgrid = NumOp.FermionGrid();
    GridBase *ugrid = NumOp.GaugeGrid();

    FermionField  X(fgrid);
    FermionField  Y(fgrid);
    FermionField  tmp(fgrid);

    GaugeField   force(ugrid);	

    FermionField DobiDdbPhi(fgrid);      // Vector A in my notes
    FermionField DoiDdDobiDdbPhi(fgrid); // Vector B in my notes
    FermionField DoidP_Phi(fgrid);    // Vector E in my notes
    FermionField DobidDddDoidP_Phi(fgrid);    // Vector F in my notes
    
    FermionField P_Phi(fgrid);
    
    // P term
    NumOp.dBoundaryBar(Phi,tmp);
    NumOp.dOmegaBarInv(tmp,DobiDdbPhi);        // Vector A
    NumOp.dBoundary(DobiDdbPhi,tmp);
    NumOp.dOmegaInv(tmp,DoiDdDobiDdbPhi);      // Vector B
    P_Phi  = Phi - DoiDdDobiDdbPhi;
    NumOp.ProjectBoundaryBar(P_Phi);
    
    // P^dag P term
    NumOp.dOmegaDagInv(P_Phi,DoidP_Phi); // Vector E
    NumOp.dBoundaryDag(DoidP_Phi,tmp);
    NumOp.dOmegaBarDagInv(tmp,DobidDddDoidP_Phi);   // Vector F
    NumOp.dBoundaryBarDag(DobidDddDoidP_Phi,tmp);

    X = DobiDdbPhi;
    Y = DobidDddDoidP_Phi;
    NumOp.DirichletFermOpD.MDeriv(force,Y,X,DaggerNo);    dSdU=force;
    NumOp.DirichletFermOpD.MDeriv(force,X,Y,DaggerYes);   dSdU=dSdU+force;

    X = DoiDdDobiDdbPhi;
    Y = DoidP_Phi;
    NumOp.DirichletFermOpD.MDeriv(force,Y,X,DaggerNo);    dSdU=dSdU+force;
    NumOp.DirichletFermOpD.MDeriv(force,X,Y,DaggerYes);   dSdU=dSdU+force;

    dSdU *= -1.0;

  };
};

NAMESPACE_END(Grid);

