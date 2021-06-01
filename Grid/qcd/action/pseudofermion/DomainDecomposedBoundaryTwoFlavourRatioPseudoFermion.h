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
template<class ImplD,class ImplF>
class DomainDecomposedBoundaryTwoFlavourRatioPseudoFermion : public Action<typename ImplD::GaugeField> {
public:
  INHERIT_IMPL_TYPES(ImplD);

private:
  SchurFactoredFermionOperator<ImplD,ImplF> & NumOp;// the basic operator
  SchurFactoredFermionOperator<ImplD,ImplF> & DenOp;// the basic operator

  RealD ActionStoppingCondition;
  RealD DerivativeStoppingCondition;
  
  FermionField Phi; // the pseudo fermion field for this trajectory

public:
  DomainDecomposedBoundaryTwoFlavourRatioPseudoFermion(SchurFactoredFermionOperator<ImplD,ImplF>  &_NumOp, 
						       SchurFactoredFermionOperator<ImplD,ImplF>  &_DenOp,
						       RealD _DerivativeTol, RealD _ActionTol)
    : NumOp(_NumOp), DenOp(_DenOp),
      Phi(_NumOp.PeriodicFermOpD.FermionGrid()),
      DerivativeStoppingCondition(_DerivativeTol),
      ActionStoppingCondition(_ActionTol)
  {};
      
  virtual std::string action_name(){return "DomainDecomposedBoundaryTwoFlavourRatioPseudoFermion";}
 
  virtual std::string LogParameters(){
    std::stringstream sstream;
    return sstream.str();
  }  
  
  virtual void refresh(const GaugeField &U, GridSerialRNG& sRNG, GridParallelRNG& pRNG)
  {
    NumOp.ImportGauge(U);
    DenOp.ImportGauge(U);

    FermionField eta(NumOp.PeriodicFermOpD.FermionGrid());
    FermionField tmp(NumOp.PeriodicFermOpD.FermionGrid());

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

    gaussian(pRNG,eta);    eta=eta*scale;
    
    NumOp.ProjectBoundaryBar(eta);
    DenOp.tol = ActionStoppingCondition;
    NumOp.tol = ActionStoppingCondition;
    DenOp.R(eta,tmp);
    NumOp.RInv(tmp,Phi);

  };

  //////////////////////////////////////////////////////
  // S = phi^dag Pdag Rdag^-1 R^-1 P phi
  //////////////////////////////////////////////////////
  virtual RealD S(const GaugeField &U) {

    NumOp.ImportGauge(U);
    DenOp.ImportGauge(U);

    FermionField X(NumOp.PeriodicFermOpD.FermionGrid());
    FermionField Y(NumOp.PeriodicFermOpD.FermionGrid());

    DenOp.tol = ActionStoppingCondition;
    NumOp.tol = ActionStoppingCondition;
    NumOp.R(Phi,Y);
    DenOp.RInv(Y,X);

    RealD action = norm2(X);
    //    std::cout << " DD boundary action is " <<action<<std::endl;

    return action;
  };

  virtual void deriv(const GaugeField &U,GaugeField & dSdU)
  {
    NumOp.ImportGauge(U);
    DenOp.ImportGauge(U);

    GridBase *fgrid = NumOp.PeriodicFermOpD.FermionGrid();
    GridBase *ugrid = NumOp.PeriodicFermOpD.GaugeGrid();

    FermionField  X(fgrid);
    FermionField  Y(fgrid);
    FermionField  tmp(fgrid);

    GaugeField   force(ugrid);	

    FermionField DobiDdbPhi(fgrid);      // Vector A in my notes
    FermionField DoiDdDobiDdbPhi(fgrid); // Vector B in my notes
    FermionField DiDdbP_Phi(fgrid);      // Vector C in my notes
    FermionField DidRinvP_Phi(fgrid);    // Vector D in my notes
    FermionField DdbdDidRinvP_Phi(fgrid);
    FermionField DoidRinvDagRinvP_Phi(fgrid);    // Vector E in my notes
    FermionField DobidDddDoidRinvDagRinvP_Phi(fgrid);    // Vector F in my notes
    
    FermionField P_Phi(fgrid);
    FermionField RinvP_Phi(fgrid);
    FermionField RinvDagRinvP_Phi(fgrid);
    FermionField PdagRinvDagRinvP_Phi(fgrid);

    //    RealD action = S(U);
    DenOp.tol = DerivativeStoppingCondition;
    NumOp.tol = DerivativeStoppingCondition;
    
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
    DenOp.ProjectBoundaryBar(RinvP_Phi); // Correct to here

 
    // R^-dagger R^-1 P term
    DenOp.DinverseDag(RinvP_Phi,DidRinvP_Phi); // Vector D
    DenOp.dBoundaryBarDag(DidRinvP_Phi,DdbdDidRinvP_Phi);
    RinvDagRinvP_Phi = RinvP_Phi - DdbdDidRinvP_Phi;
    DenOp.ProjectBoundaryBar(RinvDagRinvP_Phi);

    
    // P^dag R^-dagger R^-1 P term
    NumOp.dOmegaDagInv(RinvDagRinvP_Phi,DoidRinvDagRinvP_Phi); // Vector E
    NumOp.dBoundaryDag(DoidRinvDagRinvP_Phi,tmp);
    NumOp.dOmegaBarDagInv(tmp,DobidDddDoidRinvDagRinvP_Phi);   // Vector F
    NumOp.dBoundaryBarDag(DobidDddDoidRinvDagRinvP_Phi,tmp);
    PdagRinvDagRinvP_Phi = RinvDagRinvP_Phi- tmp;
    NumOp.ProjectBoundaryBar(PdagRinvDagRinvP_Phi);

    /*
    std::cout << "S eval  "<< action << std::endl;
    std::cout << "S - IP1 "<< innerProduct(Phi,PdagRinvDagRinvP_Phi) << std::endl;
    std::cout << "S - IP2 "<< norm2(RinvP_Phi) << std::endl;

    NumOp.R(Phi,tmp);
    tmp = tmp - P_Phi;
    std::cout << "diff1 "<<norm2(tmp) <<std::endl;
    
    
    DenOp.RInv(P_Phi,tmp);
    tmp = tmp - RinvP_Phi;
    std::cout << "diff2 "<<norm2(tmp) <<std::endl;

    DenOp.RDagInv(RinvP_Phi,tmp);
    tmp  = tmp - RinvDagRinvP_Phi;
    std::cout << "diff3 "<<norm2(tmp) <<std::endl;

    DenOp.RDag(RinvDagRinvP_Phi,tmp);
    tmp  = tmp - PdagRinvDagRinvP_Phi;
    std::cout << "diff4 "<<norm2(tmp) <<std::endl;
    */
    
    dSdU=Zero();

    X = DobiDdbPhi;
    Y = DobidDddDoidRinvDagRinvP_Phi;
    NumOp.DirichletFermOpD.MDeriv(force,Y,X,DaggerNo);    dSdU=dSdU+force;
    NumOp.DirichletFermOpD.MDeriv(force,X,Y,DaggerYes);   dSdU=dSdU+force;

    X = DoiDdDobiDdbPhi;
    Y = DoidRinvDagRinvP_Phi;
    NumOp.DirichletFermOpD.MDeriv(force,Y,X,DaggerNo);    dSdU=dSdU+force;
    NumOp.DirichletFermOpD.MDeriv(force,X,Y,DaggerYes);   dSdU=dSdU+force;

    X = DiDdbP_Phi;
    Y = DidRinvP_Phi;
    DenOp.PeriodicFermOpD.MDeriv(force,Y,X,DaggerNo);    dSdU=dSdU+force;
    DenOp.PeriodicFermOpD.MDeriv(force,X,Y,DaggerYes);   dSdU=dSdU+force;

    dSdU *= -1.0;

  };
};

NAMESPACE_END(Grid);

