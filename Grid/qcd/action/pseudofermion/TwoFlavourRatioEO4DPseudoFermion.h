/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/action/pseudofermion/TwoFlavourRatio.h

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: Peter Boyle <peterboyle@Peters-MacBook-Pro-2.local>
Author: paboyle <paboyle@ph.ed.ac.uk>

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
class TwoFlavourRatioEO4DPseudoFermionAction : public Action<typename Impl::GaugeField> {
public:
  INHERIT_IMPL_TYPES(Impl);

private:
  typedef FermionOperator<Impl> FermOp;
  FermionOperator<Impl> & NumOp;// the basic operator
  FermionOperator<Impl> & DenOp;// the basic operator

  OperatorFunction<FermionField> &DerivativeSolver;
  OperatorFunction<FermionField> &DerivativeDagSolver;
  OperatorFunction<FermionField> &ActionSolver;
  OperatorFunction<FermionField> &HeatbathSolver;

  FermionField phi4; // the pseudo fermion field for this trajectory

public:
  TwoFlavourRatioEO4DPseudoFermionAction(FermionOperator<Impl>  &_NumOp, 
					 FermionOperator<Impl>  &_DenOp, 
					 OperatorFunction<FermionField> & DS,
					 OperatorFunction<FermionField> & AS ) : 
    TwoFlavourRatioEO4DPseudoFermionAction(_NumOp,_DenOp, DS,DS,AS,AS) {};
  TwoFlavourRatioEO4DPseudoFermionAction(FermionOperator<Impl>  &_NumOp, 
					 FermionOperator<Impl>  &_DenOp, 
					 OperatorFunction<FermionField> & DS,
					 OperatorFunction<FermionField> & DDS,
					 OperatorFunction<FermionField> & AS,
					 OperatorFunction<FermionField> & HS
				       ) : NumOp(_NumOp),
					   DenOp(_DenOp),
					   DerivativeSolver(DS),
					   DerivativeDagSolver(DDS),
					   ActionSolver(AS),
					   HeatbathSolver(HS),
					   phi4(_NumOp.GaugeGrid())
  {};
      
  virtual std::string action_name(){return "TwoFlavourRatioEO4DPseudoFermionAction";}

  virtual std::string LogParameters(){
    std::stringstream sstream;
    sstream << GridLogMessage << "["<<action_name()<<"] has no parameters" << std::endl;
    return sstream.str();
  }  
      
  virtual void refresh(const GaugeField &U, GridSerialRNG &sRNG, GridParallelRNG& pRNG) {

    // P(phi) = e^{- phi^dag (V^dag M^-dag)_11  (M^-1 V)_11 phi}
    //
    // NumOp == V
    // DenOp == M
    //
    // Take phi = (V^{-1} M)_11 eta  ; eta = (M^{-1} V)_11 Phi
    //
    // P(eta) = e^{- eta^dag eta}
    //
    // e^{x^2/2 sig^2} => sig^2 = 0.5.
    // 
    // So eta should be of width sig = 1/sqrt(2) and must multiply by 0.707....
    //
    RealD scale = std::sqrt(0.5);

    FermionField eta4(NumOp.GaugeGrid());
    FermionField eta5(NumOp.FermionGrid());
    FermionField tmp(NumOp.FermionGrid());
    FermionField phi5(NumOp.FermionGrid());

    gaussian(pRNG,eta4);
    NumOp.ImportFourDimPseudoFermion(eta4,eta5);
    NumOp.ImportGauge(U);
    DenOp.ImportGauge(U);

    SchurRedBlackDiagMooeeSolve<FermionField> PrecSolve(HeatbathSolver);

    DenOp.M(eta5,tmp);               // M eta
    PrecSolve(NumOp,tmp,phi5);  // phi = V^-1 M eta
    phi5=phi5*scale;
    std::cout << GridLogMessage << "4d pf refresh "<< norm2(phi5)<<"\n";
    // Project to 4d
    NumOp.ExportFourDimPseudoFermion(phi5,phi4);
      
  };

  //////////////////////////////////////////////////////
  // S = phi^dag (V^dag M^-dag)_11  (M^-1 V)_11 phi
  //////////////////////////////////////////////////////
  virtual RealD S(const GaugeField &U) {

    NumOp.ImportGauge(U);
    DenOp.ImportGauge(U);

    FermionField Y4(NumOp.GaugeGrid());
    FermionField X(NumOp.FermionGrid());
    FermionField Y(NumOp.FermionGrid());
    FermionField phi5(NumOp.FermionGrid());
	
    MdagMLinearOperator<FermionOperator<Impl> ,FermionField> MdagMOp(DenOp);
    SchurRedBlackDiagMooeeSolve<FermionField> PrecSolve(ActionSolver);

    NumOp.ImportFourDimPseudoFermion(phi4,phi5);
    NumOp.M(phi5,X);              // X= V phi
    PrecSolve(DenOp,X,Y);    // Y= (MdagM)^-1 Mdag Vdag phi = M^-1 V phi
    NumOp.ExportFourDimPseudoFermion(Y,Y4);

    RealD action = norm2(Y4);

    return action;
  };

  //////////////////////////////////////////////////////
  // dS/du = 2 Re phi^dag (V^dag M^-dag)_11  (M^-1 d V)_11  phi
  //       - 2 Re phi^dag (dV^dag M^-dag)_11  (M^-1 dM M^-1 V)_11  phi
  //////////////////////////////////////////////////////
  virtual void deriv(const GaugeField &U,GaugeField & dSdU) {

    NumOp.ImportGauge(U);
    DenOp.ImportGauge(U);

    FermionField  X(NumOp.FermionGrid());
    FermionField  Y(NumOp.FermionGrid());
    FermionField       phi(NumOp.FermionGrid());
    FermionField      Vphi(NumOp.FermionGrid());
    FermionField  MinvVphi(NumOp.FermionGrid());
    FermionField      tmp4(NumOp.GaugeGrid());
    FermionField  MdagInvMinvVphi(NumOp.FermionGrid());

    GaugeField   force(NumOp.GaugeGrid());	

    //Y=V phi
    //X = (Mdag V phi
    //Y = (Mdag M)^-1 Mdag V phi = M^-1 V Phi
    NumOp.ImportFourDimPseudoFermion(phi4,phi);
    NumOp.M(phi,Vphi);               //  V phi
    SchurRedBlackDiagMooeeSolve<FermionField> PrecSolve(DerivativeSolver);
    PrecSolve(DenOp,Vphi,MinvVphi);// M^-1 V phi
    std::cout << GridLogMessage << "4d deriv solve "<< norm2(MinvVphi)<<"\n";

    // Projects onto the physical space and back
    NumOp.ExportFourDimPseudoFermion(MinvVphi,tmp4);
    NumOp.ImportFourDimPseudoFermion(tmp4,Y);

    SchurRedBlackDiagMooeeDagSolve<FermionField> PrecDagSolve(DerivativeDagSolver);
    // X = proj M^-dag V phi
    // Need an adjoint solve
    PrecDagSolve(DenOp,Y,MdagInvMinvVphi);
    std::cout << GridLogMessage << "4d deriv solve dag "<< norm2(MdagInvMinvVphi)<<"\n";
    
    // phi^dag (Vdag Mdag^-1) (M^-1 dV)  phi
    NumOp.MDeriv(force ,MdagInvMinvVphi , phi, DaggerNo );  dSdU=force;
  
    // phi^dag (dVdag Mdag^-1) (M^-1 V)  phi
    NumOp.MDeriv(force , phi, MdagInvMinvVphi ,DaggerYes  );  dSdU=dSdU+force;

    //    - 2 Re phi^dag (dV^dag M^-dag)_11  (M^-1 dM M^-1 V)_11  phi
    DenOp.MDeriv(force,MdagInvMinvVphi,MinvVphi,DaggerNo);   dSdU=dSdU-force;
    DenOp.MDeriv(force,MinvVphi,MdagInvMinvVphi,DaggerYes);  dSdU=dSdU-force;

    dSdU *= -1.0; 
    //dSdU = - Ta(dSdU);
    
  };
};

NAMESPACE_END(Grid);


