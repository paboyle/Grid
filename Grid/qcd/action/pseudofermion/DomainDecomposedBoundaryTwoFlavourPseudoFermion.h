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
class DomainDecomposedBoundaryTwoFlavourPseudoFermion : public Action<typename Impl::GaugeField> {
public:
  INHERIT_IMPL_TYPES(Impl);

private:
  SchurFactoredFermionOperator<Impl> & DenOp;// the basic operator

  OperatorFunction<FermionField> &DerivativeSolver;
  OperatorFunction<FermionField> &ActionSolver;

  FermionField Phi; // the pseudo fermion field for this trajectory

  RealD refresh_action;
public:
  DomainDecomposedBoundaryTwoFlavourPseudoFermion(SchurFactoredFermionOperator<Impl>  &_DenOp,
				    OperatorFunction<FermionField> & DS,
				    OperatorFunction<FermionField> & AS
				    ) : DenOp(_DenOp),
					DerivativeSolver(DS), ActionSolver(AS),
					Phi(_DenOp.FermOp.FermionGrid()) {};
      
  virtual std::string action_name(){return "DomainDecomposedBoundaryTwoFlavourPseudoFermion";}

 
  virtual std::string LogParameters(){
    std::stringstream sstream;
    return sstream.str();
  }  
  
  virtual void refresh(const GaugeField &U, GridSerialRNG& sRNG, GridParallelRNG& pRNG)
  {
    // P(phi) = e^{- phi^dag Rdag^-1 R^-1 phi}
    //
    // DenOp == R
    //
    // Take phi = R eta  ; eta = R^-1 Phi
    //
    // P(eta) = e^{- eta^dag eta}
    //
    // e^{x^2/2 sig^2} => sig^2 = 0.5.
    // 
    // So eta should be of width sig = 1/sqrt(2) and must multiply by 0.707....
    //
    RealD scale = std::sqrt(0.5);

    DenOp.ImportGauge(U);

    FermionField eta(DenOp.FermOp.FermionGrid());

    gaussian(pRNG,eta);    eta=eta*scale;
    
    DenOp.ProjectBoundaryBar(eta);
    DenOp.R(eta,Phi);

    refresh_action = norm2(eta);
  };

  //////////////////////////////////////////////////////
  // S = phi^dag Rdag^-1 R^-1 phi
  //////////////////////////////////////////////////////
  virtual RealD S(const GaugeField &U) {

    DenOp.ImportGauge(U);

    FermionField X(DenOp.FermOp.FermionGrid());

    DenOp.RInv(Phi,X);

    RealD action = norm2(X);

    return action;
  };

  virtual void deriv(const GaugeField &U,GaugeField & dSdU)
  {
    DenOp.ImportGauge(U);

    GridBase *fgrid = DenOp.FermOp.FermionGrid();
    GridBase *ugrid = DenOp.FermOp.GaugeGrid();

    FermionField  X(fgrid);
    FermionField  Y(fgrid);
    FermionField  tmp(fgrid);

    GaugeField   force(ugrid);	

    FermionField DiDdb_Phi(fgrid);      // Vector C in my notes
    FermionField DidRinv_Phi(fgrid);    // Vector D in my notes
    FermionField DdbdDidRinv_Phi(fgrid);
    
    FermionField Rinv_Phi(fgrid);
    FermionField RinvDagRinv_Phi(fgrid);

    // R^-1 term
    DenOp.dBoundaryBar(Phi,tmp);
    DenOp.Dinverse(tmp,DiDdb_Phi);            // Vector C
    Rinv_Phi = Phi - DiDdb_Phi;
    DenOp.ProjectBoundaryBar(Rinv_Phi); 
 
    // R^-dagger R^-1 term
    DenOp.DinverseDag(Rinv_Phi,DidRinv_Phi); // Vector D
    DenOp.dBoundaryBarDag(DidRinv_Phi,DdbdDidRinv_Phi);
    RinvDagRinv_Phi = Rinv_Phi - DdbdDidRinv_Phi;
    DenOp.ProjectBoundaryBar(RinvDagRinv_Phi);

    X = DiDdb_Phi;
    Y = DidRinv_Phi;
    DenOp.FermOp.MDeriv(force,Y,X,DaggerNo);    dSdU=force;
    DenOp.FermOp.MDeriv(force,X,Y,DaggerYes);   dSdU=dSdU+force;

    dSdU *= -1.0;

  };
};

NAMESPACE_END(Grid);

