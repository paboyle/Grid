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
#ifndef QCD_PSEUDOFERMION_TWO_FLAVOUR_RATIO_H
#define QCD_PSEUDOFERMION_TWO_FLAVOUR_RATIO_H

namespace Grid{
  namespace QCD{

    ///////////////////////////////////////
    // Two flavour ratio
    ///////////////////////////////////////
    template<class Impl>
    class TwoFlavourRatioPseudoFermionAction : public Action<typename Impl::GaugeField> {
    public:
      INHERIT_IMPL_TYPES(Impl);

    private:
      FermionOperator<Impl> & NumOp;// the basic operator
      FermionOperator<Impl> & DenOp;// the basic operator

      OperatorFunction<FermionField> &DerivativeSolver;
      OperatorFunction<FermionField> &ActionSolver;

      FermionField Phi; // the pseudo fermion field for this trajectory

    public:
      TwoFlavourRatioPseudoFermionAction(FermionOperator<Impl>  &_NumOp, 
					 FermionOperator<Impl>  &_DenOp, 
					 OperatorFunction<FermionField> & DS,
					 OperatorFunction<FermionField> & AS
					 ) : NumOp(_NumOp), DenOp(_DenOp), DerivativeSolver(DS), ActionSolver(AS), Phi(_NumOp.FermionGrid()) {};
      
      virtual void refresh(const GaugeField &U, GridParallelRNG& pRNG) {

	// P(phi) = e^{- phi^dag V (MdagM)^-1 Vdag phi}
	//
	// NumOp == V
	// DenOp == M
	//
	// Take phi = Vdag^{-1} Mdag eta  ; eta = Mdag^{-1} Vdag Phi
	//
	// P(eta) = e^{- eta^dag eta}
	//
	// e^{x^2/2 sig^2} => sig^2 = 0.5.
	// 
	// So eta should be of width sig = 1/sqrt(2) and must multiply by 0.707....
	//
	RealD scale = std::sqrt(0.5);

	FermionField eta(NumOp.FermionGrid());
	FermionField tmp(NumOp.FermionGrid());

	gaussian(pRNG,eta);

	NumOp.ImportGauge(U);
	DenOp.ImportGauge(U);

	// Note: this hard codes normal equations type solvers; alternate implementation needed for 
	// non-herm style solvers.
	MdagMLinearOperator<FermionOperator<Impl> ,FermionField> MdagMOp(NumOp);

	DenOp.Mdag(eta,Phi);            // Mdag eta
	tmp = zero;
	ActionSolver(MdagMOp,Phi,tmp);  // (VdagV)^-1 Mdag eta = V^-1 Vdag^-1 Mdag eta
	NumOp.M(tmp,Phi);               // Vdag^-1 Mdag eta

	Phi=Phi*scale;
	
      };

      //////////////////////////////////////////////////////
      // S = phi^dag V (Mdag M)^-1 Vdag phi
      //////////////////////////////////////////////////////
      virtual RealD S(const GaugeField &U) {

	NumOp.ImportGauge(U);
	DenOp.ImportGauge(U);

	FermionField X(NumOp.FermionGrid());
	FermionField Y(NumOp.FermionGrid());
	
	MdagMLinearOperator<FermionOperator<Impl> ,FermionField> MdagMOp(DenOp);

	NumOp.Mdag(Phi,Y);              // Y= Vdag phi
	X=zero;
	ActionSolver(MdagMOp,Y,X);      // X= (MdagM)^-1 Vdag phi
	DenOp.M(X,Y);                  // Y=  Mdag^-1 Vdag phi

	RealD action = norm2(Y);

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

	MdagMLinearOperator<FermionOperator<Impl> ,FermionField> MdagMOp(DenOp);

	FermionField  X(NumOp.FermionGrid());
	FermionField  Y(NumOp.FermionGrid());

	GaugeField   force(NumOp.GaugeGrid());	


	//Y=Vdag phi
	//X = (Mdag M)^-1 V^dag phi
	//Y = (Mdag)^-1 V^dag  phi
	NumOp.Mdag(Phi,Y);              // Y= Vdag phi
	X=zero;
	DerivativeSolver(MdagMOp,Y,X);      // X= (MdagM)^-1 Vdag phi
	DenOp.M(X,Y);                  // Y=  Mdag^-1 Vdag phi

	// phi^dag V (Mdag M)^-1 dV^dag  phi
	NumOp.MDeriv(force , X, Phi, DaggerYes );  dSdU=force;
  
	// phi^dag dV (Mdag M)^-1 V^dag  phi
	NumOp.MDeriv(force , Phi, X ,DaggerNo  );  dSdU=dSdU+force;

	//    -    phi^dag V (Mdag M)^-1 Mdag dM   (Mdag M)^-1 V^dag  phi
	//    -    phi^dag V (Mdag M)^-1 dMdag M   (Mdag M)^-1 V^dag  phi
	DenOp.MDeriv(force,Y,X,DaggerNo);   dSdU=dSdU-force;
	DenOp.MDeriv(force,X,Y,DaggerYes);  dSdU=dSdU-force;

	dSdU = - Ta(dSdU);

      };
    };
  }
}
#endif
