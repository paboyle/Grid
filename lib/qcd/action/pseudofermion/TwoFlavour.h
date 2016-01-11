    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/action/pseudofermion/TwoFlavour.h

    Copyright (C) 2015

Author: Peter Boyle <pabobyle@ph.ed.ac.uk>
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
#ifndef QCD_PSEUDOFERMION_TWO_FLAVOUR_H
#define QCD_PSEUDOFERMION_TWO_FLAVOUR_H

namespace Grid{
  namespace QCD{

    ////////////////////////////////////////////////////////////////////////
    // Two flavour pseudofermion action for any dop
    ////////////////////////////////////////////////////////////////////////
    template<class Impl>
    class TwoFlavourPseudoFermionAction : public Action<typename Impl::GaugeField> {
    public:
      INHERIT_IMPL_TYPES(Impl);

    private:
      
      FermionOperator<Impl> & FermOp;// the basic operator

      OperatorFunction<FermionField> &DerivativeSolver;

      OperatorFunction<FermionField> &ActionSolver;

      FermionField Phi; // the pseudo fermion field for this trajectory

    public:
      /////////////////////////////////////////////////
      // Pass in required objects.
      /////////////////////////////////////////////////
    TwoFlavourPseudoFermionAction(FermionOperator<Impl>  &Op, 
				  OperatorFunction<FermionField> & DS,
				  OperatorFunction<FermionField> & AS
				  ) : FermOp(Op), DerivativeSolver(DS), ActionSolver(AS), Phi(Op.FermionGrid()) {
      };
      
      //////////////////////////////////////////////////////////////////////////////////////
      // Push the gauge field in to the dops. Assume any BC's and smearing already applied
      //////////////////////////////////////////////////////////////////////////////////////
      virtual void refresh(const GaugeField &U, GridParallelRNG& pRNG) {

	// P(phi) = e^{- phi^dag (MdagM)^-1 phi}
	// Phi = Mdag eta 
	// P(eta) = e^{- eta^dag eta}
	//
	// e^{x^2/2 sig^2} => sig^2 = 0.5.
	// 
	// So eta should be of width sig = 1/sqrt(2).
	// and must multiply by 0.707....
	//
	// Chroma has this scale factor: two_flavor_monomial_w.h
	// IroIro: does not use this scale. It is absorbed by a change of vars
	//         in the Phi integral, and thus is only an irrelevant prefactor for the partition function.
	//
	RealD scale = std::sqrt(0.5);
	FermionField eta(FermOp.FermionGrid());

	gaussian(pRNG,eta);

	FermOp.ImportGauge(U);
	FermOp.Mdag(eta,Phi);

	Phi=Phi*scale;
	
      };

      //////////////////////////////////////////////////////
      // S = phi^dag (Mdag M)^-1 phi
      //////////////////////////////////////////////////////
      virtual RealD S(const GaugeField &U) {

	FermOp.ImportGauge(U);

	FermionField X(FermOp.FermionGrid());
	FermionField Y(FermOp.FermionGrid());
	
	MdagMLinearOperator<FermionOperator<Impl> ,FermionField> MdagMOp(FermOp);
	X=zero;
	ActionSolver(MdagMOp,Phi,X);
	MdagMOp.Op(X,Y);

	RealD action = norm2(Y);
	std::cout << GridLogMessage << "Pseudofermion action "<<action<<std::endl;
	return action;
      };

      //////////////////////////////////////////////////////
      // dS/du = - phi^dag  (Mdag M)^-1 [ Mdag dM + dMdag M ]  (Mdag M)^-1 phi
      //       = - phi^dag M^-1 dM (MdagM)^-1 phi -  phi^dag (MdagM)^-1 dMdag dM (Mdag)^-1 phi 
      //
      //       = - Ydag dM X  - Xdag dMdag Y
      //
      //////////////////////////////////////////////////////
      virtual void deriv(const GaugeField &U,GaugeField & dSdU) {

	FermOp.ImportGauge(U);

	FermionField X(FermOp.FermionGrid());
	FermionField Y(FermOp.FermionGrid());
	GaugeField   tmp(FermOp.GaugeGrid());

	MdagMLinearOperator<FermionOperator<Impl> ,FermionField> MdagMOp(FermOp);

	X=zero;
	DerivativeSolver(MdagMOp,Phi,X);
	MdagMOp.Op(X,Y);

	// Our conventions really make this UdSdU; We do not differentiate wrt Udag here.
	// So must take dSdU - adj(dSdU) and left multiply by mom to get dS/dt.

	FermOp.MDeriv(tmp , Y, X,DaggerNo );  dSdU=tmp;
	FermOp.MDeriv(tmp , X, Y,DaggerYes);  dSdU=dSdU+tmp;
	
	dSdU = Ta(dSdU);

      };

    };
    
  }
}

#endif
