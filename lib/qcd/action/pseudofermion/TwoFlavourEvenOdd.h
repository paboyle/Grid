    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/action/pseudofermion/TwoFlavourEvenOdd.h

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
#ifndef QCD_PSEUDOFERMION_TWO_FLAVOUR_EVEN_ODD_H
#define QCD_PSEUDOFERMION_TWO_FLAVOUR_EVEN_ODD_H

namespace Grid{
  namespace QCD{



    ////////////////////////////////////////////////////////////////////////
    // Two flavour pseudofermion action for any EO prec dop
    ////////////////////////////////////////////////////////////////////////
    template<class Impl>
    class TwoFlavourEvenOddPseudoFermionAction : public Action<typename Impl::GaugeField> {

    public:

      INHERIT_IMPL_TYPES(Impl);

    private:
      
      FermionOperator<Impl> & FermOp;// the basic operator

      OperatorFunction<FermionField> &DerivativeSolver;
      OperatorFunction<FermionField> &ActionSolver;

      FermionField PhiOdd;   // the pseudo fermion field for this trajectory
      FermionField PhiEven;  // the pseudo fermion field for this trajectory

    public:
      /////////////////////////////////////////////////
      // Pass in required objects.
      /////////////////////////////////////////////////
      TwoFlavourEvenOddPseudoFermionAction(FermionOperator<Impl>  &Op, 
					 OperatorFunction<FermionField> & DS,
					 OperatorFunction<FermionField> & AS
					   ) : 
        FermOp(Op), 
	DerivativeSolver(DS), 
	ActionSolver(AS), 
        PhiEven(Op.FermionRedBlackGrid()),
	PhiOdd(Op.FermionRedBlackGrid())
		  {};
      
      //////////////////////////////////////////////////////////////////////////////////////
      // Push the gauge field in to the dops. Assume any BC's and smearing already applied
      //////////////////////////////////////////////////////////////////////////////////////
      virtual void refresh(const GaugeField &U, GridParallelRNG& pRNG) {

	// P(phi) = e^{- phi^dag (MpcdagMpc)^-1 phi}
	// Phi = McpDag eta 
	// P(eta) = e^{- eta^dag eta}
	//
	// e^{x^2/2 sig^2} => sig^2 = 0.5.

	RealD scale = std::sqrt(0.5);

	FermionField eta    (FermOp.FermionGrid());
	FermionField etaOdd (FermOp.FermionRedBlackGrid());
	FermionField etaEven(FermOp.FermionRedBlackGrid());

	gaussian(pRNG,eta);
	pickCheckerboard(Even,etaEven,eta);
	pickCheckerboard(Odd,etaOdd,eta);

	FermOp.ImportGauge(U);
	SchurDifferentiableOperator<Impl> PCop(FermOp);
	

	PCop.MpcDag(etaOdd,PhiOdd);

	FermOp.MooeeDag(etaEven,PhiEven);

	PhiOdd =PhiOdd*scale;
	PhiEven=PhiEven*scale;
	
      };

      //////////////////////////////////////////////////////
      // S = phi^dag (Mdag M)^-1 phi  (odd)
      //   + phi^dag (Mdag M)^-1 phi  (even)
      //////////////////////////////////////////////////////
      virtual RealD S(const GaugeField &U) {

	FermOp.ImportGauge(U);

	FermionField X(FermOp.FermionRedBlackGrid());
	FermionField Y(FermOp.FermionRedBlackGrid());
	
	SchurDifferentiableOperator<Impl> PCop(FermOp);

	X=zero;
	ActionSolver(PCop,PhiOdd,X);
	PCop.Op(X,Y);
	RealD action = norm2(Y);

	// The EE factorised block; normally can replace with zero if det is constant (gauge field indept)
	// Only really clover term that creates this.
	FermOp.MooeeInvDag(PhiEven,Y);
	action = action + norm2(Y);

	std::cout << GridLogMessage << "Pseudofermion EO action "<<action<<std::endl;
	return action;
      };

      //////////////////////////////////////////////////////
      //
      // dS/du = - phi^dag  (Mdag M)^-1 [ Mdag dM + dMdag M ]  (Mdag M)^-1 phi
      //       = - phi^dag M^-1 dM (MdagM)^-1 phi -  phi^dag (MdagM)^-1 dMdag dM (Mdag)^-1 phi 
      //
      //       = - Ydag dM X  - Xdag dMdag Y
      //
      //////////////////////////////////////////////////////
      virtual void deriv(const GaugeField &U,GaugeField & dSdU) {

	FermOp.ImportGauge(U);

	FermionField X(FermOp.FermionRedBlackGrid());
	FermionField Y(FermOp.FermionRedBlackGrid());
	GaugeField tmp(FermOp.GaugeGrid());

	SchurDifferentiableOperator<Impl> Mpc(FermOp);

	// Our conventions really make this UdSdU; We do not differentiate wrt Udag here.
	// So must take dSdU - adj(dSdU) and left multiply by mom to get dS/dt.

	X=zero;
	DerivativeSolver(Mpc,PhiOdd,X);
	Mpc.Mpc(X,Y);
  	Mpc.MpcDeriv(tmp , Y, X );    dSdU=tmp;
	Mpc.MpcDagDeriv(tmp , X, Y);  dSdU=dSdU+tmp;

	// Treat the EE case. (MdagM)^-1 = Minv Minvdag
	// Deriv defaults to zero.
	//        FermOp.MooeeInvDag(PhiOdd,Y);
	//      FermOp.MooeeInv(Y,X);
	//	FermOp.MeeDeriv(tmp , Y, X,DaggerNo );    dSdU=tmp;
	//  FermOp.MeeDeriv(tmp , X, Y,DaggerYes);  dSdU=dSdU+tmp;

	assert(FermOp.ConstEE() == 1);

	/*
        FermOp.MooeeInvDag(PhiOdd,Y);
        FermOp.MooeeInv(Y,X);
  	FermOp.MeeDeriv(tmp , Y, X,DaggerNo );    dSdU=tmp;
	FermOp.MeeDeriv(tmp , X, Y,DaggerYes);  dSdU=dSdU+tmp;
	*/
	
	dSdU = Ta(dSdU);

      };

    };
    
  }
}

#endif
