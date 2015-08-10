#ifndef QCD_PSEUDOFERMION_TWO_FLAVOUR_H
#define QCD_PSEUDOFERMION_TWO_FLAVOUR_H

namespace Grid{
  namespace QCD{

    // Placeholder comments:

    ///////////////////////////////////////
    // Two flavour ratio
    ///////////////////////////////////////
    // S = phi^dag V (Mdag M)^-1 V^dag  phi
    // dS/du = phi^dag dV (Mdag M)^-1 V^dag  phi
    //       - phi^dag V (Mdag M)^-1 [ Mdag dM + dMdag M ]  (Mdag M)^-1 V^dag  phi
    //       + phi^dag V (Mdag M)^-1 dV^dag  phi

    ///////////////////////////////////////
    // One flavour rational
    ///////////////////////////////////////

    // S_f = chi^dag *  N(M^dag*M)/D(M^dag*M) * chi
    //
    // Here, M is some operator 
    // N and D makeup the rat. poly 
    //
    // Need
    // dS_f/dU = chi^dag   P/Q d[N/D]  P/Q  chi
    //
    // Here N/D \sim R_{-1/2} ~ (M^dagM)^{-1/2}
    //
    // N/D is expressed as partial fraction expansion:
    //
    //           a0 + \sum_k ak/(M^dagM + bk)
    //
    // d[N/D] is then
    //
    //          \sum_k -ak [M^dagM+bk]^{-1}  [ dM^dag M + M^dag dM ] [M^dag M + bk]^{-1}
    //
    // Need
    //
    //       Mf Phi_k = [MdagM+bk]^{-1} Phi
    //       Mf Phi   = \sum_k ak [MdagM+bk]^{-1} Phi
    //
    // With these building blocks
    //
    //       dS/dU =  \sum_k -ak Mf Phi_k^dag      [ dM^dag M + M^dag dM ] Mf Phi_k
    //        S    = innerprodReal(Phi,Mf Phi);
    
    ///////////////////////////////////////
    // One flavour rational ratio
    ///////////////////////////////////////

    // S_f = chi^dag* P(V^dag*V)/Q(V^dag*V)* N(M^dag*M)/D(M^dag*M)* P(V^dag*V)/Q(V^dag*V)* chi
    //
    // Here, M is some 5D operator and V is the Pauli-Villars field
    // N and D makeup the rat. poly of the M term and P and & makeup the rat.poly of the denom term
    //
    // Need
    // dS_f/dU =  chi^dag d[P/Q]  N/D   P/Q  chi
    //         +  chi^dag   P/Q d[N/D]  P/Q  chi
    //         +  chi^dag   P/Q   N/D d[P/Q] chi
    //
    // Here P/Q \sim R_{1/4}  ~ (V^dagV)^{1/4}
    // Here N/D \sim R_{-1/2} ~ (M^dagM)^{-1/2}
    //
    // P/Q is expressed as partial fraction expansion:
    //
    //           a0 + \sum_k ak/(V^dagV + bk)
    //
    // d[P/Q] is then
    //
    //          \sum_k -ak [V^dagV+bk]^{-1}  [ dV^dag V + V^dag dV ] [V^dag V + bk]^{-1}
    //
    // and similar for N/D.
    // 
    // Need
    //       MpvPhi_k   = [Vdag V + bk]^{-1} chi
    //
    //       MpvPhi     = {a0 +  \sum_k ak [Vdag V + bk]^{-1} }chi
    //
    //       MfMpvPhi_k = [MdagM+bk]^{-1} MpvPhi
    //      
    //       MfMpvPhi   = {a0 +  \sum_k ak [Mdag M + bk]^{-1} } MpvPhi
    //
    //       MpvMfMpvPhi_k = [Vdag V + bk]^{-1} MfMpvchi
    //
    // With these building blocks
    //
    //       dS/dU =  
    //                 \sum_k -ak MpvPhi_k^dag        [ dV^dag V + V^dag dV ] MpvMfMpvPhi_k           <- deriv on P left
    //             +   \sum_k -ak MpvMfMpvPhi_k^\dag  [ dV^dag V + V^dag dV ] MpvPhi_k
    //             +   \sum_k -ak MfMpvPhi_k^dag      [ dM^dag M + M^dag dM ] MfMpvPhi_k

    
    ////////////////////////////////////////////////////////////////////////
    // Two flavour pseudofermion action for any dop
    ////////////////////////////////////////////////////////////////////////
    template<class Impl>
    class TwoFlavourPseudoFermionAction : public Action<typename Impl::GaugeField> {

    private:
#include <qcd/action/fermion/FermionImplTypedefs.h>
      
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
      virtual void init(const GaugeField &U, GridParallelRNG& pRNG) {

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
