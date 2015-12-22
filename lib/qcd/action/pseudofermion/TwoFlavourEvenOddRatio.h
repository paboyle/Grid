#ifndef QCD_PSEUDOFERMION_TWO_FLAVOUR_EVEN_ODD_RATIO_H
#define QCD_PSEUDOFERMION_TWO_FLAVOUR_EVEN_ODD_RATIO_H

namespace Grid{
  namespace QCD{

    ///////////////////////////////////////
    // Two flavour ratio
    ///////////////////////////////////////
    template<class Impl>
    class TwoFlavourEvenOddRatioPseudoFermionAction : public Action<typename Impl::GaugeField> {
    public:
      INHERIT_IMPL_TYPES(Impl);

    private:
      FermionOperator<Impl> & NumOp;// the basic operator
      FermionOperator<Impl> & DenOp;// the basic operator

      OperatorFunction<FermionField> &DerivativeSolver;
      OperatorFunction<FermionField> &ActionSolver;

      FermionField PhiOdd;   // the pseudo fermion field for this trajectory
      FermionField PhiEven;  // the pseudo fermion field for this trajectory

    public:
      TwoFlavourEvenOddRatioPseudoFermionAction(FermionOperator<Impl>  &_NumOp, 
						FermionOperator<Impl>  &_DenOp, 
						OperatorFunction<FermionField> & DS,
						OperatorFunction<FermionField> & AS) :
      NumOp(_NumOp), 
      DenOp(_DenOp), 
      DerivativeSolver(DS), 
      ActionSolver(AS),
      PhiEven(_NumOp.FermionRedBlackGrid()),
      PhiOdd(_NumOp.FermionRedBlackGrid()) 
	{
	  conformable(_NumOp.FermionGrid(), _DenOp.FermionGrid());
	  conformable(_NumOp.FermionRedBlackGrid(), _DenOp.FermionRedBlackGrid());
	  conformable(_NumOp.GaugeGrid(), _DenOp.GaugeGrid());
	  conformable(_NumOp.GaugeRedBlackGrid(), _DenOp.GaugeRedBlackGrid());
	};
      
      virtual void refresh(const GaugeField &U, GridParallelRNG& pRNG) {

	// P(phi) = e^{- phi^dag Vpc (MpcdagMpc)^-1 Vpcdag phi}
	//
	// NumOp == V
	// DenOp == M
	//
	// Take phi_o = Vpcdag^{-1} Mpcdag eta_o  ; eta_o = Mpcdag^{-1} Vpcdag Phi
	//
	// P(eta_o) = e^{- eta_o^dag eta_o}
	//
	// e^{x^2/2 sig^2} => sig^2 = 0.5.
	// 
	RealD scale = std::sqrt(0.5);

	FermionField eta    (NumOp.FermionGrid());
	FermionField etaOdd (NumOp.FermionRedBlackGrid());
	FermionField etaEven(NumOp.FermionRedBlackGrid());
	FermionField tmp    (NumOp.FermionRedBlackGrid());

	gaussian(pRNG,eta);

	pickCheckerboard(Even,etaEven,eta);
	pickCheckerboard(Odd,etaOdd,eta);

	NumOp.ImportGauge(U);
	DenOp.ImportGauge(U);

	SchurDifferentiableOperator<Impl> Mpc(DenOp);
	SchurDifferentiableOperator<Impl> Vpc(NumOp);

	// Odd det factors
	Mpc.MpcDag(etaOdd,PhiOdd);
	tmp=zero;
	ActionSolver(Vpc,PhiOdd,tmp);
	Vpc.Mpc(tmp,PhiOdd);            

	// Even det factors
	DenOp.MooeeDag(etaEven,tmp);
	NumOp.MooeeInvDag(tmp,PhiEven);

	PhiOdd =PhiOdd*scale;
	PhiEven=PhiEven*scale;
	
      };

      //////////////////////////////////////////////////////
      // S = phi^dag V (Mdag M)^-1 Vdag phi
      //////////////////////////////////////////////////////
      virtual RealD S(const GaugeField &U) {

	NumOp.ImportGauge(U);
	DenOp.ImportGauge(U);

	SchurDifferentiableOperator<Impl> Mpc(DenOp);
	SchurDifferentiableOperator<Impl> Vpc(NumOp);

	FermionField X(NumOp.FermionRedBlackGrid());
	FermionField Y(NumOp.FermionRedBlackGrid());

	Vpc.MpcDag(PhiOdd,Y);           // Y= Vdag phi
	X=zero;
	ActionSolver(Mpc,Y,X);          // X= (MdagM)^-1 Vdag phi
	Mpc.Mpc(X,Y);                   // Y=  Mdag^-1 Vdag phi

	RealD action = norm2(Y);

	// The EE factorised block; normally can replace with zero if det is constant (gauge field indept)
	// Only really clover term that creates this. Leave the EE portion as a future to do to make most
	// rapid progresss on DWF for now.
	//
	NumOp.MooeeDag(PhiEven,X);
	DenOp.MooeeInvDag(X,Y);
	action = action + norm2(Y);

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

	SchurDifferentiableOperator<Impl> Mpc(DenOp);
	SchurDifferentiableOperator<Impl> Vpc(NumOp);

	FermionField  X(NumOp.FermionRedBlackGrid());
	FermionField  Y(NumOp.FermionRedBlackGrid());

	GaugeField   force(NumOp.GaugeGrid());	

	//Y=Vdag phi
	//X = (Mdag M)^-1 V^dag phi
	//Y = (Mdag)^-1 V^dag  phi
	Vpc.MpcDag(PhiOdd,Y);          // Y= Vdag phi
	X=zero;
	DerivativeSolver(Mpc,Y,X);     // X= (MdagM)^-1 Vdag phi
	Mpc.Mpc(X,Y);                  // Y=  Mdag^-1 Vdag phi

	// phi^dag V (Mdag M)^-1 dV^dag  phi
	Vpc.MpcDagDeriv(force , X, PhiOdd );  dSdU=force;
  
	// phi^dag dV (Mdag M)^-1 V^dag  phi
	Vpc.MpcDeriv(force , PhiOdd, X );  dSdU=dSdU+force;

	//    -    phi^dag V (Mdag M)^-1 Mdag dM   (Mdag M)^-1 V^dag  phi
	//    -    phi^dag V (Mdag M)^-1 dMdag M   (Mdag M)^-1 V^dag  phi
	Mpc.MpcDeriv(force,Y,X);   dSdU=dSdU-force;
	Mpc.MpcDagDeriv(force,X,Y);  dSdU=dSdU-force;

	// FIXME No force contribution from EvenEven assumed here
	// Needs a fix for clover.
	assert(NumOp.ConstEE() == 1);
	assert(DenOp.ConstEE() == 1);

	dSdU = -Ta(dSdU);

      };
    };
  }
}
#endif
