#ifndef QCD_PSEUDOFERMION_TWO_FLAVOUR_EVEN_ODD_H
#define QCD_PSEUDOFERMION_TWO_FLAVOUR_EVEN_ODD_H

namespace Grid{
  namespace QCD{

    template<class Matrix,class FermionField>
      class SchurDifferentiableOperator :  public SchurDiagMooeeOperator<Matrix,FermionField> 
      {
    public:
	SchurDifferentiableOperator (Matrix &Mat) : SchurDiagMooeeOperator<Matrix,FermionField>(Mat) {};

	void MpcDeriv(LatticeGaugeField &Force,const FermionField &U,const FermionField &V) {
	
	  GridBase *fgrid   = this->_Mat.FermionGrid();
	  GridBase *fcbgrid = this->_Mat.FermionRedBlackGrid();
	  GridBase *ugrid = this->_Mat.GaugeGrid();
	  GridBase *ucbgrid = this->_Mat.GaugeRedBlackGrid();

	  Real coeff = 1.0;
	  FermionField tmp1(fcbgrid);
	  FermionField tmp2(fcbgrid);

	  conformable(fcbgrid,U._grid);
	  conformable(fcbgrid,V._grid);

	  // Assert the checkerboard?? or code for either
	  assert(U.checkerboard==Odd);
	  assert(V.checkerboard==V.checkerboard);

	  LatticeGaugeField ForceO(ucbgrid);
	  LatticeGaugeField ForceE(ucbgrid);

	  //  X^dag Der_oe MeeInv Meo Y
	  // Use Mooee as nontrivial but gauge field indept
	  this->_Mat.Meooe   (V,tmp1);      // odd->even -- implicit -0.5 factor to be applied
	  this->_Mat.MooeeInv(tmp1,tmp2);   // even->even 
	  this->_Mat.MoeDeriv(ForceO,U,tmp2,DaggerNo);
	  
	  //  Accumulate X^dag M_oe MeeInv Der_eo Y
	  this->_Mat.MeooeDag   (U,tmp1);    // even->odd -- implicit -0.5 factor to be applied
	  this->_Mat.MooeeInvDag(tmp1,tmp2); // even->even 
	  this->_Mat.MeoDeriv(ForceE,tmp2,V,DaggerNo);
	  
	  setCheckerboard(Force,ForceE); 
	  setCheckerboard(Force,ForceO);
	  Force=-Force;
	}


	void MpcDagDeriv(LatticeGaugeField &Force,const FermionField &U,const FermionField &V) {
	
	  GridBase *fgrid   = this->_Mat.FermionGrid();
	  GridBase *fcbgrid = this->_Mat.FermionRedBlackGrid();
	  GridBase *ugrid   = this->_Mat.GaugeGrid();
	  GridBase *ucbgrid = this->_Mat.GaugeRedBlackGrid();

	  Real coeff = 1.0;
	  FermionField tmp1(fcbgrid);
	  FermionField tmp2(fcbgrid);

	  conformable(fcbgrid,U._grid);
	  conformable(fcbgrid,V._grid);

	  // Assert the checkerboard?? or code for either
	  assert(V.checkerboard==Odd);
	  assert(V.checkerboard==V.checkerboard);

	  LatticeGaugeField ForceO(ucbgrid);
	  LatticeGaugeField ForceE(ucbgrid);

	  //  X^dag Der_oe MeeInv Meo Y
	  // Use Mooee as nontrivial but gauge field indept
	  this->_Mat.MeooeDag   (V,tmp1);      // odd->even -- implicit -0.5 factor to be applied
	  this->_Mat.MooeeInvDag(tmp1,tmp2);   // even->even 
	  this->_Mat.MoeDeriv(ForceO,U,tmp2,DaggerYes);
	  
	  //  Accumulate X^dag M_oe MeeInv Der_eo Y
	  this->_Mat.Meooe   (U,tmp1);    // even->odd -- implicit -0.5 factor to be applied
	  this->_Mat.MooeeInv(tmp1,tmp2); // even->even 
	  this->_Mat.MeoDeriv(ForceE,tmp2,V,DaggerYes);

	  setCheckerboard(Force,ForceE); 
	  setCheckerboard(Force,ForceO);
	  Force=-Force;
	}

    };


    ////////////////////////////////////////////////////////////////////////
    // Two flavour pseudofermion action for any EO prec dop
    ////////////////////////////////////////////////////////////////////////
    template<class GaugeField,class MatrixField,class FermionField>
      class TwoFlavourEvenOddPseudoFermionAction : public Action<GaugeField> {

    private:
      
      FermionOperator<FermionField,GaugeField> & FermOp;// the basic operator

      OperatorFunction<FermionField> &DerivativeSolver;

      OperatorFunction<FermionField> &ActionSolver;

      FermionField PhiOdd;   // the pseudo fermion field for this trajectory
      FermionField PhiEven;  // the pseudo fermion field for this trajectory

    public:
      /////////////////////////////////////////////////
      // Pass in required objects.
      /////////////////////////////////////////////////
      TwoFlavourEvenOddPseudoFermionAction(FermionOperator<FermionField,GaugeField>  &Op, 
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
      virtual void init(const GaugeField &U, GridParallelRNG& pRNG) {

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

	SchurDifferentiableOperator<FermionOperator<FermionField,GaugeField>,FermionField> PCop(FermOp);
	
	FermOp.ImportGauge(U);

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
	
	SchurDifferentiableOperator<FermionOperator<FermionField,GaugeField>,FermionField> PCop(FermOp);

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

	SchurDifferentiableOperator<FermionOperator<FermionField,GaugeField>,FermionField> PCop(FermOp);

	X=zero;
	DerivativeSolver(PCop,PhiOdd,X);
	PCop.Op(X,Y);

	// Our conventions really make this UdSdU; We do not differentiate wrt Udag here.
	// So must take dSdU - adj(dSdU) and left multiply by mom to get dS/dt.

  	PCop.MpcDeriv(tmp , Y, X );    dSdU=tmp;
	PCop.MpcDagDeriv(tmp , X, Y);  dSdU=dSdU+tmp;

	// Treat the EE case. (MdagM)^-1 = Minv Minvdag
	// Deriv defaults to zero.
        FermOp.MooeeInvDag(PhiOdd,Y);
        FermOp.MooeeInv(Y,X);
  	FermOp.MeeDeriv(tmp , Y, X,DaggerNo );    dSdU=tmp;
	FermOp.MeeDeriv(tmp , X, Y,DaggerYes);  dSdU=dSdU+tmp;
	
	dSdU = Ta(dSdU);

      };

    };
    
  }
}

#endif
