    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/action/pseudofermion/GeneralEvenOddRationalRatio.h

    Copyright (C) 2015

    Author: Christopher Kelly <ckelly@bnl.gov>
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
#ifndef QCD_PSEUDOFERMION_GENERAL_EVEN_ODD_RATIONAL_RATIO_H
#define QCD_PSEUDOFERMION_GENERAL_EVEN_ODD_RATIONAL_RATIO_H

NAMESPACE_BEGIN(Grid);

    /////////////////////////////////////////////////////////
    // Generic rational approximation for ratios of operators
    /////////////////////////////////////////////////////////

    /* S_f = -log( det(  [M^dag M]/[V^dag V] )^{1/inv_pow}  )
           = chi^dag ( [M^dag M]/[V^dag V] )^{-1/inv_pow} chi\
	   = chi^dag ( [V^dag V]^{-1/2} [M^dag M] [V^dag V]^{-1/2} )^{-1/inv_pow} chi\
	   = chi^dag [V^dag V]^{1/(2*inv_pow)} [M^dag M]^{-1/inv_pow} [V^dag V]^{-1/(2*inv_pow)} chi\

	   S_f = chi^dag* P(V^dag*V)/Q(V^dag*V)* N(M^dag*M)/D(M^dag*M)* P(V^dag*V)/Q(V^dag*V)* chi       
    
       BIG WARNING:	   
       Here V^dag V is referred to in this code as the "numerator" operator and M^dag M is the *denominator* operator.
       this refers to their position in the pseudofermion action, which is the *inverse* of what appears in the determinant
       Thus for DWF the numerator operator is the Pauli-Villars operator

       Here P/Q \sim R_{1/(2*inv_pow)}  ~ (V^dagV)^{1/(2*inv_pow)}  
       Here N/D \sim R_{-1/inv_pow} ~ (M^dagM)^{-1/inv_pow}  
    */
      
    template<class Impl>
    class GeneralEvenOddRatioRationalPseudoFermionAction : public Action<typename Impl::GaugeField> {
    public:

      INHERIT_IMPL_TYPES(Impl);

      typedef RationalActionParams Params;
      Params param;

      MultiShiftFunction ApproxPower   ;  //rational approx for X^{1/inv_pow}
      MultiShiftFunction ApproxNegPower;  //rational approx for X^{-1/inv_pow}
      MultiShiftFunction ApproxHalfPower;   //rational approx for X^{1/(2*inv_pow)}
      MultiShiftFunction ApproxNegHalfPower; //rational approx for X^{-1/(2*inv_pow)}

    private:
     
      FermionOperator<Impl> & NumOp;// the basic operator
      FermionOperator<Impl> & DenOp;// the basic operator
      FermionField PhiEven; // the pseudo fermion field for this trajectory
      FermionField PhiOdd; // the pseudo fermion field for this trajectory

    public:

      GeneralEvenOddRatioRationalPseudoFermionAction(FermionOperator<Impl>  &_NumOp, 
						     FermionOperator<Impl>  &_DenOp, 
						     Params & p
						     ) : 
	NumOp(_NumOp), 
	DenOp(_DenOp), 
	PhiOdd (_NumOp.FermionRedBlackGrid()),
	PhiEven(_NumOp.FermionRedBlackGrid()),
	param(p) 
      {
	AlgRemez remez(param.lo,param.hi,param.precision);

	int inv_pow = param.inv_pow;
	int _2_inv_pow = 2*inv_pow;

	// MdagM^(+- 1/inv_pow)
	std::cout<<GridLogMessage << "Generating degree "<<param.degree<<" for x^(1/" << inv_pow << ")"<<std::endl;
	remez.generateApprox(param.degree,1,inv_pow);
	ApproxPower.Init(remez,param.tolerance,false);
	ApproxNegPower.Init(remez,param.tolerance,true);

	// VdagV^(+- 1/(2*inv_pow))
	std::cout<<GridLogMessage << "Generating degree "<<param.degree<<" for x^(1/" << _2_inv_pow << ")"<<std::endl;
	remez.generateApprox(param.degree,1,_2_inv_pow);
   	ApproxHalfPower.Init(remez,param.tolerance,false);
	ApproxNegHalfPower.Init(remez,param.tolerance,true);
      };

      virtual std::string action_name(){return "GeneralEvenOddRatioRationalPseudoFermionAction";}

      virtual std::string LogParameters(){
	std::stringstream sstream;
	sstream << GridLogMessage << "["<<action_name()<<"] Power          : 1/" << param.inv_pow <<  std::endl;
	sstream << GridLogMessage << "["<<action_name()<<"] Low            :" << param.lo <<  std::endl;
	sstream << GridLogMessage << "["<<action_name()<<"] High           :" << param.hi <<  std::endl;
	sstream << GridLogMessage << "["<<action_name()<<"] Max iterations :" << param.MaxIter <<  std::endl;
	sstream << GridLogMessage << "["<<action_name()<<"] Tolerance      :" << param.tolerance <<  std::endl;
	sstream << GridLogMessage << "["<<action_name()<<"] Degree         :" << param.degree <<  std::endl;
	sstream << GridLogMessage << "["<<action_name()<<"] Precision      :" << param.precision <<  std::endl;
	return sstream.str();
      }
      
      
      virtual void refresh(const GaugeField &U, GridParallelRNG& pRNG) {

	// S_f = chi^dag* P(V^dag*V)/Q(V^dag*V)* N(M^dag*M)/D(M^dag*M)* P(V^dag*V)/Q(V^dag*V)* chi       
	//
	// P(phi) = e^{- phi^dag (VdagV)^1/(2*inv_pow) (MdagM)^-1/inv_pow (VdagV)^1/(2*inv_pow) phi}
	//        = e^{- phi^dag  (VdagV)^1/(2*inv_pow) (MdagM)^-1/(2*inv_pow) (MdagM)^-1/(2*inv_pow)  (VdagV)^1/(2*inv_pow) phi}
	//
	// Phi =  (VdagV)^-1/(2*inv_pow) Mdag^{1/(2*inv_pow)} eta 
	//
	// P(eta) = e^{- eta^dag eta}
	//	
	// General gaussian random draws from   e^{x^2/(2 sig^2)} => sig^2 = 0.5.
	// 
	// So eta should be of width sig = 1/sqrt(2).

	RealD scale = std::sqrt(0.5);

	FermionField eta(NumOp.FermionGrid());
	FermionField etaOdd (NumOp.FermionRedBlackGrid());
	FermionField etaEven(NumOp.FermionRedBlackGrid());
	FermionField     tmp(NumOp.FermionRedBlackGrid());

	gaussian(pRNG,eta);	eta=eta*scale;

	pickCheckerboard(Even,etaEven,eta);
	pickCheckerboard(Odd,etaOdd,eta);

	NumOp.ImportGauge(U);
	DenOp.ImportGauge(U);


	// MdagM^1/(2*inv_pow) eta
	SchurDifferentiableOperator<Impl> MdagM(DenOp);
	ConjugateGradientMultiShift<FermionField> msCG_M(param.MaxIter,ApproxHalfPower);
	msCG_M(MdagM,etaOdd,tmp);

	// VdagV^-1/(2*inv_pow) MdagM^1/(2*inv_pow) eta
	SchurDifferentiableOperator<Impl> VdagV(NumOp);
	ConjugateGradientMultiShift<FermionField> msCG_V(param.MaxIter,ApproxNegHalfPower);
	msCG_V(VdagV,tmp,PhiOdd);

	assert(NumOp.ConstEE() == 1);
	assert(DenOp.ConstEE() == 1);
	PhiEven = Zero();
	
      };

      //////////////////////////////////////////////////////
      // S_f = chi^dag* P(V^dag*V)/Q(V^dag*V)* N(M^dag*M)/D(M^dag*M)* P(V^dag*V)/Q(V^dag*V)* chi       
      //////////////////////////////////////////////////////
      virtual RealD S(const GaugeField &U) {

	NumOp.ImportGauge(U);
	DenOp.ImportGauge(U);

	FermionField X(NumOp.FermionRedBlackGrid());
	FermionField Y(NumOp.FermionRedBlackGrid());

	// VdagV^1/(2*inv_pow) Phi
	SchurDifferentiableOperator<Impl> VdagV(NumOp);
	ConjugateGradientMultiShift<FermionField> msCG_V(param.MaxIter,ApproxHalfPower);
	msCG_V(VdagV,PhiOdd,X);

	// MdagM^-1/(2*inv_pow) VdagV^1/(2*inv_pow) Phi
	SchurDifferentiableOperator<Impl> MdagM(DenOp);
	ConjugateGradientMultiShift<FermionField> msCG_M(param.MaxIter,ApproxNegHalfPower);
	msCG_M(MdagM,X,Y);

	// Randomly apply rational bounds checks.
	if ( (rand()%param.BoundsCheckFreq)==0 ) { 
	  FermionField gauss(NumOp.FermionRedBlackGrid());
	  gauss = PhiOdd;
	  HighBoundCheck(MdagM,gauss,param.hi);

	  InversePowerBoundsCheck(param.inv_pow,param.MaxIter,param.tolerance*100,MdagM,gauss,ApproxNegPower);
	}

	//  Phidag VdagV^1/(2*inv_pow) MdagM^-1/(2*inv_pow)  MdagM^-1/(2*inv_pow) VdagV^1/(2*inv_pow) Phi
	RealD action = norm2(Y);

	return action;
      };

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
      //       MpvPhi     = {a0 +  \sum_k ak [Vdag V + bk]^{-1} }chi   
      //   
      //       MfMpvPhi_k = [MdagM+bk]^{-1} MpvPhi  
      //       MfMpvPhi   = {a0 +  \sum_k ak [Mdag M + bk]^{-1} } MpvPhi
      // 
      //       MpvMfMpvPhi_k = [Vdag V + bk]^{-1} MfMpvchi   
      //  

      virtual void deriv(const GaugeField &U,GaugeField & dSdU) {

	const int n_f  = ApproxNegPower.poles.size();
	const int n_pv = ApproxHalfPower.poles.size();

	std::vector<FermionField> MpvPhi_k     (n_pv,NumOp.FermionRedBlackGrid());
	std::vector<FermionField> MpvMfMpvPhi_k(n_pv,NumOp.FermionRedBlackGrid());
	std::vector<FermionField> MfMpvPhi_k   (n_f ,NumOp.FermionRedBlackGrid());

	FermionField      MpvPhi(NumOp.FermionRedBlackGrid());
	FermionField    MfMpvPhi(NumOp.FermionRedBlackGrid());
	FermionField MpvMfMpvPhi(NumOp.FermionRedBlackGrid());
	FermionField           Y(NumOp.FermionRedBlackGrid());

	GaugeField   tmp(NumOp.GaugeGrid());

	NumOp.ImportGauge(U);
	DenOp.ImportGauge(U);

	SchurDifferentiableOperator<Impl> VdagV(NumOp);
	SchurDifferentiableOperator<Impl> MdagM(DenOp);

	ConjugateGradientMultiShift<FermionField> msCG_V(param.MaxIter,ApproxHalfPower);
	ConjugateGradientMultiShift<FermionField> msCG_M(param.MaxIter,ApproxNegPower);

	msCG_V(VdagV,PhiOdd,MpvPhi_k,MpvPhi);
	msCG_M(MdagM,MpvPhi,MfMpvPhi_k,MfMpvPhi);
	msCG_V(VdagV,MfMpvPhi,MpvMfMpvPhi_k,MpvMfMpvPhi);

	RealD ak;

	dSdU = Zero();

	// With these building blocks  
	//  
	//       dS/dU = 
	//                 \sum_k -ak MfMpvPhi_k^dag      [ dM^dag M + M^dag dM ] MfMpvPhi_k         (1)
	//             +   \sum_k -ak MpvMfMpvPhi_k^\dag  [ dV^dag V + V^dag dV ] MpvPhi_k           (2)
	//                        -ak MpvPhi_k^dag        [ dV^dag V + V^dag dV ] MpvMfMpvPhi_k      (3)

	//(1)
	for(int k=0;k<n_f;k++){
	  ak = ApproxNegPower.residues[k];
	  MdagM.Mpc(MfMpvPhi_k[k],Y);
	  MdagM.MpcDagDeriv(tmp , MfMpvPhi_k[k], Y );  dSdU=dSdU+ak*tmp;
	  MdagM.MpcDeriv(tmp , Y, MfMpvPhi_k[k] );  dSdU=dSdU+ak*tmp;
	}
	
	//(2)
	//(3)
	for(int k=0;k<n_pv;k++){

          ak = ApproxHalfPower.residues[k];
	  
	  VdagV.Mpc(MpvPhi_k[k],Y);
	  VdagV.MpcDagDeriv(tmp,MpvMfMpvPhi_k[k],Y); dSdU=dSdU+ak*tmp;
	  VdagV.MpcDeriv   (tmp,Y,MpvMfMpvPhi_k[k]);  dSdU=dSdU+ak*tmp;     
	  
	  VdagV.Mpc(MpvMfMpvPhi_k[k],Y);                // V as we take Ydag 
	  VdagV.MpcDeriv   (tmp,Y, MpvPhi_k[k]); dSdU=dSdU+ak*tmp;
	  VdagV.MpcDagDeriv(tmp,MpvPhi_k[k], Y); dSdU=dSdU+ak*tmp;

	}

	//dSdU = Ta(dSdU);

      };
    };

NAMESPACE_END(Grid);

#endif
