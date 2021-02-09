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

      //For action evaluation
      MultiShiftFunction ApproxPowerAction   ;  //rational approx for X^{1/inv_pow}
      MultiShiftFunction ApproxNegPowerAction;  //rational approx for X^{-1/inv_pow}
      MultiShiftFunction ApproxHalfPowerAction;   //rational approx for X^{1/(2*inv_pow)}
      MultiShiftFunction ApproxNegHalfPowerAction; //rational approx for X^{-1/(2*inv_pow)}

      //For the MD integration
      MultiShiftFunction ApproxPowerMD   ;  //rational approx for X^{1/inv_pow}
      MultiShiftFunction ApproxNegPowerMD;  //rational approx for X^{-1/inv_pow}
      MultiShiftFunction ApproxHalfPowerMD;   //rational approx for X^{1/(2*inv_pow)}
      MultiShiftFunction ApproxNegHalfPowerMD; //rational approx for X^{-1/(2*inv_pow)}

    private:
     
      FermionOperator<Impl> & NumOp;// the basic operator
      FermionOperator<Impl> & DenOp;// the basic operator
      FermionField PhiEven; // the pseudo fermion field for this trajectory
      FermionField PhiOdd; // the pseudo fermion field for this trajectory

      //Generate the approximation to x^{1/inv_pow} (->approx)   and x^{-1/inv_pow} (-> approx_inv)  by an approx_degree degree rational approximation
      //CG_tolerance is used to issue a warning if the approximation error is larger than the tolerance of the CG and is otherwise just stored in the MultiShiftFunction for use by the multi-shift
      static void generateApprox(MultiShiftFunction &approx, MultiShiftFunction &approx_inv, int inv_pow, int approx_degree, double CG_tolerance, AlgRemez &remez){
	std::cout<<GridLogMessage << "Generating degree "<< approx_degree<<" approximation for x^(1/" << inv_pow << ")"<<std::endl;
	double error = remez.generateApprox(approx_degree,1,inv_pow);	
	if(error > CG_tolerance)
	  std::cout<<GridLogMessage << "WARNING: Remez approximation has a larger error " << error << " than the CG tolerance " << CG_tolerance << "! Try increasing the number of poles" << std::endl;
	
	approx.Init(remez, CG_tolerance,false);
	approx_inv.Init(remez, CG_tolerance,true);
      }


    protected:
      static constexpr bool Numerator = true;
      static constexpr bool Denominator = false;

      //Allow derived classes to override the multishift CG
      virtual void multiShiftInverse(bool numerator, const MultiShiftFunction &approx, const Integer MaxIter, const FermionField &in, FermionField &out){
	SchurDifferentiableOperator<Impl> schurOp(numerator ? NumOp : DenOp);
	ConjugateGradientMultiShift<FermionField> msCG(MaxIter, approx);
	msCG(schurOp,in, out);
      }
      virtual void multiShiftInverse(bool numerator, const MultiShiftFunction &approx, const Integer MaxIter, const FermionField &in, std::vector<FermionField> &out_elems, FermionField &out){
	SchurDifferentiableOperator<Impl> schurOp(numerator ? NumOp : DenOp);
	ConjugateGradientMultiShift<FermionField> msCG(MaxIter, approx);
	msCG(schurOp,in, out_elems, out);
      }
      //Allow derived classes to override the gauge import
      virtual void ImportGauge(const GaugeField &U){
	NumOp.ImportGauge(U);
	DenOp.ImportGauge(U);
      }
      
    public:

      GeneralEvenOddRatioRationalPseudoFermionAction(FermionOperator<Impl>  &_NumOp, 
						     FermionOperator<Impl>  &_DenOp, 
						     const Params & p
						     ) : 
	NumOp(_NumOp), 
	DenOp(_DenOp), 
	PhiOdd (_NumOp.FermionRedBlackGrid()),
	PhiEven(_NumOp.FermionRedBlackGrid()),
	param(p) 
      {
	std::cout<<GridLogMessage << action_name() << " initialize: starting" << std::endl;
	AlgRemez remez(param.lo,param.hi,param.precision);

	//Generate approximations for action eval
	generateApprox(ApproxPowerAction, ApproxNegPowerAction, param.inv_pow, param.action_degree, param.action_tolerance, remez);
	generateApprox(ApproxHalfPowerAction, ApproxNegHalfPowerAction, 2*param.inv_pow, param.action_degree, param.action_tolerance, remez);

	//Generate approximations for MD
	if(param.md_degree != param.action_degree){ //note the CG tolerance is unrelated to the stopping condition of the Remez algorithm
	  generateApprox(ApproxPowerMD, ApproxNegPowerMD, param.inv_pow, param.md_degree, param.md_tolerance, remez);
	  generateApprox(ApproxHalfPowerMD, ApproxNegHalfPowerMD, 2*param.inv_pow, param.md_degree, param.md_tolerance, remez);
	}else{
	  std::cout<<GridLogMessage << "Using same rational approximations for MD as for action evaluation" << std::endl;
	  ApproxPowerMD = ApproxPowerAction; 
	  ApproxNegPowerMD = ApproxNegPowerAction;
	  for(int i=0;i<ApproxPowerMD.tolerances.size();i++)
	    ApproxNegPowerMD.tolerances[i] = ApproxPowerMD.tolerances[i] = param.md_tolerance; //used for multishift

	  ApproxHalfPowerMD = ApproxHalfPowerAction;
	  ApproxNegHalfPowerMD = ApproxNegHalfPowerAction;
	  for(int i=0;i<ApproxPowerMD.tolerances.size();i++)
	    ApproxNegHalfPowerMD.tolerances[i] = ApproxHalfPowerMD.tolerances[i] = param.md_tolerance;
	}

	std::cout<<GridLogMessage << action_name() << " initialize: complete" << std::endl;
      };

      virtual std::string action_name(){return "GeneralEvenOddRatioRationalPseudoFermionAction";}

      virtual std::string LogParameters(){
	std::stringstream sstream;
	sstream << GridLogMessage << "["<<action_name()<<"] Power              : 1/" << param.inv_pow <<  std::endl;
	sstream << GridLogMessage << "["<<action_name()<<"] Low                :" << param.lo <<  std::endl;
	sstream << GridLogMessage << "["<<action_name()<<"] High               :" << param.hi <<  std::endl;
	sstream << GridLogMessage << "["<<action_name()<<"] Max iterations     :" << param.MaxIter <<  std::endl;
	sstream << GridLogMessage << "["<<action_name()<<"] Tolerance (Action) :" << param.action_tolerance <<  std::endl;
	sstream << GridLogMessage << "["<<action_name()<<"] Degree (Action)    :" << param.action_degree <<  std::endl;
	sstream << GridLogMessage << "["<<action_name()<<"] Tolerance (MD)     :" << param.md_tolerance <<  std::endl;
	sstream << GridLogMessage << "["<<action_name()<<"] Degree (MD)        :" << param.md_degree <<  std::endl;
	sstream << GridLogMessage << "["<<action_name()<<"] Precision          :" << param.precision <<  std::endl;
	return sstream.str();
      }

      //Access the fermion field
      const FermionField &getPhiOdd() const{ return PhiOdd; }
      
      virtual void refresh(const GaugeField &U, GridParallelRNG& pRNG) {

	// S_f = chi^dag* P(V^dag*V)/Q(V^dag*V)* N(M^dag*M)/D(M^dag*M)* P(V^dag*V)/Q(V^dag*V)* chi       
	//
	// P(phi) = e^{- phi^dag (VdagV)^1/(2*inv_pow) (MdagM)^-1/inv_pow (VdagV)^1/(2*inv_pow) phi}
	//        = e^{- phi^dag  (VdagV)^1/(2*inv_pow) (MdagM)^-1/(2*inv_pow) (MdagM)^-1/(2*inv_pow)  (VdagV)^1/(2*inv_pow) phi}
	//
	// Phi =  (VdagV)^-1/(2*inv_pow) Mdag^{1/(2*inv_pow)} eta 
	
	std::cout<<GridLogMessage << action_name() << " refresh: starting" << std::endl;

	FermionField eta(NumOp.FermionGrid());
	FermionField etaOdd (NumOp.FermionRedBlackGrid());
	FermionField etaEven(NumOp.FermionRedBlackGrid());
	FermionField     tmp(NumOp.FermionRedBlackGrid());

	// P(eta) \propto e^{- eta^dag eta}
	//	
	// The gaussian function draws from  P(x) \propto e^{- x^2 / 2 }    [i.e. sigma=1]
	// Thus eta = x/sqrt{2} = x * sqrt(1/2)
	RealD scale = std::sqrt(0.5);
	gaussian(pRNG,eta);	eta=eta*scale;

	pickCheckerboard(Even,etaEven,eta);
	pickCheckerboard(Odd,etaOdd,eta);

	ImportGauge(U);

	// MdagM^1/(2*inv_pow) eta
	std::cout<<GridLogMessage << action_name() << " refresh: doing (M^dag M)^{1/" << 2*param.inv_pow << "} eta" << std::endl;
	multiShiftInverse(Denominator, ApproxHalfPowerAction, param.MaxIter, etaOdd, tmp);

	// VdagV^-1/(2*inv_pow) MdagM^1/(2*inv_pow) eta
	std::cout<<GridLogMessage << action_name() << " refresh: doing (V^dag V)^{-1/" << 2*param.inv_pow << "} ( (M^dag M)^{1/" << 2*param.inv_pow << "} eta)" << std::endl;
	multiShiftInverse(Numerator, ApproxNegHalfPowerAction, param.MaxIter, tmp, PhiOdd);
		
	assert(NumOp.ConstEE() == 1);
	assert(DenOp.ConstEE() == 1);
	PhiEven = Zero();
	std::cout<<GridLogMessage << action_name() << " refresh: starting" << std::endl;
      };

      //////////////////////////////////////////////////////
      // S_f = chi^dag* P(V^dag*V)/Q(V^dag*V)* N(M^dag*M)/D(M^dag*M)* P(V^dag*V)/Q(V^dag*V)* chi       
      //////////////////////////////////////////////////////
      virtual RealD S(const GaugeField &U) {
	std::cout<<GridLogMessage << action_name() << " compute action: starting" << std::endl;
	ImportGauge(U);

	FermionField X(NumOp.FermionRedBlackGrid());
	FermionField Y(NumOp.FermionRedBlackGrid());

	// VdagV^1/(2*inv_pow) Phi
	std::cout<<GridLogMessage << action_name() << " compute action: doing (V^dag V)^{1/" << 2*param.inv_pow << "} Phi" << std::endl;
	multiShiftInverse(Numerator, ApproxHalfPowerAction, param.MaxIter, PhiOdd,X);

	// MdagM^-1/(2*inv_pow) VdagV^1/(2*inv_pow) Phi
	std::cout<<GridLogMessage << action_name() << " compute action: doing (M^dag M)^{-1/" << 2*param.inv_pow << "} ( (V^dag V)^{1/" << 2*param.inv_pow << "} Phi)" << std::endl;
	multiShiftInverse(Denominator, ApproxNegHalfPowerAction, param.MaxIter, X,Y);

	// Randomly apply rational bounds checks.
	int rcheck = rand();
	CartesianCommunicator::BroadcastWorld(0,(void *)&rcheck,sizeof(int)); //make sure all nodes have the same number or you will sporadically hang and spend days trying to find out why (trust me - CK)

	if ( param.BoundsCheckFreq != 0 && (rcheck % param.BoundsCheckFreq)==0 ) { 
	  std::cout<<GridLogMessage << action_name() << " compute action: doing bounds check" << std::endl;
	  FermionField gauss(NumOp.FermionRedBlackGrid());
	  gauss = PhiOdd;
	  SchurDifferentiableOperator<Impl> MdagM(DenOp);
	  std::cout<<GridLogMessage << action_name() << " compute action: checking high bounds" << std::endl;
	  HighBoundCheck(MdagM,gauss,param.hi);
	  std::cout<<GridLogMessage << action_name() << " compute action: full approximation" << std::endl;
	  InversePowerBoundsCheck(param.inv_pow,param.MaxIter,param.action_tolerance*100,MdagM,gauss,ApproxNegPowerAction);
	  std::cout<<GridLogMessage << action_name() << " compute action: bounds check complete" << std::endl;
	}

	//  Phidag VdagV^1/(2*inv_pow) MdagM^-1/(2*inv_pow)  MdagM^-1/(2*inv_pow) VdagV^1/(2*inv_pow) Phi
	RealD action = norm2(Y);
	std::cout<<GridLogMessage << action_name() << " compute action: complete" << std::endl;

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
	std::cout<<GridLogMessage << action_name() << " deriv: starting" << std::endl;
	const int n_f  = ApproxNegPowerMD.poles.size();
	const int n_pv = ApproxHalfPowerMD.poles.size();

	std::vector<FermionField> MpvPhi_k     (n_pv,NumOp.FermionRedBlackGrid());
	std::vector<FermionField> MpvMfMpvPhi_k(n_pv,NumOp.FermionRedBlackGrid());
	std::vector<FermionField> MfMpvPhi_k   (n_f ,NumOp.FermionRedBlackGrid());

	FermionField      MpvPhi(NumOp.FermionRedBlackGrid());
	FermionField    MfMpvPhi(NumOp.FermionRedBlackGrid());
	FermionField MpvMfMpvPhi(NumOp.FermionRedBlackGrid());
	FermionField           Y(NumOp.FermionRedBlackGrid());

	GaugeField   tmp(NumOp.GaugeGrid());

	ImportGauge(U);

	std::cout<<GridLogMessage << action_name() << " deriv: doing (V^dag V)^{1/" << 2*param.inv_pow << "} Phi" << std::endl;
	multiShiftInverse(Numerator, ApproxHalfPowerMD, param.MaxIter, PhiOdd,MpvPhi_k,MpvPhi);

	std::cout<<GridLogMessage << action_name() << " deriv: doing (M^dag M)^{-1/" << param.inv_pow << "} ( (V^dag V)^{1/" << 2*param.inv_pow << "} Phi)" << std::endl;
	multiShiftInverse(Denominator, ApproxNegPowerMD, param.MaxIter, MpvPhi,MfMpvPhi_k,MfMpvPhi);

	std::cout<<GridLogMessage << action_name() << " deriv: doing (V^dag V)^{1/" << 2*param.inv_pow << "} ( (M^dag M)^{-1/" << param.inv_pow << "} (V^dag V)^{1/" << 2*param.inv_pow << "} Phi)" << std::endl;
	multiShiftInverse(Numerator, ApproxHalfPowerMD, param.MaxIter, MfMpvPhi,MpvMfMpvPhi_k,MpvMfMpvPhi);
		

	SchurDifferentiableOperator<Impl> MdagM(DenOp);
	SchurDifferentiableOperator<Impl> VdagV(NumOp);


	RealD ak;

	dSdU = Zero();

	// With these building blocks  
	//  
	//       dS/dU = 
	//                 \sum_k -ak MfMpvPhi_k^dag      [ dM^dag M + M^dag dM ] MfMpvPhi_k         (1)
	//             +   \sum_k -ak MpvMfMpvPhi_k^\dag  [ dV^dag V + V^dag dV ] MpvPhi_k           (2)
	//                        -ak MpvPhi_k^dag        [ dV^dag V + V^dag dV ] MpvMfMpvPhi_k      (3)

	//(1)	
	std::cout<<GridLogMessage << action_name() << " deriv: doing dS/dU part (1)" << std::endl;
	for(int k=0;k<n_f;k++){
	  ak = ApproxNegPowerMD.residues[k];
	  MdagM.Mpc(MfMpvPhi_k[k],Y);
	  MdagM.MpcDagDeriv(tmp , MfMpvPhi_k[k], Y );  dSdU=dSdU+ak*tmp;
	  MdagM.MpcDeriv(tmp , Y, MfMpvPhi_k[k] );  dSdU=dSdU+ak*tmp;
	}
	
	//(2)
	//(3)
	std::cout<<GridLogMessage << action_name() << " deriv: doing dS/dU part (2)+(3)" << std::endl;
	for(int k=0;k<n_pv;k++){

          ak = ApproxHalfPowerMD.residues[k];
	  
	  VdagV.Mpc(MpvPhi_k[k],Y);
	  VdagV.MpcDagDeriv(tmp,MpvMfMpvPhi_k[k],Y); dSdU=dSdU+ak*tmp;
	  VdagV.MpcDeriv   (tmp,Y,MpvMfMpvPhi_k[k]);  dSdU=dSdU+ak*tmp;     
	  
	  VdagV.Mpc(MpvMfMpvPhi_k[k],Y);                // V as we take Ydag 
	  VdagV.MpcDeriv   (tmp,Y, MpvPhi_k[k]); dSdU=dSdU+ak*tmp;
	  VdagV.MpcDagDeriv(tmp,MpvPhi_k[k], Y); dSdU=dSdU+ak*tmp;

	}

	//dSdU = Ta(dSdU);
	std::cout<<GridLogMessage << action_name() << " deriv: complete" << std::endl;
      };
    };

NAMESPACE_END(Grid);

#endif
