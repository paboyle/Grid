    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/action/pseudofermion/OneFlavourEvenOddRationalRatio.h

    Copyright (C) 2015

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
#ifndef QCD_PSEUDOFERMION_ONE_FLAVOUR_EVEN_ODD_RATIONAL_RATIO_H
#define QCD_PSEUDOFERMION_ONE_FLAVOUR_EVEN_ODD_RATIONAL_RATIO_H

namespace Grid{
  namespace QCD{

    ///////////////////////////////////////
    // One flavour rational
    ///////////////////////////////////////

    // S_f = chi^dag* P(V^dag*V)/Q(V^dag*V)* N(M^dag*M)/D(M^dag*M)* P(V^dag*V)/Q(V^dag*V)* chi       
    //
    // Here P/Q \sim R_{1/4}  ~ (V^dagV)^{1/4}  
    // Here N/D \sim R_{-1/2} ~ (M^dagM)^{-1/2}  
  
    template<class Impl>
    class OneFlavourEvenOddRatioRationalPseudoFermionAction : public Action<typename Impl::GaugeField> {
    public:

      INHERIT_IMPL_TYPES(Impl);

      typedef OneFlavourRationalParams Params;
      Params param;

      MultiShiftFunction PowerHalf   ;
      MultiShiftFunction PowerNegHalf;
      MultiShiftFunction PowerQuarter;
      MultiShiftFunction PowerNegQuarter;

    private:
     
      FermionOperator<Impl> & NumOp;// the basic operator
      FermionOperator<Impl> & DenOp;// the basic operator
      FermionField PhiEven; // the pseudo fermion field for this trajectory
      FermionField PhiOdd; // the pseudo fermion field for this trajectory

    public:

      OneFlavourEvenOddRatioRationalPseudoFermionAction(FermionOperator<Impl>  &_NumOp, 
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

	// MdagM^(+- 1/2)
	std::cout<<GridLogMessage << "Generating degree "<<param.degree<<" for x^(1/2)"<<std::endl;
	remez.generateApprox(param.degree,1,2);
	PowerHalf.Init(remez,param.tolerance,false);
	PowerNegHalf.Init(remez,param.tolerance,true);

	// MdagM^(+- 1/4)
	std::cout<<GridLogMessage << "Generating degree "<<param.degree<<" for x^(1/4)"<<std::endl;
	remez.generateApprox(param.degree,1,4);
   	PowerQuarter.Init(remez,param.tolerance,false);
	PowerNegQuarter.Init(remez,param.tolerance,true);
      };
      
      virtual void refresh(const GaugeField &U, GridParallelRNG& pRNG) {

	// S_f = chi^dag* P(V^dag*V)/Q(V^dag*V)* N(M^dag*M)/D(M^dag*M)* P(V^dag*V)/Q(V^dag*V)* chi       
	//
	// P(phi) = e^{- phi^dag (VdagV)^1/4 (MdagM)^-1/2 (VdagV)^1/4 phi}
	//        = e^{- phi^dag  (VdagV)^1/4 (MdagM)^-1/4 (MdagM)^-1/4  (VdagV)^1/4 phi}
	//
	// Phi =  (VdagV)^-1/4 Mdag^{1/4} eta 
	//
	// P(eta) = e^{- eta^dag eta}
	//
	// e^{x^2/2 sig^2} => sig^2 = 0.5.
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


	// MdagM^1/4 eta
	SchurDifferentiableOperator<Impl> MdagM(DenOp);
	ConjugateGradientMultiShift<FermionField> msCG_M(param.MaxIter,PowerQuarter);
	msCG_M(MdagM,etaOdd,tmp);

	// VdagV^-1/4 MdagM^1/4 eta
	SchurDifferentiableOperator<Impl> VdagV(NumOp);
	ConjugateGradientMultiShift<FermionField> msCG_V(param.MaxIter,PowerNegQuarter);
	msCG_V(VdagV,tmp,PhiOdd);

	assert(NumOp.ConstEE() == 1);
	assert(DenOp.ConstEE() == 1);
	PhiEven = zero;
	
      };

      //////////////////////////////////////////////////////
      // S_f = chi^dag* P(V^dag*V)/Q(V^dag*V)* N(M^dag*M)/D(M^dag*M)* P(V^dag*V)/Q(V^dag*V)* chi       
      //////////////////////////////////////////////////////
      virtual RealD S(const GaugeField &U) {

	NumOp.ImportGauge(U);
	DenOp.ImportGauge(U);

	FermionField X(NumOp.FermionRedBlackGrid());
	FermionField Y(NumOp.FermionRedBlackGrid());

	// VdagV^1/4 Phi
	SchurDifferentiableOperator<Impl> VdagV(NumOp);
	ConjugateGradientMultiShift<FermionField> msCG_V(param.MaxIter,PowerQuarter);
	msCG_V(VdagV,PhiOdd,X);

	// MdagM^-1/4 VdagV^1/4 Phi
	SchurDifferentiableOperator<Impl> MdagM(DenOp);
	ConjugateGradientMultiShift<FermionField> msCG_M(param.MaxIter,PowerNegQuarter);
	msCG_M(MdagM,X,Y);

	//  Phidag VdagV^1/4 MdagM^-1/4  MdagM^-1/4 VdagV^1/4 Phi
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

	const int n_f  = PowerNegHalf.poles.size();
	const int n_pv = PowerQuarter.poles.size();

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

	ConjugateGradientMultiShift<FermionField> msCG_V(param.MaxIter,PowerQuarter);
	ConjugateGradientMultiShift<FermionField> msCG_M(param.MaxIter,PowerNegHalf);

	msCG_V(VdagV,PhiOdd,MpvPhi_k,MpvPhi);
	msCG_M(MdagM,MpvPhi,MfMpvPhi_k,MfMpvPhi);
	msCG_V(VdagV,MfMpvPhi,MpvMfMpvPhi_k,MpvMfMpvPhi);

	RealD ak;

	dSdU = zero;

	// With these building blocks  
	//  
	//       dS/dU = 
	//                 \sum_k -ak MfMpvPhi_k^dag      [ dM^dag M + M^dag dM ] MfMpvPhi_k         (1)
	//             +   \sum_k -ak MpvMfMpvPhi_k^\dag  [ dV^dag V + V^dag dV ] MpvPhi_k           (2)
	//                        -ak MpvPhi_k^dag        [ dV^dag V + V^dag dV ] MpvMfMpvPhi_k      (3)

	//(1)
	for(int k=0;k<n_f;k++){
	  ak = PowerNegHalf.residues[k];
	  MdagM.Mpc(MfMpvPhi_k[k],Y);
	  MdagM.MpcDagDeriv(tmp , MfMpvPhi_k[k], Y );  dSdU=dSdU+ak*tmp;
	  MdagM.MpcDeriv(tmp , Y, MfMpvPhi_k[k] );  dSdU=dSdU+ak*tmp;
	}
	
	//(2)
	//(3)
	for(int k=0;k<n_pv;k++){

          ak = PowerQuarter.residues[k];
	  
	  VdagV.Mpc(MpvPhi_k[k],Y);
	  VdagV.MpcDagDeriv(tmp,MpvMfMpvPhi_k[k],Y); dSdU=dSdU+ak*tmp;
	  VdagV.MpcDeriv   (tmp,Y,MpvMfMpvPhi_k[k]);  dSdU=dSdU+ak*tmp;     
	  
	  VdagV.Mpc(MpvMfMpvPhi_k[k],Y);                // V as we take Ydag 
	  VdagV.MpcDeriv   (tmp,Y, MpvPhi_k[k]); dSdU=dSdU+ak*tmp;
	  VdagV.MpcDagDeriv(tmp,MpvPhi_k[k], Y); dSdU=dSdU+ak*tmp;

	}

	dSdU = Ta(dSdU);

      };
    };
  }
}


#endif
