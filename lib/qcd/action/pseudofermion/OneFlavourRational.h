    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/action/pseudofermion/OneFlavourRational.h

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
#ifndef QCD_PSEUDOFERMION_ONE_FLAVOUR_RATIONAL_H
#define QCD_PSEUDOFERMION_ONE_FLAVOUR_RATIONAL_H

namespace Grid{
  namespace QCD{

    ///////////////////////////////////////
    // One flavour rational
    ///////////////////////////////////////

    // S_f = chi^dag *  N(M^dag*M)/D(M^dag*M) * chi
    //
    // Here, M is some operator 
    // N and D makeup the rat. poly 
    //
  
    template<class Impl>
    class OneFlavourRationalPseudoFermionAction : public Action<typename Impl::GaugeField> {
    public:
      INHERIT_IMPL_TYPES(Impl);

      typedef OneFlavourRationalParams Params;
      Params param;

      MultiShiftFunction PowerHalf   ;
      MultiShiftFunction PowerNegHalf;
      MultiShiftFunction PowerQuarter;
      MultiShiftFunction PowerNegQuarter;

    private:
     
      FermionOperator<Impl> & FermOp;// the basic operator

      // NOT using "Nroots"; IroIro is -- perhaps later, but this wasn't good for us historically
      // and hasenbusch works better

      FermionField Phi; // the pseudo fermion field for this trajectory

    public:

      OneFlavourRationalPseudoFermionAction(FermionOperator<Impl>  &Op, 
					    Params & p
					    ) : FermOp(Op), Phi(Op.FermionGrid()), param(p) 
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

	// P(phi) = e^{- phi^dag (MdagM)^-1/2 phi}
	//        = e^{- phi^dag (MdagM)^-1/4 (MdagM)^-1/4 phi}
	// Phi = Mdag^{1/4} eta 
	// P(eta) = e^{- eta^dag eta}
	//
	// e^{x^2/2 sig^2} => sig^2 = 0.5.
	// 
	// So eta should be of width sig = 1/sqrt(2).

	RealD scale = std::sqrt(0.5);

	FermionField eta(FermOp.FermionGrid());

	gaussian(pRNG,eta);

	FermOp.ImportGauge(U);

	// mutishift CG
	MdagMLinearOperator<FermionOperator<Impl> ,FermionField> MdagMOp(FermOp);
	ConjugateGradientMultiShift<FermionField> msCG(param.MaxIter,PowerQuarter);
	msCG(MdagMOp,eta,Phi);

	Phi=Phi*scale;
	
      };

      //////////////////////////////////////////////////////
      // S = phi^dag (Mdag M)^-1/2 phi
      //////////////////////////////////////////////////////
      virtual RealD S(const GaugeField &U) {

	FermOp.ImportGauge(U);

	FermionField Y(FermOp.FermionGrid());
	
	MdagMLinearOperator<FermionOperator<Impl> ,FermionField> MdagMOp(FermOp);

	ConjugateGradientMultiShift<FermionField> msCG(param.MaxIter,PowerNegQuarter);

	msCG(MdagMOp,Phi,Y);

	RealD action = norm2(Y);
	std::cout << GridLogMessage << "Pseudofermion action FIXME -- is -1/4 solve or -1/2 solve faster??? "<<action<<std::endl;
	return action;
      };

      //////////////////////////////////////////////////////
      // Need
      // dS_f/dU = chi^dag   d[N/D]  chi
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
      //       Mf Phi_k = [MdagM+bk]^{-1} Phi
      //       Mf Phi   = \sum_k ak [MdagM+bk]^{-1} Phi
      //
      // With these building blocks
      //
      //       dS/dU =  \sum_k -ak Mf Phi_k^dag      [ dM^dag M + M^dag dM ] Mf Phi_k
      //        S    = innerprodReal(Phi,Mf Phi);
      //////////////////////////////////////////////////////
      virtual void deriv(const GaugeField &U,GaugeField & dSdU) {

	const int Npole = PowerNegHalf.poles.size();

	std::vector<FermionField> MPhi_k (Npole,FermOp.FermionGrid());

	FermionField X(FermOp.FermionGrid());
	FermionField Y(FermOp.FermionGrid());

	GaugeField   tmp(FermOp.GaugeGrid());

	FermOp.ImportGauge(U);

	MdagMLinearOperator<FermionOperator<Impl> ,FermionField> MdagMOp(FermOp);

	ConjugateGradientMultiShift<FermionField> msCG(param.MaxIter,PowerNegHalf);

	msCG(MdagMOp,Phi,MPhi_k);

	dSdU = zero;
	for(int k=0;k<Npole;k++){

	  RealD ak = PowerNegHalf.residues[k];

	  X  = MPhi_k[k];

	  FermOp.M(X,Y);

	  FermOp.MDeriv(tmp , Y, X,DaggerNo );  dSdU=dSdU+ak*tmp;
	  FermOp.MDeriv(tmp , X, Y,DaggerYes);  dSdU=dSdU+ak*tmp;

	}

	dSdU = Ta(dSdU);

      };
    };
  }
}


#endif
