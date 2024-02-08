/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/action/pseudofermion/ExactOneFlavourRatio.h

Copyright (C) 2017

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: David Murphy <dmurphy@phys.columbia.edu>

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

/////////////////////////////////////////////////////////////////
// Implementation of exact one flavour algorithm (EOFA)         //
// using fermion classes defined in:                           //
//    Grid/qcd/action/fermion/DomainWallEOFAFermion.h (Shamir) //
//    Grid/qcd/action/fermion/MobiusEOFAFermion.h (Mobius)     //
// arXiv: 1403.1683, 1706.05843                                //
/////////////////////////////////////////////////////////////////

#ifndef QCD_PSEUDOFERMION_EXACT_ONE_FLAVOUR_RATIO_H
#define QCD_PSEUDOFERMION_EXACT_ONE_FLAVOUR_RATIO_H

NAMESPACE_BEGIN(Grid);

  ///////////////////////////////////////////////////////////////
  // Exact one flavour implementation of DWF determinant ratio //
  ///////////////////////////////////////////////////////////////

  //Note: using mixed prec CG for the heatbath solver in this action class will not work
  //      because the L, R operators must have their shift coefficients updated throughout the heatbath step
  //      You will find that the heatbath solver simply won't converge.
  //      To use mixed precision here use the ExactOneFlavourRatioMixedPrecHeatbathPseudoFermionAction variant below
  template<class Impl>
  class ExactOneFlavourRatioPseudoFermionAction : public Action<typename Impl::GaugeField>
  {
    public:
      INHERIT_IMPL_TYPES(Impl);
      typedef OneFlavourRationalParams Params;
      Params param;
      MultiShiftFunction PowerNegHalf;

    private:
      bool use_heatbath_forecasting;
      AbstractEOFAFermion<Impl>& Lop; // the basic LH operator
      AbstractEOFAFermion<Impl>& Rop; // the basic RH operator
      SchurRedBlackDiagMooeeSolve<FermionField> SolverHBL;
      SchurRedBlackDiagMooeeSolve<FermionField> SolverHBR;
      SchurRedBlackDiagMooeeSolve<FermionField> SolverL;
      SchurRedBlackDiagMooeeSolve<FermionField> SolverR;
      SchurRedBlackDiagMooeeSolve<FermionField> DerivativeSolverL;
      SchurRedBlackDiagMooeeSolve<FermionField> DerivativeSolverR;
      FermionField Phi; // the pseudofermion field for this trajectory

      RealD norm2_eta; //|eta|^2 where eta is the random gaussian field used to generate the pseudofermion field
      bool initial_action; //true for the first call to S after refresh, for which the identity S = |eta|^2 holds provided the rational approx is good
    public:

      //Used in the heatbath, refresh the shift coefficients of the L (LorR=0) or R (LorR=1) operator
      virtual void heatbathRefreshShiftCoefficients(int LorR, RealD to){
	AbstractEOFAFermion<Impl>&op = LorR == 0 ? Lop : Rop;
	op.RefreshShiftCoefficients(to);
      }


      //Use the same solver for L,R in all cases
      ExactOneFlavourRatioPseudoFermionAction(AbstractEOFAFermion<Impl>& _Lop, 
					      AbstractEOFAFermion<Impl>& _Rop,
					      OperatorFunction<FermionField>& CG, 
					      Params& p, 
					      bool use_fc=false) 
	: ExactOneFlavourRatioPseudoFermionAction(_Lop,_Rop,CG,CG,CG,CG,CG,CG,p,use_fc) {};

      //Use the same solver for L,R in the heatbath but different solvers elsewhere
      ExactOneFlavourRatioPseudoFermionAction(AbstractEOFAFermion<Impl>& _Lop, 
					      AbstractEOFAFermion<Impl>& _Rop,
					      OperatorFunction<FermionField>& HeatbathCG,
					      OperatorFunction<FermionField>& ActionCGL, OperatorFunction<FermionField>& ActionCGR, 
					      OperatorFunction<FermionField>& DerivCGL , OperatorFunction<FermionField>& DerivCGR, 
					      Params& p, 
					      bool use_fc=false)
	: ExactOneFlavourRatioPseudoFermionAction(_Lop,_Rop,HeatbathCG,HeatbathCG, ActionCGL, ActionCGR, DerivCGL,DerivCGR,p,use_fc) {};

      //Use different solvers for L,R in all cases
      ExactOneFlavourRatioPseudoFermionAction(AbstractEOFAFermion<Impl>& _Lop, 
					      AbstractEOFAFermion<Impl>& _Rop,
					      OperatorFunction<FermionField>& HeatbathCGL, OperatorFunction<FermionField>& HeatbathCGR,
					      OperatorFunction<FermionField>& ActionCGL, OperatorFunction<FermionField>& ActionCGR, 
					      OperatorFunction<FermionField>& DerivCGL , OperatorFunction<FermionField>& DerivCGR, 
					      Params& p, 
					      bool use_fc=false) : 
        Lop(_Lop), 
	Rop(_Rop), 
	SolverHBL(HeatbathCGL,false,true), SolverHBR(HeatbathCGR,false,true),
	SolverL(ActionCGL, false, true), SolverR(ActionCGR, false, true), 
	DerivativeSolverL(DerivCGL, false, true), DerivativeSolverR(DerivCGR, false, true), 
	Phi(_Lop.FermionGrid()), 
	param(p), 
	use_heatbath_forecasting(use_fc),
	initial_action(false)
      {
        AlgRemez remez(param.lo, param.hi, param.precision);

        // MdagM^(+- 1/2)
        std::cout << GridLogMessage << "Generating degree " << param.degree << " for x^(-1/2)" << std::endl;
        remez.generateApprox(param.degree, 1, 2);
        PowerNegHalf.Init(remez, param.tolerance, true);
      };

      const FermionField &getPhi() const{ return Phi; }

      //Set the pseudofermion Phi field for testing or other purposes
      void setPhi(const FermionField &phi_in){
	Phi = phi_in;
      }

      virtual std::string action_name() { return "ExactOneFlavourRatioPseudoFermionAction"; }

      virtual std::string LogParameters() {
        std::stringstream sstream;
        sstream << GridLogMessage << "[" << action_name() << "] Low            :" << param.lo << std::endl;
        sstream << GridLogMessage << "[" << action_name() << "] High           :" << param.hi << std::endl;
        sstream << GridLogMessage << "[" << action_name() << "] Max iterations :" << param.MaxIter << std::endl;
        sstream << GridLogMessage << "[" << action_name() << "] Tolerance      :" << param.tolerance << std::endl;
        sstream << GridLogMessage << "[" << action_name() << "] Degree         :" << param.degree << std::endl;
        sstream << GridLogMessage << "[" << action_name() << "] Precision      :" << param.precision << std::endl;
        return sstream.str();
      }

      // Spin projection
      void spProj(const FermionField& in, FermionField& out, int sign, int Ls)
      {
        if(sign == 1){ for(int s=0; s<Ls; ++s){ axpby_ssp_pplus(out, 0.0, in, 1.0, in, s, s); } }
        else{ for(int s=0; s<Ls; ++s){ axpby_ssp_pminus(out, 0.0, in, 1.0, in, s, s); } }
      }

      virtual void refresh(const GaugeField &U, GridSerialRNG &sRNG, GridParallelRNG& pRNG) {
        // P(eta_o) = e^{- eta_o^dag eta_o}
        //
        // e^{x^2/2 sig^2} => sig^2 = 0.5.
        // 
        RealD scale = std::sqrt(0.5);

        FermionField eta    (Lop.FermionGrid());
        gaussian(pRNG,eta); eta = eta * scale;

	refresh(U,eta);
      }

      // EOFA heatbath: see Eqn. (29) of arXiv:1706.05843
      // We generate a Gaussian noise vector \eta, and then compute
      //  \Phi = M_{\rm EOFA}^{-1/2} * \eta
      // using a rational approximation to the inverse square root
      //
      // As a check of rational require \Phi^dag M_{EOFA} \Phi == eta^dag M^-1/2^dag M M^-1/2 eta = eta^dag eta
      //
     void refresh(const GaugeField &U, const FermionField &eta) {
        Lop.ImportGauge(U);
        Rop.ImportGauge(U);

        FermionField CG_src      (Lop.FermionGrid());
        FermionField CG_soln     (Lop.FermionGrid());
        FermionField Forecast_src(Lop.FermionGrid());
        std::vector<FermionField> tmp(2, Lop.FermionGrid());

        // Use chronological inverter to forecast solutions across poles
        std::vector<FermionField> prev_solns;
        if(use_heatbath_forecasting){ prev_solns.reserve(param.degree); }
        ChronoForecast<AbstractEOFAFermion<Impl>, FermionField> Forecast;

        // \Phi = ( \alpha_{0} + \sum_{k=1}^{N_{p}} \alpha_{l} * \gamma_{l} ) * \eta
        RealD N(PowerNegHalf.norm);
        for(int k=0; k<param.degree; ++k){ N += PowerNegHalf.residues[k] / ( 1.0 + PowerNegHalf.poles[k] ); }
        Phi = eta * N;

        // LH terms:
        // \Phi = \Phi + k \sum_{k=1}^{N_{p}} P_{-} \Omega_{-}^{\dagger} ( H(mf)
        //          - \gamma_{l} \Delta_{-}(mf,mb) P_{-} )^{-1} \Omega_{-} P_{-} \eta
        RealD gamma_l(0.0);
        spProj(eta, tmp[0], -1, Lop.Ls);
        Lop.Omega(tmp[0], tmp[1], -1, 0);
        G5R5(CG_src, tmp[1]);
        tmp[1] = Zero();
        for(int k=0; k<param.degree; ++k){
          gamma_l = 1.0 / ( 1.0 + PowerNegHalf.poles[k] );
          heatbathRefreshShiftCoefficients(0, -gamma_l);
          if(use_heatbath_forecasting){ // Forecast CG guess using solutions from previous poles
            Lop.Mdag(CG_src, Forecast_src);
            CG_soln = Forecast(Lop, Forecast_src, prev_solns);
            SolverHBL(Lop, CG_src, CG_soln);
            prev_solns.push_back(CG_soln);
          } else {
            CG_soln = Zero(); // Just use zero as the initial guess
	    SolverHBL(Lop, CG_src, CG_soln);
          }
          Lop.Dtilde(CG_soln, tmp[0]); // We actually solved Cayley preconditioned system: transform back
          tmp[1] = tmp[1] + ( PowerNegHalf.residues[k]*gamma_l*gamma_l*Lop.k ) * tmp[0];
        }
        Lop.Omega(tmp[1], tmp[0], -1, 1);
        spProj(tmp[0], tmp[1], -1, Lop.Ls);
        Phi = Phi + tmp[1];

        // RH terms:
        // \Phi = \Phi - k \sum_{k=1}^{N_{p}} P_{+} \Omega_{+}^{\dagger} ( H(mb)
        //          + \gamma_{l} \Delta_{+}(mf,mb) P_{+} )^{-1} \Omega_{+} P_{+} \eta
        spProj(eta, tmp[0], 1, Rop.Ls);
        Rop.Omega(tmp[0], tmp[1], 1, 0);
        G5R5(CG_src, tmp[1]);
        tmp[1] = Zero();
        if(use_heatbath_forecasting){ prev_solns.clear(); } // empirically, LH solns don't help for RH solves
        for(int k=0; k<param.degree; ++k){
          gamma_l = 1.0 / ( 1.0 + PowerNegHalf.poles[k] );
	  heatbathRefreshShiftCoefficients(1, -gamma_l*PowerNegHalf.poles[k]);
          if(use_heatbath_forecasting){
            Rop.Mdag(CG_src, Forecast_src);
            CG_soln = Forecast(Rop, Forecast_src, prev_solns);
            SolverHBR(Rop, CG_src, CG_soln);
            prev_solns.push_back(CG_soln);
          } else {
            CG_soln = Zero();
            SolverHBR(Rop, CG_src, CG_soln);
          }
          Rop.Dtilde(CG_soln, tmp[0]); // We actually solved Cayley preconditioned system: transform back
          tmp[1] = tmp[1] - ( PowerNegHalf.residues[k]*gamma_l*gamma_l*Rop.k ) * tmp[0];
        }
        Rop.Omega(tmp[1], tmp[0], 1, 1);
        spProj(tmp[0], tmp[1], 1, Rop.Ls);
        Phi = Phi + tmp[1];

        // Reset shift coefficients for energy and force evals
	heatbathRefreshShiftCoefficients(0, 0.0);
	heatbathRefreshShiftCoefficients(1, -1.0);

	//Mark that the next call to S is the first after refresh
	initial_action = true;


	// Bounds check
	RealD EtaDagEta = norm2(eta);
	norm2_eta = EtaDagEta;

	//	RealD PhiDagMPhi= norm2(eta);

      };

      void Meofa(const GaugeField& U,const FermionField &in, FermionField & out) 
      {
        Lop.ImportGauge(U);
        Rop.ImportGauge(U);

        FermionField spProj_in(Lop.FermionGrid());
        std::vector<FermionField> tmp(2, Lop.FermionGrid());
	out = in;
	
        // LH term: S = S - k <\Phi| P_{-} \Omega_{-}^{\dagger} H(mf)^{-1} \Omega_{-} P_{-} |\Phi>
        spProj(in, spProj_in, -1, Lop.Ls);
        Lop.Omega(spProj_in, tmp[0], -1, 0);
        G5R5(tmp[1], tmp[0]);
        tmp[0] = Zero();
        SolverL(Lop, tmp[1], tmp[0]);
        Lop.Dtilde(tmp[0], tmp[1]); // We actually solved Cayley preconditioned system: transform back
        Lop.Omega(tmp[1], tmp[0], -1, 1);
	spProj(tmp[0], tmp[1], -1, Lop.Ls);

	out = out -  Lop.k * tmp[1];

        // RH term: S = S + k <\Phi| P_{+} \Omega_{+}^{\dagger} ( H(mb)
        //               - \Delta_{+}(mf,mb) P_{+} )^{-1} \Omega_{+} P_{+} |\Phi>
        spProj(in, spProj_in, 1, Rop.Ls);
        Rop.Omega(spProj_in, tmp[0], 1, 0);
        G5R5(tmp[1], tmp[0]);
        tmp[0] = Zero();
        SolverR(Rop, tmp[1], tmp[0]);
        Rop.Dtilde(tmp[0], tmp[1]);
        Rop.Omega(tmp[1], tmp[0], 1, 1);
	spProj(tmp[0], tmp[1], 1, Rop.Ls);

        out = out + Rop.k * tmp[1];
      }

      //Due to the structure of EOFA, it is no more expensive to compute the inverse of Meofa
      //To ensure correctness we can simply reuse the heatbath code but use the rational approx
      //f(x) = 1/x   which corresponds to alpha_0=0,  alpha_1=1,  beta_1=0 => gamma_1=1
      void MeofaInv(const GaugeField &U, const FermionField &in, FermionField &out) {
        Lop.ImportGauge(U);
        Rop.ImportGauge(U);

        FermionField CG_src      (Lop.FermionGrid());
        FermionField CG_soln     (Lop.FermionGrid());
        std::vector<FermionField> tmp(2, Lop.FermionGrid());

        // \Phi = ( \alpha_{0} + \sum_{k=1}^{N_{p}} \alpha_{l} * \gamma_{l} ) * \eta
	// = 1 * \eta
        out = in;

        // LH terms:
        // \Phi = \Phi + k \sum_{k=1}^{N_{p}} P_{-} \Omega_{-}^{\dagger} ( H(mf)
        //          - \gamma_{l} \Delta_{-}(mf,mb) P_{-} )^{-1} \Omega_{-} P_{-} \eta
        spProj(in, tmp[0], -1, Lop.Ls);
        Lop.Omega(tmp[0], tmp[1], -1, 0);
        G5R5(CG_src, tmp[1]);
        {
          heatbathRefreshShiftCoefficients(0, -1.); //-gamma_1 = -1.

	  CG_soln = Zero(); // Just use zero as the initial guess
	  SolverHBL(Lop, CG_src, CG_soln);

          Lop.Dtilde(CG_soln, tmp[0]); // We actually solved Cayley preconditioned system: transform back
          tmp[1] = Lop.k * tmp[0];
        }
        Lop.Omega(tmp[1], tmp[0], -1, 1);
        spProj(tmp[0], tmp[1], -1, Lop.Ls);
        out = out + tmp[1];

        // RH terms:
        // \Phi = \Phi - k \sum_{k=1}^{N_{p}} P_{+} \Omega_{+}^{\dagger} ( H(mb)
        //          - \beta_l\gamma_{l} \Delta_{+}(mf,mb) P_{+} )^{-1} \Omega_{+} P_{+} \eta
        spProj(in, tmp[0], 1, Rop.Ls);
        Rop.Omega(tmp[0], tmp[1], 1, 0);
        G5R5(CG_src, tmp[1]);
        {
	  heatbathRefreshShiftCoefficients(1, 0.); //-gamma_1 * beta_1 = 0

	  CG_soln = Zero();
	  SolverHBR(Rop, CG_src, CG_soln);

          Rop.Dtilde(CG_soln, tmp[0]); // We actually solved Cayley preconditioned system: transform back
          tmp[1] = - Rop.k * tmp[0];
        }
        Rop.Omega(tmp[1], tmp[0], 1, 1);
        spProj(tmp[0], tmp[1], 1, Rop.Ls);
        out = out + tmp[1];

        // Reset shift coefficients for energy and force evals
	heatbathRefreshShiftCoefficients(0, 0.0);
	heatbathRefreshShiftCoefficients(1, -1.0);
      };




      // EOFA action: see Eqn. (10) of arXiv:1706.05843
      virtual RealD S(const GaugeField& U)
      {
        Lop.ImportGauge(U);
        Rop.ImportGauge(U);

        FermionField spProj_Phi(Lop.FermionGrid());
        std::vector<FermionField> tmp(2, Lop.FermionGrid());

        // S = <\Phi|\Phi>
        RealD action(norm2(Phi));

        // LH term: S = S - k <\Phi| P_{-} \Omega_{-}^{\dagger} H(mf)^{-1} \Omega_{-} P_{-} |\Phi>
        spProj(Phi, spProj_Phi, -1, Lop.Ls);
        Lop.Omega(spProj_Phi, tmp[0], -1, 0);
        G5R5(tmp[1], tmp[0]);
        tmp[0] = Zero();
        SolverL(Lop, tmp[1], tmp[0]);
        Lop.Dtilde(tmp[0], tmp[1]); // We actually solved Cayley preconditioned system: transform back
        Lop.Omega(tmp[1], tmp[0], -1, 1);
        action -= Lop.k * innerProduct(spProj_Phi, tmp[0]).real();

        // RH term: S = S + k <\Phi| P_{+} \Omega_{+}^{\dagger} ( H(mb)
        //               - \Delta_{+}(mf,mb) P_{+} )^{-1} \Omega_{+} P_{+} |\Phi>
        spProj(Phi, spProj_Phi, 1, Rop.Ls);
        Rop.Omega(spProj_Phi, tmp[0], 1, 0);
        G5R5(tmp[1], tmp[0]);
        tmp[0] = Zero();
        SolverR(Rop, tmp[1], tmp[0]);
        Rop.Dtilde(tmp[0], tmp[1]);
        Rop.Omega(tmp[1], tmp[0], 1, 1);
        action += Rop.k * innerProduct(spProj_Phi, tmp[0]).real();

	if(initial_action){
	  //For the first call to S after refresh,  S = |eta|^2. We can use this to ensure the rational approx is good
	  RealD diff = action - norm2_eta;

	  //S_init = eta^dag M^{-1/2} M M^{-1/2} eta
	  //S_init - eta^dag eta =  eta^dag ( M^{-1/2} M M^{-1/2} - 1 ) eta

	  //If approximate solution
	  //S_init - eta^dag eta =  eta^dag ( [M^{-1/2}+\delta M^{-1/2}] M [M^{-1/2}+\delta M^{-1/2}] - 1 ) eta
	  //               \approx  eta^dag ( \delta M^{-1/2} M^{1/2} + M^{1/2}\delta M^{-1/2} ) eta
	  // We divide out |eta|^2 to remove source scaling but the tolerance on this check should still be somewhat higher than the actual approx tolerance
	  RealD test = fabs(diff)/norm2_eta; //test the quality of the rational approx

	  std::cout << GridLogMessage << action_name() << " initial action " << action << " expect " << norm2_eta << "; diff " << diff << std::endl;
	  std::cout << GridLogMessage << action_name() << "[ eta^dag ( M^{-1/2} M M^{-1/2} - 1 ) eta ]/|eta^2| = " << test << "  expect 0 (tol " << param.BoundsCheckTol << ")" << std::endl;

	  assert( ( test < param.BoundsCheckTol ) && " Initial action check failed" );
	  initial_action = false;
	}

        return action;
      };

      // EOFA pseudofermion force: see Eqns. (34)-(36) of arXiv:1706.05843
      virtual void deriv(const GaugeField& U, GaugeField& dSdU)
      {
        Lop.ImportGauge(U);
        Rop.ImportGauge(U);

        FermionField spProj_Phi      (Lop.FermionGrid());
        FermionField Omega_spProj_Phi(Lop.FermionGrid());
        FermionField CG_src          (Lop.FermionGrid());
        FermionField Chi             (Lop.FermionGrid());
        FermionField g5_R5_Chi       (Lop.FermionGrid());

        GaugeField force(Lop.GaugeGrid());

	/////////////////////////////////////////////
	// PAB: 
	//   Optional single precision derivative ?
	/////////////////////////////////////////////

        // LH: dSdU = k \chi_{L}^{\dagger} \gamma_{5} R_{5} ( \partial_{x,\mu} D_{w} ) \chi_{L}
        //     \chi_{L} = H(mf)^{-1} \Omega_{-} P_{-} \Phi
        spProj(Phi, spProj_Phi, -1, Lop.Ls);
        Lop.Omega(spProj_Phi, Omega_spProj_Phi, -1, 0);
        G5R5(CG_src, Omega_spProj_Phi);
        spProj_Phi = Zero();
        DerivativeSolverL(Lop, CG_src, spProj_Phi);
        Lop.Dtilde(spProj_Phi, Chi);
        G5R5(g5_R5_Chi, Chi);
        Lop.MDeriv(force, g5_R5_Chi, Chi, DaggerNo);
        dSdU = -Lop.k * force;

        // RH: dSdU = dSdU - k \chi_{R}^{\dagger} \gamma_{5} R_{5} ( \partial_{x,\mu} D_{w} ) \chi_{}
        //     \chi_{R} = ( H(mb) - \Delta_{+}(mf,mb) P_{+} )^{-1} \Omega_{+} P_{+} \Phi
        spProj(Phi, spProj_Phi, 1, Rop.Ls);
        Rop.Omega(spProj_Phi, Omega_spProj_Phi, 1, 0);
        G5R5(CG_src, Omega_spProj_Phi);
        spProj_Phi = Zero();
        DerivativeSolverR(Rop, CG_src, spProj_Phi);
        Rop.Dtilde(spProj_Phi, Chi);
        G5R5(g5_R5_Chi, Chi);
        Lop.MDeriv(force, g5_R5_Chi, Chi, DaggerNo);
        dSdU = dSdU + Rop.k * force;
      };
  };

  template<class ImplD, class ImplF>
  class ExactOneFlavourRatioMixedPrecHeatbathPseudoFermionAction : public ExactOneFlavourRatioPseudoFermionAction<ImplD>{
  public:
    INHERIT_IMPL_TYPES(ImplD);
    typedef OneFlavourRationalParams Params;

  private:
    AbstractEOFAFermion<ImplF>& LopF; // the basic LH operator
    AbstractEOFAFermion<ImplF>& RopF; // the basic RH operator

  public:
    
    virtual std::string action_name() { return "ExactOneFlavourRatioMixedPrecHeatbathPseudoFermionAction"; }
    
    //Used in the heatbath, refresh the shift coefficients of the L (LorR=0) or R (LorR=1) operator
    virtual void heatbathRefreshShiftCoefficients(int LorR, RealD to){
      AbstractEOFAFermion<ImplF> &op = LorR == 0 ? LopF : RopF;
      op.RefreshShiftCoefficients(to);
      this->ExactOneFlavourRatioPseudoFermionAction<ImplD>::heatbathRefreshShiftCoefficients(LorR,to);
    }
    
    ExactOneFlavourRatioMixedPrecHeatbathPseudoFermionAction(AbstractEOFAFermion<ImplF>& _LopF, 
							     AbstractEOFAFermion<ImplF>& _RopF,
							     AbstractEOFAFermion<ImplD>& _LopD, 
							     AbstractEOFAFermion<ImplD>& _RopD,
							     OperatorFunction<FermionField>& HeatbathCGL, OperatorFunction<FermionField>& HeatbathCGR,
							     OperatorFunction<FermionField>& ActionCGL, OperatorFunction<FermionField>& ActionCGR, 
							     OperatorFunction<FermionField>& DerivCGL , OperatorFunction<FermionField>& DerivCGR, 
							     Params& p, 
							     bool use_fc=false) : 
    LopF(_LopF), RopF(_RopF), ExactOneFlavourRatioPseudoFermionAction<ImplD>(_LopD, _RopD, HeatbathCGL, HeatbathCGR, ActionCGL, ActionCGR, DerivCGL, DerivCGR, p, use_fc){}
  };


NAMESPACE_END(Grid);

#endif
