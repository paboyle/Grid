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

namespace Grid{
namespace QCD{

  ///////////////////////////////////////////////////////////////
  // Exact one flavour implementation of DWF determinant ratio //
  ///////////////////////////////////////////////////////////////

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
      SchurRedBlackDiagMooeeSolve<FermionField> Solver;
      FermionField Phi; // the pseudofermion field for this trajectory

    public:
      ExactOneFlavourRatioPseudoFermionAction(AbstractEOFAFermion<Impl>& _Lop, AbstractEOFAFermion<Impl>& _Rop,
        OperatorFunction<FermionField>& S, Params& p, bool use_fc=false) : Lop(_Lop), Rop(_Rop), Solver(S),
        Phi(_Lop.FermionGrid()), param(p), use_heatbath_forecasting(use_fc)
      {
        AlgRemez remez(param.lo, param.hi, param.precision);

        // MdagM^(+- 1/2)
        std::cout << GridLogMessage << "Generating degree " << param.degree << " for x^(-1/2)" << std::endl;
        remez.generateApprox(param.degree, 1, 2);
        PowerNegHalf.Init(remez, param.tolerance, true);
      };

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

      // EOFA heatbath: see Eqn. (29) of arXiv:1706.05843
      // We generate a Gaussian noise vector \eta, and then compute
      //  \Phi = M_{\rm EOFA}^{-1/2} * \eta
      // using a rational approximation to the inverse square root
      virtual void refresh(const GaugeField& U, GridParallelRNG& pRNG)
      {
        Lop.ImportGauge(U);
        Rop.ImportGauge(U);

        FermionField eta         (Lop.FermionGrid());
        FermionField CG_src      (Lop.FermionGrid());
        FermionField CG_soln     (Lop.FermionGrid());
        FermionField Forecast_src(Lop.FermionGrid());
        std::vector<FermionField> tmp(2, Lop.FermionGrid());

        // Use chronological inverter to forecast solutions across poles
        std::vector<FermionField> prev_solns;
        if(use_heatbath_forecasting){ prev_solns.reserve(param.degree); }
        ChronoForecast<AbstractEOFAFermion<Impl>, FermionField> Forecast;

        // Seed with Gaussian noise vector (var = 0.5)
        RealD scale = std::sqrt(0.5);
        gaussian(pRNG,eta);
        eta = eta * scale;
        printf("Heatbath source vector: <\\eta|\\eta> = %1.15e\n", norm2(eta));

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
        tmp[1] = zero;
        for(int k=0; k<param.degree; ++k){
          gamma_l = 1.0 / ( 1.0 + PowerNegHalf.poles[k] );
          Lop.RefreshShiftCoefficients(-gamma_l);
          if(use_heatbath_forecasting){ // Forecast CG guess using solutions from previous poles
            Lop.Mdag(CG_src, Forecast_src);
            CG_soln = Forecast(Lop, Forecast_src, prev_solns);
            Solver(Lop, CG_src, CG_soln);
            prev_solns.push_back(CG_soln);
          } else {
            CG_soln = zero; // Just use zero as the initial guess
            Solver(Lop, CG_src, CG_soln);
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
        tmp[1] = zero;
        if(use_heatbath_forecasting){ prev_solns.clear(); } // empirically, LH solns don't help for RH solves
        for(int k=0; k<param.degree; ++k){
          gamma_l = 1.0 / ( 1.0 + PowerNegHalf.poles[k] );
          Rop.RefreshShiftCoefficients(-gamma_l*PowerNegHalf.poles[k]);
          if(use_heatbath_forecasting){
            Rop.Mdag(CG_src, Forecast_src);
            CG_soln = Forecast(Rop, Forecast_src, prev_solns);
            Solver(Rop, CG_src, CG_soln);
            prev_solns.push_back(CG_soln);
          } else {
            CG_soln = zero;
            Solver(Rop, CG_src, CG_soln);
          }
          Rop.Dtilde(CG_soln, tmp[0]); // We actually solved Cayley preconditioned system: transform back
          tmp[1] = tmp[1] - ( PowerNegHalf.residues[k]*gamma_l*gamma_l*Rop.k ) * tmp[0];
        }
        Rop.Omega(tmp[1], tmp[0], 1, 1);
        spProj(tmp[0], tmp[1], 1, Rop.Ls);
        Phi = Phi + tmp[1];

        // Reset shift coefficients for energy and force evals
        Lop.RefreshShiftCoefficients(0.0);
        Rop.RefreshShiftCoefficients(-1.0);
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
        tmp[0] = zero;
        Solver(Lop, tmp[1], tmp[0]);
        Lop.Dtilde(tmp[0], tmp[1]); // We actually solved Cayley preconditioned system: transform back
        Lop.Omega(tmp[1], tmp[0], -1, 1);
        action -= Lop.k * innerProduct(spProj_Phi, tmp[0]).real();

        // RH term: S = S + k <\Phi| P_{+} \Omega_{+}^{\dagger} ( H(mb)
        //               - \Delta_{+}(mf,mb) P_{+} )^{-1} \Omega_{-} P_{-} |\Phi>
        spProj(Phi, spProj_Phi, 1, Rop.Ls);
        Rop.Omega(spProj_Phi, tmp[0], 1, 0);
        G5R5(tmp[1], tmp[0]);
        tmp[0] = zero;
        Solver(Rop, tmp[1], tmp[0]);
        Rop.Dtilde(tmp[0], tmp[1]);
        Rop.Omega(tmp[1], tmp[0], 1, 1);
        action += Rop.k * innerProduct(spProj_Phi, tmp[0]).real();

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

        // LH: dSdU = k \chi_{L}^{\dagger} \gamma_{5} R_{5} ( \partial_{x,\mu} D_{w} ) \chi_{L}
        //     \chi_{L} = H(mf)^{-1} \Omega_{-} P_{-} \Phi
        spProj(Phi, spProj_Phi, -1, Lop.Ls);
        Lop.Omega(spProj_Phi, Omega_spProj_Phi, -1, 0);
        G5R5(CG_src, Omega_spProj_Phi);
        spProj_Phi = zero;
        Solver(Lop, CG_src, spProj_Phi);
        Lop.Dtilde(spProj_Phi, Chi);
        G5R5(g5_R5_Chi, Chi);
        Lop.MDeriv(force, g5_R5_Chi, Chi, DaggerNo);
        dSdU = Lop.k * force;

        // RH: dSdU = dSdU - k \chi_{R}^{\dagger} \gamma_{5} R_{5} ( \partial_{x,\mu} D_{w} ) \chi_{}
        //     \chi_{R} = ( H(mb) - \Delta_{+}(mf,mb) P_{+} )^{-1} \Omega_{+} P_{+} \Phi
        spProj(Phi, spProj_Phi, 1, Rop.Ls);
        Rop.Omega(spProj_Phi, Omega_spProj_Phi, 1, 0);
        G5R5(CG_src, Omega_spProj_Phi);
        spProj_Phi = zero;
        Solver(Rop, CG_src, spProj_Phi);
        Rop.Dtilde(spProj_Phi, Chi);
        G5R5(g5_R5_Chi, Chi);
        Lop.MDeriv(force, g5_R5_Chi, Chi, DaggerNo);
        dSdU = dSdU - Rop.k * force;
      };
  };
}}

#endif
