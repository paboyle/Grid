/*************************************************************************************
 
 Grid physics library, www.github.com/paboyle/Grid
 
 Source file: ./lib/qcd/utils/BaryonUtils.h
 
 Copyright (C) 2019
 
 Author: Felix Erben <felix.erben@ed.ac.uk>
 Author: Raoul Hodgson <raoul.hodgson@ed.ac.uk>

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
#pragma once
#include <Grid/Eigen/unsupported/CXX11/Tensor>

NAMESPACE_BEGIN(Grid);

template <typename FImpl>
class BaryonUtils 
{
public:
  typedef typename FImpl::ComplexField ComplexField;
  typedef typename FImpl::FermionField FermionField;
  typedef typename FImpl::PropagatorField PropagatorField;

  typedef Lattice<iSpinMatrix<typename FImpl::Simd>> SpinMatrixField;

  private: 
  template <class mobj, class robj> accelerator_inline
  static void BaryonSite(const mobj &D1,
         const mobj &D2,
         const mobj &D3,
         const Gamma GammaA_left,
         const Gamma GammaB_left,
         const Gamma GammaA_right,
         const Gamma GammaB_right,
         const int parity,
         const int wick_contractions,
         robj &result);
  template <class mobj, class robj> accelerator_inline
  static void BaryonSiteMatrix(const mobj &D1,
         const mobj &D2,
         const mobj &D3,
         const Gamma GammaA_left,
         const Gamma GammaB_left,
         const Gamma GammaA_right,
         const Gamma GammaB_right,
         const int wick_contractions,
         robj &result);
  public:
  static void WickContractions(std::string qi, 
                 std::string qf, 
                 int &wick_contractions);
  static void ContractBaryons(const PropagatorField &q1_left,
         const PropagatorField &q2_left,
         const PropagatorField &q3_left,
         const Gamma GammaA_left,
         const Gamma GammaB_left,
         const Gamma GammaA_right,
         const Gamma GammaB_right,
         const int wick_contractions,
         const int parity,
         ComplexField &baryon_corr);
  static void ContractBaryonsMatrix(const PropagatorField &q1_left,
         const PropagatorField &q2_left,
         const PropagatorField &q3_left,
         const Gamma GammaA_left,
         const Gamma GammaB_left,
         const Gamma GammaA_right,
         const Gamma GammaB_right,
         const int wick_contractions,
         SpinMatrixField &baryon_corr);
  template <class mobj, class robj>
  static void ContractBaryonsSliced(const mobj &D1,
         const mobj &D2,
         const mobj &D3,
         const Gamma GammaA_left,
         const Gamma GammaB_left,
         const Gamma GammaA_right,
         const Gamma GammaB_right,
         const int wick_contractions,
         const int parity,
         const int nt,
         robj &result);
  template <class mobj, class robj>
  static void ContractBaryonsSlicedMatrix(const mobj &D1,
         const mobj &D2,
         const mobj &D3,
         const Gamma GammaA_left,
         const Gamma GammaB_left,
         const Gamma GammaA_right,
         const Gamma GammaB_right,
         const int wick_contractions,
         const int nt,
         robj &result);
  private:
  template <class mobj, class mobj2, class robj> accelerator_inline
  static void BaryonGamma3ptGroup1Site(
           const mobj &Dq1_ti,
           const mobj2 &Dq2_spec,
           const mobj2 &Dq3_spec,
           const mobj &Dq4_tf,
           const Gamma GammaJ,
           const Gamma GammaBi,
           const Gamma GammaBf,
           int wick_contraction,
           robj &result);

  template <class mobj, class mobj2, class robj> accelerator_inline
  static void BaryonGamma3ptGroup2Site(
           const mobj2 &Dq1_spec,
           const mobj &Dq2_ti,
           const mobj2 &Dq3_spec,
           const mobj &Dq4_tf,
           const Gamma GammaJ,
           const Gamma GammaBi,
           const Gamma GammaBf,
           int wick_contraction,
           robj &result);

  template <class mobj, class mobj2, class robj> accelerator_inline
  static void BaryonGamma3ptGroup3Site(
           const mobj2 &Dq1_spec,
           const mobj2 &Dq2_spec,
           const mobj &Dq3_ti,
           const mobj &Dq4_tf,
           const Gamma GammaJ,
           const Gamma GammaBi,
           const Gamma GammaBf,
           int wick_contraction,
           robj &result);
  public:
  template <class mobj>
  static void BaryonGamma3pt(
           const PropagatorField &q_ti,
           const mobj &Dq_spec1,
           const mobj &Dq_spec2,
           const PropagatorField &q_tf,
           int group,
           int wick_contraction,
           const Gamma GammaJ,
           const Gamma GammaBi,
           const Gamma GammaBf,
           SpinMatrixField &stn_corr);
  private: 
  template <class mobj, class mobj2, class robj> accelerator_inline
  static void SigmaToNucleonQ1EyeSite(const mobj &Dq_loop,
             const mobj2 &Du_spec,
             const mobj &Dd_tf,
             const mobj &Ds_ti,
             const Gamma Gamma_H,
             const Gamma GammaB_sigma,
             const Gamma GammaB_nucl,
             robj &result);
  template <class mobj, class mobj2, class robj> accelerator_inline
  static void SigmaToNucleonQ1NonEyeSite(const mobj &Du_ti,
             const mobj &Du_tf,
             const mobj2 &Du_spec,
             const mobj &Dd_tf,
             const mobj &Ds_ti,
             const Gamma Gamma_H,
             const Gamma GammaB_sigma,
             const Gamma GammaB_nucl,
             robj &result);


  template <class mobj, class mobj2, class robj> accelerator_inline
  static void SigmaToNucleonQ2EyeSite(const mobj &Dq_loop,
             const mobj2 &Du_spec,
             const mobj &Dd_tf,
             const mobj &Ds_ti,
             const Gamma Gamma_H,
             const Gamma GammaB_sigma,
             const Gamma GammaB_nucl,
             robj &result);
  template <class mobj, class mobj2, class robj> accelerator_inline
  static void SigmaToNucleonQ2NonEyeSite(const mobj &Du_ti,
	     const mobj &Du_tf,
	     const mobj2 &Du_spec,
	     const mobj &Dd_tf,
	     const mobj &Ds_ti,
	     const Gamma Gamma_H,
	     const Gamma GammaB_sigma,
	     const Gamma GammaB_nucl,
	     robj &result);
  template <class mobj, class mobj2, class robj> accelerator_inline
  static void XiToSigmaQ1EyeSite(const mobj &Dq_loop,
	     const mobj2 &Dd_spec,
	     const mobj2 &Ds_spec,
	     const mobj &Dd_tf,
	     const mobj &Ds_ti,
	     const Gamma Gamma_H,
	     const Gamma GammaB_sigma,
	     const Gamma GammaB_nucl,
	     robj &result);
  template <class mobj, class mobj2, class robj> accelerator_inline
  static void XiToSigmaQ2EyeSite(const mobj &Dq_loop,
	     const mobj2 &Dd_spec,
	     const mobj2 &Ds_spec,
	     const mobj &Dd_tf,
	     const mobj &Ds_ti,
	     const Gamma Gamma_H,
	     const Gamma GammaB_sigma,
	     const Gamma GammaB_nucl,
	     robj &result);
  public:
  template <class mobj>
  static void SigmaToNucleonEye(const PropagatorField &qq_loop,
         const mobj &Du_spec,
         const PropagatorField &qd_tf,
         const PropagatorField &qs_ti,
         const Gamma Gamma_H,
         const Gamma GammaB_sigma,
         const Gamma GammaB_nucl,
         const std::string op,
         SpinMatrixField &stn_corr);
  template <class mobj>
  static void SigmaToNucleonNonEye(const PropagatorField &qq_ti,
	 const PropagatorField &qq_tf,
	 const mobj &Du_spec,
	 const PropagatorField &qd_tf,
	 const PropagatorField &qs_ti,
	 const Gamma Gamma_H,
	 const Gamma GammaB_sigma,
	 const Gamma GammaB_nucl,
         const std::string op,
	 SpinMatrixField &stn_corr);
  template <class mobj>
  static void XiToSigmaEye(const PropagatorField &qq_loop,
	 const mobj &Dd_spec,
	 const mobj &Ds_spec,
	 const PropagatorField &qd_tf,
	 const PropagatorField &qs_ti,
	 const Gamma Gamma_H,
	 const Gamma GammaB_sigma,
	 const Gamma GammaB_nucl,
         const std::string op,
	 SpinMatrixField &xts_corr);
};
//This computes a baryon contraction on a lattice site, including the spin-trace of the correlation matrix
template <class FImpl>
template <class mobj, class robj> accelerator_inline
void BaryonUtils<FImpl>::BaryonSite(const mobj &D1,
                const mobj &D2,
                const mobj &D3,
                const Gamma GammaA_i,
                const Gamma GammaB_i,
                const Gamma GammaA_f,
                const Gamma GammaB_f,
                const int parity,
                const int wick_contraction,
                robj &result)
{

  Gamma g4(Gamma::Algebra::GammaT); //needed for parity P_\pm = 0.5*(1 \pm \gamma_4)
    
  auto D1_GAi =  D1 * GammaA_i;
  auto D1_GAi_g4 = D1_GAi * g4;
  auto D1_GAi_P = 0.5*(D1_GAi + (Real)parity * D1_GAi_g4);
  auto GAf_D1_GAi_P = GammaA_f * D1_GAi_P;
  auto GBf_D1_GAi_P = GammaB_f * D1_GAi_P;

  auto D2_GBi = D2 * GammaB_i;
  auto GBf_D2_GBi = GammaB_f * D2_GBi;
  auto GAf_D2_GBi = GammaA_f * D2_GBi;

  auto GBf_D3 = GammaB_f * D3;
  auto GAf_D3 = GammaA_f * D3;

  Real ee;

  for (int ie_f=0; ie_f < 6 ; ie_f++){
    int a_f    = (ie_f < 3 ? ie_f       : (6-ie_f)%3 ); //epsilon[ie_n][0]; //a
    int b_f    = (ie_f < 3 ? (ie_f+1)%3 : (8-ie_f)%3 ); //epsilon[ie_n][1]; //b
    int c_f    = (ie_f < 3 ? (ie_f+2)%3 : (7-ie_f)%3 ); //epsilon[ie_n][2]; //c
    int eSgn_f = (ie_f < 3 ? 1 : -1);
    for (int ie_i=0; ie_i < 6 ; ie_i++){
      int a_i = (ie_i < 3 ? ie_i       : (6-ie_i)%3 ); //epsilon[ie_s][0]; //a'
      int b_i = (ie_i < 3 ? (ie_i+1)%3 : (8-ie_i)%3 ); //epsilon[ie_s][1]; //b'
      int c_i = (ie_i < 3 ? (ie_i+2)%3 : (7-ie_i)%3 ); //epsilon[ie_s][2]; //c'
      int eSgn_i = (ie_i < 3 ? 1 : -1);

      ee = Real(eSgn_f * eSgn_i); //epsilon_sgn[ie_n] * epsilon_sgn[ie_s];
      //This is the \delta_{456}^{123} part
      if (wick_contraction & 1){
        for (int rho=0; rho<Ns; rho++){
          auto GAf_D1_GAi_P_rr_cc = GAf_D1_GAi_P()(rho,rho)(c_f,c_i);
          for (int alpha_f=0; alpha_f<Ns; alpha_f++){
          for (int beta_i=0; beta_i<Ns; beta_i++){
            result()()() += ee  * GAf_D1_GAi_P_rr_cc
                                        * D2_GBi    ()(alpha_f,beta_i)(a_f,a_i)
                                        * GBf_D3    ()(alpha_f,beta_i)(b_f,b_i);
          }}
        }
      }   
      //This is the \delta_{456}^{231} part
      if (wick_contraction & 2){
        for (int rho=0; rho<Ns; rho++){
        for (int alpha_f=0; alpha_f<Ns; alpha_f++){
          auto D1_GAi_P_ar_ac = D1_GAi_P()(alpha_f,rho)(a_f,c_i);
          for (int beta_i=0; beta_i<Ns; beta_i++){
            result()()() += ee  * D1_GAi_P_ar_ac
                                        * GBf_D2_GBi    ()(alpha_f,beta_i)(b_f,a_i)
                                        * GAf_D3        ()(rho,beta_i)(c_f,b_i);
          }
        }}
      }   
      //This is the \delta_{456}^{312} part
      if (wick_contraction & 4){
        for (int rho=0; rho<Ns; rho++){
        for (int alpha_f=0; alpha_f<Ns; alpha_f++){
          auto GBf_D1_GAi_P_ar_bc = GBf_D1_GAi_P()(alpha_f,rho)(b_f,c_i);
          for (int beta_i=0; beta_i<Ns; beta_i++){
            result()()() += ee  * GBf_D1_GAi_P_ar_bc
                                        * GAf_D2_GBi    ()(rho,beta_i)(c_f,a_i)
                                        * D3            ()(alpha_f,beta_i)(a_f,b_i);
          }
        }}
      }   
      //This is the \delta_{456}^{132} part
      if (wick_contraction & 8){
        for (int rho=0; rho<Ns; rho++){
          auto GAf_D1_GAi_P_rr_cc = GAf_D1_GAi_P()(rho,rho)(c_f,c_i);
          for (int alpha_f=0; alpha_f<Ns; alpha_f++){
          for (int beta_i=0; beta_i<Ns; beta_i++){
            result()()() -= ee  * GAf_D1_GAi_P_rr_cc
                                        * GBf_D2_GBi    ()(alpha_f,beta_i)(b_f,a_i)
                                        * D3            ()(alpha_f,beta_i)(a_f,b_i);
          }}
        }
      }   
      //This is the \delta_{456}^{321} part
      if (wick_contraction & 16){
        for (int rho=0; rho<Ns; rho++){
        for (int alpha_f=0; alpha_f<Ns; alpha_f++){
          auto GBf_D1_GAi_P_ar_bc = GBf_D1_GAi_P()(alpha_f,rho)(b_f,c_i);
          for (int beta_i=0; beta_i<Ns; beta_i++){
            result()()() -= ee  * GBf_D1_GAi_P_ar_bc
                                        * D2_GBi    ()(alpha_f,beta_i)(a_f,a_i)
                                        * GAf_D3    ()(rho,beta_i)(c_f,b_i);
          }
        }}
      }   
      //This is the \delta_{456}^{213} part
      if (wick_contraction & 32){
        for (int rho=0; rho<Ns; rho++){
        for (int alpha_f=0; alpha_f<Ns; alpha_f++){
          auto D1_GAi_P_ar_ac = D1_GAi_P()(alpha_f,rho)(a_f,c_i);
          for (int beta_i=0; beta_i<Ns; beta_i++){
            result()()() -= ee  * D1_GAi_P_ar_ac
                                        * GAf_D2_GBi    ()(rho,beta_i)(c_f,a_i)
                                        * GBf_D3        ()(alpha_f,beta_i)(b_f,b_i);
          }
        }}
      }
    }
  }
}

//New version without parity projection or trace
template <class FImpl>
template <class mobj, class robj> accelerator_inline
void BaryonUtils<FImpl>::BaryonSiteMatrix(const mobj &D1,
                const mobj &D2,
                const mobj &D3,
                const Gamma GammaA_i,
                const Gamma GammaB_i,
                const Gamma GammaA_f,
                const Gamma GammaB_f,
                const int wick_contraction,
                robj &result)
{

  auto D1_GAi =  D1 * GammaA_i;
  auto GAf_D1_GAi = GammaA_f * D1_GAi;
  auto GBf_D1_GAi = GammaB_f * D1_GAi;

  auto D2_GBi = D2 * GammaB_i;
  auto GBf_D2_GBi = GammaB_f * D2_GBi;
  auto GAf_D2_GBi = GammaA_f * D2_GBi;

  auto GBf_D3 = GammaB_f * D3;
  auto GAf_D3 = GammaA_f * D3;

  Real ee;

  for (int ie_f=0; ie_f < 6 ; ie_f++){
    int a_f    = (ie_f < 3 ? ie_f       : (6-ie_f)%3 ); //epsilon[ie_n][0]; //a
    int b_f    = (ie_f < 3 ? (ie_f+1)%3 : (8-ie_f)%3 ); //epsilon[ie_n][1]; //b
    int c_f    = (ie_f < 3 ? (ie_f+2)%3 : (7-ie_f)%3 ); //epsilon[ie_n][2]; //c
    int eSgn_f = (ie_f < 3 ? 1 : -1);
    for (int ie_i=0; ie_i < 6 ; ie_i++){
      int a_i = (ie_i < 3 ? ie_i       : (6-ie_i)%3 ); //epsilon[ie_s][0]; //a'
      int b_i = (ie_i < 3 ? (ie_i+1)%3 : (8-ie_i)%3 ); //epsilon[ie_s][1]; //b'
      int c_i = (ie_i < 3 ? (ie_i+2)%3 : (7-ie_i)%3 ); //epsilon[ie_s][2]; //c'
      int eSgn_i = (ie_i < 3 ? 1 : -1);

      ee = Real(eSgn_f * eSgn_i); //epsilon_sgn[ie_n] * epsilon_sgn[ie_s];
      //This is the \delta_{456}^{123} part
      if (wick_contraction & 1){
        for (int rho_i=0; rho_i<Ns; rho_i++){
        for (int rho_f=0; rho_f<Ns; rho_f++){
          auto GAf_D1_GAi_rr_cc = GAf_D1_GAi()(rho_f,rho_i)(c_f,c_i);
          for (int alpha_f=0; alpha_f<Ns; alpha_f++){
          for (int beta_i=0; beta_i<Ns; beta_i++){
            result()(rho_f,rho_i)() += ee  * GAf_D1_GAi_rr_cc
                                        * D2_GBi    ()(alpha_f,beta_i)(a_f,a_i)
                                        * GBf_D3    ()(alpha_f,beta_i)(b_f,b_i);
          }}
        }}
      }   
      //This is the \delta_{456}^{231} part
      if (wick_contraction & 2){
        for (int rho_i=0; rho_i<Ns; rho_i++){
        for (int alpha_f=0; alpha_f<Ns; alpha_f++){
          auto D1_GAi_ar_ac = D1_GAi()(alpha_f,rho_i)(a_f,c_i);
          for (int beta_i=0; beta_i<Ns; beta_i++){
            auto GBf_D2_GBi_ab_ba = GBf_D2_GBi ()(alpha_f,beta_i)(b_f,a_i);
            for (int rho_f=0; rho_f<Ns; rho_f++){
              result()(rho_f,rho_i)() += ee  * D1_GAi_ar_ac
                                        * GBf_D2_GBi_ab_ba
                                        * GAf_D3        ()(rho_f,beta_i)(c_f,b_i);
            }
          }
        }}
      }   
      //This is the \delta_{456}^{312} part
      if (wick_contraction & 4){
        for (int rho_i=0; rho_i<Ns; rho_i++){
        for (int alpha_f=0; alpha_f<Ns; alpha_f++){
          auto GBf_D1_GAi_ar_bc = GBf_D1_GAi()(alpha_f,rho_i)(b_f,c_i);
          for (int beta_i=0; beta_i<Ns; beta_i++){
            auto D3_ab_ab = D3 ()(alpha_f,beta_i)(a_f,b_i);
            for (int rho_f=0; rho_f<Ns; rho_f++){
              result()(rho_f,rho_i)() += ee  * GBf_D1_GAi_ar_bc
                                        * GAf_D2_GBi    ()(rho_f,beta_i)(c_f,a_i)
                                        * D3_ab_ab;
            }
          }
        }}
      }   
      //This is the \delta_{456}^{132} part
      if (wick_contraction & 8){
        for (int rho_i=0; rho_i<Ns; rho_i++){
        for (int rho_f=0; rho_f<Ns; rho_f++){
          auto GAf_D1_GAi_rr_cc = GAf_D1_GAi()(rho_f,rho_i)(c_f,c_i);
          for (int alpha_f=0; alpha_f<Ns; alpha_f++){
          for (int beta_i=0; beta_i<Ns; beta_i++){
            result()(rho_f,rho_i)() -= ee  * GAf_D1_GAi_rr_cc
                                        * GBf_D2_GBi    ()(alpha_f,beta_i)(b_f,a_i)
                                        * D3            ()(alpha_f,beta_i)(a_f,b_i);
          }}
        }}
      }   
      //This is the \delta_{456}^{321} part
      if (wick_contraction & 16){
        for (int rho_i=0; rho_i<Ns; rho_i++){
        for (int alpha_f=0; alpha_f<Ns; alpha_f++){
          auto GBf_D1_GAi_ar_bc = GBf_D1_GAi()(alpha_f,rho_i)(b_f,c_i);
          for (int beta_i=0; beta_i<Ns; beta_i++){
            auto D2_GBi_ab_aa = D2_GBi()(alpha_f,beta_i)(a_f,a_i);
            for (int rho_f=0; rho_f<Ns; rho_f++){
              result()(rho_f,rho_i)() -= ee  * GBf_D1_GAi_ar_bc
                                        * D2_GBi_ab_aa
                                        * GAf_D3    ()(rho_f,beta_i)(c_f,b_i);
            }
          }
        }}
      }   
      //This is the \delta_{456}^{213} part
      if (wick_contraction & 32){
        for (int rho_i=0; rho_i<Ns; rho_i++){
        for (int alpha_f=0; alpha_f<Ns; alpha_f++){
          auto D1_GAi_ar_ac = D1_GAi()(alpha_f,rho_i)(a_f,c_i);
          for (int beta_i=0; beta_i<Ns; beta_i++){
            auto GBf_D3_ab_bb = GBf_D3()(alpha_f,beta_i)(b_f,b_i);
            for (int rho_f=0; rho_f<Ns; rho_f++){
              result()(rho_f,rho_i)() -= ee  * D1_GAi_ar_ac
                                        * GAf_D2_GBi    ()(rho_f,beta_i)(c_f,a_i)
                                        * GBf_D3_ab_bb;
            }
          }
        }}
      }
    }
  }
}

/* Computes which wick contractions should be performed for a    *
 * baryon 2pt function given the initial and finals state quark  *
 * flavours.                                                     *
 * The array wick_contractions must be of length 6               */
template<class FImpl>
void BaryonUtils<FImpl>::WickContractions(std::string qi, std::string qf, int &wick_contractions) {
    assert(qi.size() == 3 && qf.size() == 3 && "Only sets of 3 quarks accepted.");
    const int epsilon[6][3] = {{0,1,2},{1,2,0},{2,0,1},{0,2,1},{2,1,0},{1,0,2}};
    wick_contractions=0;
    for (int ie=0; ie < 6 ; ie++) {
        wick_contractions += ( ( qi[0] == qf[epsilon[ie][0]] 
                                    && qi[1] == qf[epsilon[ie][1]] 
                                    && qi[2] == qf[epsilon[ie][2]]) ? 1 : 0) << ie;
    }
}

/* The array wick_contractions must be of length 6. The order     * 
 * corresponds to the to that shown in the Hadrons documentation  *
 * at https://aportelli.github.io/Hadrons-doc/#/mcontraction      *
 * This can be computed from the quark flavours using the         *
 * Wick_Contractions function above                               */
template<class FImpl>
void BaryonUtils<FImpl>::ContractBaryons(const PropagatorField &q1_left,
             const PropagatorField &q2_left,
             const PropagatorField &q3_left,
             const Gamma GammaA_left,
             const Gamma GammaB_left,
             const Gamma GammaA_right,
             const Gamma GammaB_right,
             const int wick_contractions,
             const int parity,
             ComplexField &baryon_corr)
{

  assert(Ns==4 && "Baryon code only implemented for N_spin = 4");
  assert(Nc==3 && "Baryon code only implemented for N_colour = 3");
 
  assert(parity==1 || parity == -1 && "Parity must be +1 or -1");

  GridBase *grid = q1_left.Grid();
  
  autoView(vbaryon_corr , baryon_corr , AcceleratorWrite);
  autoView( v1          , q1_left     , AcceleratorRead);
  autoView( v2          , q2_left     , AcceleratorRead);
  autoView( v3          , q3_left     , AcceleratorRead);

  Real bytes =0.;
  bytes += grid->oSites() * (432.*sizeof(vComplex) + 126.*sizeof(int) + 36.*sizeof(Real));
  for (int ie=0; ie < 6 ; ie++){
    if(ie==0 or ie==3){
       bytes += ( wick_contractions & (1 << ie) ) ? grid->oSites() * (4.*sizeof(int) + 4752.*sizeof(vComplex)) : 0.;
    } else{
       bytes += ( wick_contractions & (1 << ie) ) ? grid->oSites() * (64.*sizeof(int) + 5184.*sizeof(vComplex)) : 0.;
    }
  }
  Real t=0.;
  t =-usecond();

  accelerator_for(ss, grid->oSites(), grid->Nsimd(), {
    auto D1 = v1(ss);
    auto D2 = v2(ss);
    auto D3 = v3(ss);
    typedef decltype(coalescedRead(vbaryon_corr[0])) cVec;
    cVec result=Zero();
    BaryonSite(D1,D2,D3,GammaA_left,GammaB_left,GammaA_right,GammaB_right,parity,wick_contractions,result);
    coalescedWrite(vbaryon_corr[ss],result);
  }  );//end loop over lattice sites

  t += usecond();

  std::cout << GridLogDebug << std::setw(10) << bytes/t*1.0e6/1024/1024/1024 << " GB/s " << std::endl;
}

template<class FImpl>
void BaryonUtils<FImpl>::ContractBaryonsMatrix(const PropagatorField &q1_left,
             const PropagatorField &q2_left,
             const PropagatorField &q3_left,
             const Gamma GammaA_left,
             const Gamma GammaB_left,
             const Gamma GammaA_right,
             const Gamma GammaB_right,
             const int wick_contractions,
             SpinMatrixField &baryon_corr)
{

  assert(Ns==4 && "Baryon code only implemented for N_spin = 4");
  assert(Nc==3 && "Baryon code only implemented for N_colour = 3");

  GridBase *grid = q1_left.Grid();

  autoView(vbaryon_corr , baryon_corr , AcceleratorWrite);
  autoView( v1          , q1_left     , AcceleratorRead);
  autoView( v2          , q2_left     , AcceleratorRead);
  autoView( v3          , q3_left     , AcceleratorRead);

  accelerator_for(ss, grid->oSites(), grid->Nsimd(), {
    auto D1 = v1(ss);
    auto D2 = v2(ss);
    auto D3 = v3(ss);
    typedef decltype(coalescedRead(vbaryon_corr[0])) spinor;
    spinor result=Zero();
    BaryonSiteMatrix(D1,D2,D3,GammaA_left,GammaB_left,GammaA_right,GammaB_right,wick_contractions,result);
    coalescedWrite(vbaryon_corr[ss],result);
  }  );//end loop over lattice sites
}

/* The array wick_contractions must be of length 6. The order     * 
 * corresponds to the to that shown in the Hadrons documentation  *
 * at https://aportelli.github.io/Hadrons-doc/#/mcontraction      *
 * This can also be computed from the quark flavours using the    *
 * Wick_Contractions function above                               */
template <class FImpl>
template <class mobj, class robj>
void BaryonUtils<FImpl>::ContractBaryonsSliced(const mobj &D1,
             const mobj &D2,
             const mobj &D3,
                         const Gamma GammaA_left,
                         const Gamma GammaB_left,
                         const Gamma GammaA_right,
                         const Gamma GammaB_right,
             const int wick_contractions,
             const int parity,
             const int nt,
             robj &result)
{

  assert(Ns==4 && "Baryon code only implemented for N_spin = 4");
  assert(Nc==3 && "Baryon code only implemented for N_colour = 3");
 
  assert(parity==1 || parity == -1 && "Parity must be +1 or -1");

  for (int t=0; t<nt; t++) {
    BaryonSite(D1[t],D2[t],D3[t],GammaA_left,GammaB_left,GammaA_right,GammaB_right,parity,wick_contractions,result[t]);
  }
}

template <class FImpl>
template <class mobj, class robj>
void BaryonUtils<FImpl>::ContractBaryonsSlicedMatrix(const mobj &D1,
             const mobj &D2,
             const mobj &D3,
             const Gamma GammaA_left,
             const Gamma GammaB_left,
             const Gamma GammaA_right,
             const Gamma GammaB_right,
             const int wick_contractions,
             const int nt,
             robj &result)
{

  assert(Ns==4 && "Baryon code only implemented for N_spin = 4");
  assert(Nc==3 && "Baryon code only implemented for N_colour = 3");

  for (int t=0; t<nt; t++) {
    BaryonSiteMatrix(D1[t],D2[t],D3[t],GammaA_left,GammaB_left,GammaA_right,GammaB_right,wick_contractions,result[t]);
  }
}

/***********************************************************************
 * End of Baryon 2pt-function code.                                    *
 *                                                                     *
 * The following code is for baryonGamma3pt function                   *
 **********************************************************************/

/* Dq1_ti is a quark line from t_i to t_J
 * Dq2_spec is a quark line from t_i to t_f
 * Dq3_spec is a quark line from t_i to t_f
 * Dq4_tf is a quark line from t_f to t_J */
template<class FImpl>
template <class mobj, class mobj2, class robj> accelerator_inline
void BaryonUtils<FImpl>::BaryonGamma3ptGroup1Site(
                        const mobj &Dq1_ti,
                        const mobj2 &Dq2_spec,
                        const mobj2 &Dq3_spec,
                        const mobj &Dq4_tf,
                        const Gamma GammaJ,
                        const Gamma GammaBi,
                        const Gamma GammaBf,
                        int wick_contraction,
                        robj &result)
{
  Gamma g5(Gamma::Algebra::Gamma5); 

  auto adjD4          = g5 * adj(Dq4_tf) * g5 ;
  auto adjD4_g_D1     = adjD4 * GammaJ * Dq1_ti;
  auto Gf_adjD4_g_D1  = GammaBf * adjD4_g_D1;
  auto D2_Gi          = Dq2_spec * GammaBi;
  auto Gf_D2_Gi       = GammaBf * D2_Gi;
  auto Gf_D3          = GammaBf * Dq3_spec;  

  Real ee;

  for (int ie_f=0; ie_f < 6 ; ie_f++){
    int a_f    = (ie_f < 3 ? ie_f       : (6-ie_f)%3 ); //epsilon[ie_n][0]; //a
    int b_f    = (ie_f < 3 ? (ie_f+1)%3 : (8-ie_f)%3 ); //epsilon[ie_n][1]; //b
    int c_f    = (ie_f < 3 ? (ie_f+2)%3 : (7-ie_f)%3 ); //epsilon[ie_n][2]; //c
    int eSgn_f = (ie_f < 3 ? 1 : -1);
    for (int ie_i=0; ie_i < 6 ; ie_i++){
      int a_i = (ie_i < 3 ? ie_i       : (6-ie_i)%3 ); //epsilon[ie_s][0]; //a'
      int b_i = (ie_i < 3 ? (ie_i+1)%3 : (8-ie_i)%3 ); //epsilon[ie_s][1]; //b'
      int c_i = (ie_i < 3 ? (ie_i+2)%3 : (7-ie_i)%3 ); //epsilon[ie_s][2]; //c'
      int eSgn_i = (ie_i < 3 ? 1 : -1);

      ee = Real(eSgn_f * eSgn_i);

      for (int alpha_f=0; alpha_f<Ns; alpha_f++){
      for (int beta_i=0; beta_i<Ns; beta_i++){
        auto D2_Gi_ab_aa        = D2_Gi     ()(alpha_f,beta_i)(a_f,a_i);
        auto Gf_D3_ab_bb        = Gf_D3     ()(alpha_f,beta_i)(b_f,b_i);
        auto Gf_D2_Gi_ab_ba     = Gf_D2_Gi  ()(alpha_f,beta_i)(b_f,a_i);
        auto Dq3_spec_ab_ab     = Dq3_spec  ()(alpha_f,beta_i)(a_f,b_i);

        for (int gamma_i=0; gamma_i<Ns; gamma_i++){
          auto ee_adjD4_g_D1_ag_ac        = ee * adjD4_g_D1   ()(alpha_f,gamma_i)(a_f,c_i);
          auto ee_Gf_adjD4_g_D1_ag_bc     = ee * Gf_adjD4_g_D1()(alpha_f,gamma_i)(b_f,c_i);
          for (int gamma_f=0; gamma_f<Ns; gamma_f++){
            auto ee_adjD4_g_D1_gg_cc        = ee * adjD4_g_D1   ()(gamma_f,gamma_i)(c_f,c_i);
            auto Dq3_spec_gb_cb             = Dq3_spec          ()(gamma_f,beta_i)(c_f,b_i);
            auto D2_Gi_gb_ca                = D2_Gi             ()(gamma_f,beta_i)(c_f,a_i);


            if(wick_contraction == 1) { // Do contraction I1
              result()(gamma_f,gamma_i)() -= ee_adjD4_g_D1_gg_cc
                                                        * D2_Gi_ab_aa
                                                        * Gf_D3_ab_bb;
            }
            if(wick_contraction == 2) { // Do contraction I2
              result()(gamma_f,gamma_i)() -= ee_adjD4_g_D1_ag_ac
                                                        * Gf_D2_Gi_ab_ba
                                                        * Dq3_spec_gb_cb;
            }
            if(wick_contraction == 3) { // Do contraction I3
              result()(gamma_f,gamma_i)() -= ee_Gf_adjD4_g_D1_ag_bc
                                                        * D2_Gi_gb_ca
                                                        * Dq3_spec_ab_ab;
            }
            if(wick_contraction == 4) { // Do contraction I4
              result()(gamma_f,gamma_i)() += ee_adjD4_g_D1_gg_cc
                                                        * Gf_D2_Gi_ab_ba
                                                        * Dq3_spec_ab_ab;
            }
            if(wick_contraction == 5) { // Do contraction I5
              result()(gamma_f,gamma_i)() += ee_Gf_adjD4_g_D1_ag_bc
                                                        * D2_Gi_ab_aa
                                                        * Dq3_spec_gb_cb;
            }
            if(wick_contraction == 6) { // Do contraction I6
              result()(gamma_f,gamma_i)() += ee_adjD4_g_D1_ag_ac
                                                        * D2_Gi_gb_ca
                                                        * Gf_D3_ab_bb;
            }
          }
        }
      }}
    }
  }
}

/* Dq1_spec is a quark line from t_i to t_f
 * Dq2_ti is a quark line from t_i to t_J
 * Dq3_spec is a quark line from t_i to t_f
 * Dq4_tf is a quark line from t_f to t_J */
template<class FImpl>
template <class mobj, class mobj2, class robj> accelerator_inline
void BaryonUtils<FImpl>::BaryonGamma3ptGroup2Site(
                        const mobj2 &Dq1_spec,
                        const mobj &Dq2_ti,
                        const mobj2 &Dq3_spec,
                        const mobj &Dq4_tf,
                        const Gamma GammaJ,
                        const Gamma GammaBi,
                        const Gamma GammaBf,
                        int wick_contraction,
                        robj &result)
{
  Gamma g5(Gamma::Algebra::Gamma5); 

  auto adjD4_g_D2_Gi      = g5 * adj(Dq4_tf) * g5 * GammaJ * Dq2_ti * GammaBi;
  auto Gf_adjD4_g_D2_Gi   = GammaBf * adjD4_g_D2_Gi;
  auto Gf_D1              = GammaBf * Dq1_spec;
  auto Gf_D3              = GammaBf * Dq3_spec;

  Real ee;

  for (int ie_f=0; ie_f < 6 ; ie_f++){
    int a_f    = (ie_f < 3 ? ie_f       : (6-ie_f)%3 ); //epsilon[ie_n][0]; //a
    int b_f    = (ie_f < 3 ? (ie_f+1)%3 : (8-ie_f)%3 ); //epsilon[ie_n][1]; //b
    int c_f    = (ie_f < 3 ? (ie_f+2)%3 : (7-ie_f)%3 ); //epsilon[ie_n][2]; //c
    int eSgn_f = (ie_f < 3 ? 1 : -1);
    for (int ie_i=0; ie_i < 6 ; ie_i++){
      int a_i = (ie_i < 3 ? ie_i       : (6-ie_i)%3 ); //epsilon[ie_s][0]; //a'
      int b_i = (ie_i < 3 ? (ie_i+1)%3 : (8-ie_i)%3 ); //epsilon[ie_s][1]; //b'
      int c_i = (ie_i < 3 ? (ie_i+2)%3 : (7-ie_i)%3 ); //epsilon[ie_s][2]; //c'
      int eSgn_i = (ie_i < 3 ? 1 : -1);

      ee = Real(eSgn_f * eSgn_i); //epsilon_sgn[ie_n] * epsilon_sgn[ie_s];

      for (int alpha_f=0; alpha_f<Ns; alpha_f++){
      for (int beta_i=0; beta_i<Ns; beta_i++){
        auto adjD4_g_D2_Gi_ab_aa        = adjD4_g_D2_Gi     ()(alpha_f,beta_i)(a_f,a_i);
        auto Gf_D3_ab_bb                = Gf_D3             ()(alpha_f,beta_i)(b_f,b_i);
        auto Gf_adjD4_g_D2_Gi_ab_ba     = Gf_adjD4_g_D2_Gi  ()(alpha_f,beta_i)(b_f,a_i);
        auto Dq3_spec_ab_ab             = Dq3_spec          ()(alpha_f,beta_i)(a_f,b_i);

        for (int gamma_i=0; gamma_i<Ns; gamma_i++){ 
          auto ee_Dq1_spec_ag_ac      = ee * Dq1_spec     ()(alpha_f,gamma_i)(a_f,c_i);
          auto ee_Gf_D1_ag_bc         = ee * Gf_D1        ()(alpha_f,gamma_i)(b_f,c_i);
          for (int gamma_f=0; gamma_f<Ns; gamma_f++){
            auto ee_Dq1_spec_gg_cc      = ee * Dq1_spec     ()(gamma_f,gamma_i)(c_f,c_i);
            auto Dq3_spec_gb_cb         = Dq3_spec          ()(gamma_f,beta_i)(c_f,b_i);
            auto adjD4_g_D2_Gi_gb_ca    = adjD4_g_D2_Gi     ()(gamma_f,beta_i)(c_f,a_i);

            if(wick_contraction == 1) { // Do contraction II1
              result()(gamma_f,gamma_i)() -= ee_Dq1_spec_gg_cc
                                                        * adjD4_g_D2_Gi_ab_aa
                                                        * Gf_D3_ab_bb;
            }
            if(wick_contraction == 2) { // Do contraction II2
              result()(gamma_f,gamma_i)() -= ee_Dq1_spec_ag_ac
                                                        * Gf_adjD4_g_D2_Gi_ab_ba
                                                        * Dq3_spec_gb_cb;
            }
            if(wick_contraction == 3) { // Do contraction II3
              result()(gamma_f,gamma_i)() -= ee_Gf_D1_ag_bc
                                                        * adjD4_g_D2_Gi_gb_ca
                                                        * Dq3_spec_ab_ab;
            }
            if(wick_contraction == 4) { // Do contraction II4
              result()(gamma_f,gamma_i)() += ee_Dq1_spec_gg_cc
                                                        * Gf_adjD4_g_D2_Gi_ab_ba
                                                        * Dq3_spec_ab_ab;
            }
            if(wick_contraction == 5) { // Do contraction II5
              result()(gamma_f,gamma_i)() += ee_Gf_D1_ag_bc
                                                        * adjD4_g_D2_Gi_ab_aa
                                                        * Dq3_spec_gb_cb;
            }
            if(wick_contraction == 6) { // Do contraction II6
              result()(gamma_f,gamma_i)() += ee_Dq1_spec_ag_ac
                                                        * adjD4_g_D2_Gi_gb_ca
                                                        * Gf_D3_ab_bb;
            }
          }
        }
      }}
    }
  }
}

/* Dq1_spec is a quark line from t_i to t_f
 * Dq2_spec is a quark line from t_i to t_f
 * Dq3_ti is a quark line from t_i to t_J
 * Dq4_tf is a quark line from t_f to t_J */
template<class FImpl>
template <class mobj, class mobj2, class robj> accelerator_inline
void BaryonUtils<FImpl>::BaryonGamma3ptGroup3Site(
                        const mobj2 &Dq1_spec,
                        const mobj2 &Dq2_spec,
                        const mobj &Dq3_ti,
                        const mobj &Dq4_tf,
                        const Gamma GammaJ,
                        const Gamma GammaBi,
                        const Gamma GammaBf,
                        int wick_contraction,
                        robj &result)
{
  Gamma g5(Gamma::Algebra::Gamma5);

  auto adjD4_g_D3     = g5 * adj(Dq4_tf) * g5 * GammaJ * Dq3_ti;
  auto Gf_adjD4_g_D3  = GammaBf * adjD4_g_D3;
  auto Gf_D1          = GammaBf * Dq1_spec;
  auto D2_Gi          = Dq2_spec * GammaBi;
  auto Gf_D2_Gi       = GammaBf * D2_Gi;

  Real ee;

  for (int ie_f=0; ie_f < 6 ; ie_f++){
    int a_f    = (ie_f < 3 ? ie_f       : (6-ie_f)%3 ); //epsilon[ie_n][0]; //a
    int b_f    = (ie_f < 3 ? (ie_f+1)%3 : (8-ie_f)%3 ); //epsilon[ie_n][1]; //b
    int c_f    = (ie_f < 3 ? (ie_f+2)%3 : (7-ie_f)%3 ); //epsilon[ie_n][2]; //c
    int eSgn_f = (ie_f < 3 ? 1 : -1);
    for (int ie_i=0; ie_i < 6 ; ie_i++){
      int a_i = (ie_i < 3 ? ie_i       : (6-ie_i)%3 ); //epsilon[ie_s][0]; //a'
      int b_i = (ie_i < 3 ? (ie_i+1)%3 : (8-ie_i)%3 ); //epsilon[ie_s][1]; //b'
      int c_i = (ie_i < 3 ? (ie_i+2)%3 : (7-ie_i)%3 ); //epsilon[ie_s][2]; //c'
      int eSgn_i = (ie_i < 3 ? 1 : -1);

      ee = Real(eSgn_f * eSgn_i); //epsilon_sgn[ie_n] * epsilon_sgn[ie_s];

      for (int alpha_f=0; alpha_f<Ns; alpha_f++){
      for (int beta_i=0; beta_i<Ns; beta_i++){
        auto D2_Gi_ab_aa            = D2_Gi         ()(alpha_f,beta_i)(a_f,a_i);
        auto Gf_adjD4_g_D3_ab_bb    = Gf_adjD4_g_D3 ()(alpha_f,beta_i)(b_f,b_i);
        auto Gf_D2_Gi_ab_ba         = Gf_D2_Gi      ()(alpha_f,beta_i)(b_f,a_i);
        auto adjD4_g_D3_ab_ab       = adjD4_g_D3    ()(alpha_f,beta_i)(a_f,b_i);

        for (int gamma_i=0; gamma_i<Ns; gamma_i++) {
          auto ee_Dq1_spec_ag_ac  = ee * Dq1_spec ()(alpha_f,gamma_i)(a_f,c_i);
          auto ee_Gf_D1_ag_bc     = ee * Gf_D1    ()(alpha_f,gamma_i)(b_f,c_i);
          for (int gamma_f=0; gamma_f<Ns; gamma_f++) {
            auto ee_Dq1_spec_gg_cc  = ee * Dq1_spec ()(gamma_f,gamma_i)(c_f,c_i);
            auto adjD4_g_D3_gb_cb   = adjD4_g_D3    ()(gamma_f,beta_i)(c_f,b_i);
            auto D2_Gi_gb_ca        = D2_Gi         ()(gamma_f,beta_i)(c_f,a_i);

            if(wick_contraction == 1) { // Do contraction III1
              result()(gamma_f,gamma_i)() -= ee_Dq1_spec_gg_cc
                                                        * D2_Gi_ab_aa
                                                        * Gf_adjD4_g_D3_ab_bb;
            }
            if(wick_contraction == 2) { // Do contraction III2
              result()(gamma_f,gamma_i)() -= ee_Dq1_spec_ag_ac
                                                        * Gf_D2_Gi_ab_ba
                                                        * adjD4_g_D3_gb_cb;
            }
            if(wick_contraction == 3) { // Do contraction III3
              result()(gamma_f,gamma_i)() -= ee_Gf_D1_ag_bc
                                                        * D2_Gi_gb_ca
                                                        * adjD4_g_D3_ab_ab;
            }
            if(wick_contraction == 4) { // Do contraction III4
              result()(gamma_f,gamma_i)() += ee_Dq1_spec_gg_cc
                                                        * Gf_D2_Gi_ab_ba
                                                        * adjD4_g_D3_ab_ab;
            }
            if(wick_contraction == 5) { // Do contraction III5
              result()(gamma_f,gamma_i)() += ee_Gf_D1_ag_bc
                                                        * D2_Gi_ab_aa
                                                        * adjD4_g_D3_gb_cb;
            }
            if(wick_contraction == 6) { // Do contraction III6
              result()(gamma_f,gamma_i)() += ee_Dq1_spec_ag_ac
                                                        * D2_Gi_gb_ca
                                                        * Gf_adjD4_g_D3_ab_bb;
            }
          }
        }
      }}
    }
  }
}

/* The group indicates which inital state quarks the current is  * 
 * connected to. It must be in the range 1-3.                    *
 * The wick_contraction must be in the range 1-6 correspond to   *
 * the contractions given in the Hadrons documentation at        *
 * https://aportelli.github.io/Hadrons-doc/#/mcontraction        */
template<class FImpl>
template <class mobj>
void BaryonUtils<FImpl>::BaryonGamma3pt(
                        const PropagatorField &q_ti,
                        const mobj &Dq_spec1,
                        const mobj &Dq_spec2,
                        const PropagatorField &q_tf,
                        int group,
                        int wick_contraction,
                        const Gamma GammaJ,
                        const Gamma GammaBi,
                        const Gamma GammaBf,
                        SpinMatrixField &stn_corr)
{
  assert(Ns==4 && "Baryon code only implemented for N_spin = 4");
  assert(Nc==3 && "Baryon code only implemented for N_colour = 3");

  GridBase *grid = q_tf.Grid();

  autoView( vcorr , stn_corr , AcceleratorWrite);
  autoView( vq_ti , q_ti     , AcceleratorRead);
  autoView( vq_tf , q_tf     , AcceleratorRead);

  Vector<mobj> my_Dq_spec{Dq_spec1,Dq_spec2};
  mobj * Dq_spec_p = &my_Dq_spec[0];

  if (group == 1) {
    accelerator_for(ss, grid->oSites(), grid->Nsimd(), {
      auto Dq_ti = vq_ti(ss);
      auto Dq_tf = vq_tf(ss);
      typedef decltype(coalescedRead(vcorr[0])) spinor;
      spinor result=Zero();
      BaryonGamma3ptGroup1Site(Dq_ti,Dq_spec_p[0],Dq_spec_p[1],Dq_tf,GammaJ,GammaBi,GammaBf,wick_contraction,result);
      coalescedWrite(vcorr[ss],coalescedRead(vcorr[ss])+result); 
    });//end loop over lattice sites
  } else if (group == 2) {
    accelerator_for(ss, grid->oSites(), grid->Nsimd(), {
      auto Dq_ti = vq_ti(ss);
      auto Dq_tf = vq_tf(ss);
      typedef decltype(coalescedRead(vcorr[0])) spinor;
      spinor result=Zero();
      BaryonGamma3ptGroup2Site(Dq_spec_p[0],Dq_ti,Dq_spec_p[1],Dq_tf,GammaJ,GammaBi,GammaBf,wick_contraction,result); 
      coalescedWrite(vcorr[ss],coalescedRead(vcorr[ss])+result); 
    });//end loop over lattice sites
  } else if (group == 3) {
    accelerator_for(ss, grid->oSites(), grid->Nsimd(), {
      auto Dq_ti = vq_ti(ss);
      auto Dq_tf = vq_tf(ss);
      typedef decltype(coalescedRead(vcorr[0])) spinor;
      spinor result=Zero();
      BaryonGamma3ptGroup3Site(Dq_spec_p[0],Dq_spec_p[1],Dq_ti,Dq_tf,GammaJ,GammaBi,GammaBf,wick_contraction,result); 
      coalescedWrite(vcorr[ss],coalescedRead(vcorr[ss])+result); 
    });//end loop over lattice sites
  }

}


/***********************************************************************
 * End of BaryonGamma3pt-function code.                                *
 *                                     *
 * The following code is for Sigma -> N rare hypeon decays             *
 **********************************************************************/

/* Dq_loop is a quark line from t_H to t_H
 * Du_spec is a quark line from t_i to t_f
 * Dd_tf is a quark line from t_f to t_H
 * Ds_ti is a quark line from t_i to t_H */
template <class FImpl>
template <class mobj, class mobj2, class robj> accelerator_inline
void BaryonUtils<FImpl>::SigmaToNucleonQ1EyeSite(const mobj &Dq_loop,
             const mobj2 &Du_spec,
             const mobj &Dd_tf,
             const mobj &Ds_ti,
             const Gamma Gamma_H,
             const Gamma GammaB_sigma,
             const Gamma GammaB_nucl,
             robj &result)
{

  Gamma g5(Gamma::Algebra::Gamma5); 

  auto adjDd_GH_Ds      = g5 * adj(Dd_tf) * g5 * Gamma_H * Ds_ti;
  auto Gn_adjDd_GH_Ds   = GammaB_nucl * adjDd_GH_Ds;
  auto Du_Gs            = Du_spec * GammaB_sigma;
  auto Dq_GH            = Dq_loop * Gamma_H;
  auto Tr_Dq_GH         = trace(Dq_GH)()()();

  Real ee;

  for (int ie_n=0; ie_n < 6 ; ie_n++){
    int a_n    = (ie_n < 3 ? ie_n       : (6-ie_n)%3 ); //epsilon[ie_n][0]; //a
    int b_n    = (ie_n < 3 ? (ie_n+1)%3 : (8-ie_n)%3 ); //epsilon[ie_n][1]; //b
    int c_n    = (ie_n < 3 ? (ie_n+2)%3 : (7-ie_n)%3 ); //epsilon[ie_n][2]; //c
    int eSgn_n = (ie_n < 3 ? 1 : -1);
    for (int ie_s=0; ie_s < 6 ; ie_s++){
      int a_s = (ie_s < 3 ? ie_s       : (6-ie_s)%3 ); //epsilon[ie_s][0]; //a'
      int b_s = (ie_s < 3 ? (ie_s+1)%3 : (8-ie_s)%3 ); //epsilon[ie_s][1]; //b'
      int c_s = (ie_s < 3 ? (ie_s+2)%3 : (7-ie_s)%3 ); //epsilon[ie_s][2]; //c'
      int eSgn_s = (ie_s < 3 ? 1 : -1);

      ee = Real(eSgn_n * eSgn_s); 
      for (int alpha_n=0; alpha_n<Ns; alpha_n++){
      for (int beta_s=0;  beta_s<Ns;  beta_s++){

        auto Gn_adjDd_GH_Ds_ab_bb = Gn_adjDd_GH_Ds ()(alpha_n, beta_s)(b_n,b_s);

        for (int gamma_s=0; gamma_s<Ns; gamma_s++){
        for (int gamma_n=0; gamma_n<Ns; gamma_n++){
          result()(gamma_n,gamma_s)() += ee   * Gn_adjDd_GH_Ds_ab_bb
                                            * Du_spec         ()(gamma_n,gamma_s)(c_n,c_s)
                                            * Du_Gs           ()(alpha_n, beta_s)(a_n,a_s) 
                                            * Tr_Dq_GH;

          result()(gamma_n,gamma_s)() -= ee   * Gn_adjDd_GH_Ds_ab_bb
                                            * Du_spec         ()(alpha_n,gamma_s)(a_n,c_s)
                                            * Du_Gs           ()(gamma_n, beta_s)(c_n,a_s) 
                                            * Tr_Dq_GH;
        }}
      }}
    }
  }
}

/* Du_ti is a quark line from t_i to t_H
 * Du_tf is a quark line from t_f to t_H
 * Du_spec is a quark line from t_i to t_f
 * Dd_tf is a quark line from t_f to t_H
 * Ds_ti is a quark line from t_i to t_H */
template <class FImpl>
template <class mobj, class mobj2, class robj> accelerator_inline
void BaryonUtils<FImpl>::SigmaToNucleonQ1NonEyeSite(const mobj &Du_ti,
             const mobj &Du_tf,
             const mobj2 &Du_spec,
             const mobj &Dd_tf,
             const mobj &Ds_ti,
             const Gamma Gamma_H,
             const Gamma GammaB_sigma,
             const Gamma GammaB_nucl,
             robj &result)
{

  Gamma g5(Gamma::Algebra::Gamma5); 

  auto Du_Gs          = Du_spec * GammaB_sigma;
  auto adjDd_GH_Ds = g5 * adj(Dd_tf) * g5 * Gamma_H * Ds_ti;
  auto Gn_adjDd_GH_Ds = GammaB_nucl * adjDd_GH_Ds;
  auto adjDu_GH_Du    = g5 * adj(Du_tf) * g5 * Gamma_H * Du_ti;
  auto adjDu_GH_Du_Gs = adjDu_GH_Du * GammaB_sigma;

  Real ee;

  for (int ie_n=0; ie_n < 6 ; ie_n++){
    int a_n    = (ie_n < 3 ? ie_n       : (6-ie_n)%3 ); //epsilon[ie_n][0]; //a
    int b_n    = (ie_n < 3 ? (ie_n+1)%3 : (8-ie_n)%3 ); //epsilon[ie_n][1]; //b
    int c_n    = (ie_n < 3 ? (ie_n+2)%3 : (7-ie_n)%3 ); //epsilon[ie_n][2]; //c
    int eSgn_n = (ie_n < 3 ? 1 : -1);
    for (int ie_s=0; ie_s < 6 ; ie_s++){
      int a_s = (ie_s < 3 ? ie_s       : (6-ie_s)%3 ); //epsilon[ie_s][0]; //a'
      int b_s = (ie_s < 3 ? (ie_s+1)%3 : (8-ie_s)%3 ); //epsilon[ie_s][1]; //b'
      int c_s = (ie_s < 3 ? (ie_s+2)%3 : (7-ie_s)%3 ); //epsilon[ie_s][2]; //c'
      int eSgn_s = (ie_s < 3 ? 1 : -1);

      ee = Real(eSgn_n * eSgn_s); //epsilon_sgn[ie_n] * epsilon_sgn[ie_s];

      for (int alpha_n=0; alpha_n<Ns; alpha_n++){
      for (int beta_s=0;  beta_s<Ns;  beta_s++){

        auto Gn_adjDd_GH_Ds_ab_bb = Gn_adjDd_GH_Ds  ()(alpha_n, beta_s)(b_n,b_s);

        for (int gamma_s=0; gamma_s<Ns; gamma_s++){
        for (int gamma_n=0; gamma_n<Ns; gamma_n++){

          result()(gamma_n,gamma_s)() += ee * Gn_adjDd_GH_Ds_ab_bb
                                          * adjDu_GH_Du     ()(alpha_n,gamma_s)(a_n,c_s)
                                          * Du_Gs           ()(gamma_n, beta_s)(c_n,a_s);

          result()(gamma_n,gamma_s)() += ee * Gn_adjDd_GH_Ds_ab_bb
                                          * adjDu_GH_Du_Gs  ()(gamma_n, beta_s)(c_n,a_s)
                                          * Du_spec         ()(alpha_n,gamma_s)(a_n,c_s);

          result()(gamma_n,gamma_s)() -= ee * Gn_adjDd_GH_Ds_ab_bb
                                          * adjDu_GH_Du_Gs  ()(alpha_n, beta_s)(a_n,a_s)
                                          * Du_spec         ()(gamma_n,gamma_s)(c_n,c_s);

          result()(gamma_n,gamma_s)() -= ee * Gn_adjDd_GH_Ds_ab_bb
                                          * adjDu_GH_Du     ()(gamma_n,gamma_s)(c_n,c_s)
                                          * Du_Gs           ()(alpha_n, beta_s)(a_n,a_s);
        }}
      }}
    }
  }
}

//Equivalent to "One-trace"
/* Dq_loop is a quark line from t_H to t_H
 * Du_spec is a quark line from t_i to t_f
 * Dd_tf is a quark line from t_f to t_H
 * Ds_ti is a quark line from t_i to t_H */
template <class FImpl>
template <class mobj, class mobj2, class robj> accelerator_inline
void BaryonUtils<FImpl>::SigmaToNucleonQ2EyeSite(const mobj &Dq_loop,
             const mobj2 &Du_spec,
             const mobj &Dd_tf,
             const mobj &Ds_ti,
             const Gamma Gamma_H,
             const Gamma GammaB_sigma,
             const Gamma GammaB_nucl,
             robj &result)
{

  Gamma g5(Gamma::Algebra::Gamma5); 

  auto adjDd_GH_Duloop_GH_Ds = g5 * adj(Dd_tf) * g5 * Gamma_H * Dq_loop * Gamma_H * Ds_ti;
  auto Gn_adjDd_GH_Duloop_GH_Ds = GammaB_nucl * adjDd_GH_Duloop_GH_Ds;
  auto Du_Gs = Du_spec * GammaB_sigma;

  Real ee;

  for (int ie_n=0; ie_n < 6 ; ie_n++){
    int a_n    = (ie_n < 3 ? ie_n       : (6-ie_n)%3 ); //epsilon[ie_n][0]; //a
    int b_n    = (ie_n < 3 ? (ie_n+1)%3 : (8-ie_n)%3 ); //epsilon[ie_n][1]; //b
    int c_n    = (ie_n < 3 ? (ie_n+2)%3 : (7-ie_n)%3 ); //epsilon[ie_n][2]; //c
    int eSgn_n = (ie_n < 3 ? 1 : -1);
    for (int ie_s=0; ie_s < 6 ; ie_s++){
      int a_s = (ie_s < 3 ? ie_s       : (6-ie_s)%3 ); //epsilon[ie_s][0]; //a'
      int b_s = (ie_s < 3 ? (ie_s+1)%3 : (8-ie_s)%3 ); //epsilon[ie_s][1]; //b'
      int c_s = (ie_s < 3 ? (ie_s+2)%3 : (7-ie_s)%3 ); //epsilon[ie_s][2]; //c'
      int eSgn_s = (ie_s < 3 ? 1 : -1);

      ee = Real(eSgn_n * eSgn_s); //epsilon_sgn[ie_n] * epsilon_sgn[ie_s];

      for (int alpha_n=0; alpha_n<Ns; alpha_n++){
      for (int beta_s=0; beta_s<Ns; beta_s++){

        auto Gn_adjDd_GH_Duloop_GH_Ds_ab_bb = Gn_adjDd_GH_Duloop_GH_Ds ()(alpha_n,beta_s)(b_n,b_s);

        for (int gamma_s=0; gamma_s<Ns; gamma_s++){
        for (int gamma_n=0; gamma_n<Ns; gamma_n++){

          result()(gamma_n,gamma_s)() -= ee   * Du_spec         ()(gamma_n,gamma_s)(c_n,c_s)
                                            * Du_Gs           ()(alpha_n,beta_s)(a_n,a_s) 
                                            * Gn_adjDd_GH_Duloop_GH_Ds_ab_bb;

          result()(gamma_n,gamma_s)() += ee   * Du_Gs         ()(alpha_n,gamma_s)(a_n,c_s)
                                            * Du_spec       ()(gamma_n,beta_s)(c_n,a_s) 
                                            * Gn_adjDd_GH_Duloop_GH_Ds_ab_bb;
        }}
      }}
    }
  }
}

/* Du_ti is a quark line from t_i to t_H
 * Du_tf is a quark line from t_f to t_H
 * Du_spec is a quark line from t_i to t_f
 * Dd_tf is a quark line from t_f to t_H
 * Ds_ti is a quark line from t_i to t_H */
template <class FImpl>
template <class mobj, class mobj2, class robj> accelerator_inline
void BaryonUtils<FImpl>::SigmaToNucleonQ2NonEyeSite(const mobj &Du_ti,
             const mobj &Du_tf,
             const mobj2 &Du_spec,
             const mobj &Dd_tf,
             const mobj &Ds_ti,
             const Gamma Gamma_H,
             const Gamma GammaB_sigma,
             const Gamma GammaB_nucl,
             robj &result)
{

  Gamma g5(Gamma::Algebra::Gamma5); 

  auto Du_Gs              = Du_spec * GammaB_sigma;
  auto adjDu_GH_Ds        = g5 * adj(Du_tf) * g5 * Gamma_H * Ds_ti;
  auto adjDd_GH_Du        = g5 * adj(Dd_tf) * g5 * Gamma_H * Du_ti;
  auto Gn_adjDd_GH_Du     = GammaB_nucl * adjDd_GH_Du; // for some reason I needed to split this into two lines to avoid the compilation error 'error: identifier "Grid::Gamma::mul" is undefined in device code'

  auto Gn_adjDd_GH_Du_Gs  = Gn_adjDd_GH_Du * GammaB_sigma;

  Real ee;

  for (int ie_n=0; ie_n < 6 ; ie_n++){
    int a_n    = (ie_n < 3 ? ie_n       : (6-ie_n)%3 ); //epsilon[ie_n][0]; //a
    int b_n    = (ie_n < 3 ? (ie_n+1)%3 : (8-ie_n)%3 ); //epsilon[ie_n][1]; //b
    int c_n    = (ie_n < 3 ? (ie_n+2)%3 : (7-ie_n)%3 ); //epsilon[ie_n][2]; //c
    int eSgn_n = (ie_n < 3 ? 1 : -1);
    for (int ie_s=0; ie_s < 6 ; ie_s++){
      int a_s = (ie_s < 3 ? ie_s       : (6-ie_s)%3 ); //epsilon[ie_s][0]; //a'
      int b_s = (ie_s < 3 ? (ie_s+1)%3 : (8-ie_s)%3 ); //epsilon[ie_s][1]; //b'
      int c_s = (ie_s < 3 ? (ie_s+2)%3 : (7-ie_s)%3 ); //epsilon[ie_s][2]; //c'
      int eSgn_s = (ie_s < 3 ? 1 : -1);

      ee = Real(eSgn_n * eSgn_s); //epsilon_sgn[ie_n] * epsilon_sgn[ie_s];

      for (int alpha_n=0; alpha_n<Ns; alpha_n++){
      for (int beta_s=0;  beta_s<Ns;   beta_s++){

        auto adjDu_GH_Ds_ab_ab = adjDu_GH_Ds()(alpha_n, beta_s)(a_n,b_s);
        auto Gn_adjDd_GH_Du_Gs_ab_ba = Gn_adjDd_GH_Du_Gs()(alpha_n, beta_s)(b_n,a_s);

        for (int gamma_s=0; gamma_s<Ns; gamma_s++){
          auto Gn_adjDd_GH_Du_ag_bc = Gn_adjDd_GH_Du()(alpha_n,gamma_s)(b_n,c_s);
          for (int gamma_n=0; gamma_n<Ns; gamma_n++){
            auto adjDu_GH_Ds_gb_cb = adjDu_GH_Ds()(gamma_n, beta_s)(c_n,b_s);

            result()(gamma_n,gamma_s)() += ee * adjDu_GH_Ds_ab_ab
                                          * Gn_adjDd_GH_Du_Gs_ab_ba
                                          * Du_spec()(gamma_n,gamma_s)(c_n,c_s);

            result()(gamma_n,gamma_s)() -= ee * adjDu_GH_Ds_gb_cb
                                          * Gn_adjDd_GH_Du_Gs_ab_ba
                                          * Du_spec()(alpha_n,gamma_s)(a_n,c_s);

            result()(gamma_n,gamma_s)() += ee * adjDu_GH_Ds_gb_cb
                                          * Gn_adjDd_GH_Du_ag_bc
                                          * Du_Gs()(alpha_n, beta_s)(a_n,a_s);

            result()(gamma_n,gamma_s)() -= ee * adjDu_GH_Ds_ab_ab
                                          * Gn_adjDd_GH_Du_ag_bc
                                          * Du_Gs()(gamma_n, beta_s)(c_n,a_s);
          }
        }
      }}
    }
  }

}

template<class FImpl>
template <class mobj>
void BaryonUtils<FImpl>::SigmaToNucleonEye(const PropagatorField &qq_loop,
             const mobj &Du_spec,
             const PropagatorField &qd_tf,
             const PropagatorField &qs_ti,
             const Gamma Gamma_H,
             const Gamma GammaB_sigma,
             const Gamma GammaB_nucl,
             const std::string op,
             SpinMatrixField &stn_corr)
{

  assert(Ns==4 && "Baryon code only implemented for N_spin = 4");
  assert(Nc==3 && "Baryon code only implemented for N_colour = 3");

  GridBase *grid = qs_ti.Grid();

  autoView( vcorr   , stn_corr , AcceleratorWrite);
  autoView( vq_loop , qq_loop  , AcceleratorRead);
  autoView( vd_tf   , qd_tf    , AcceleratorRead);
  autoView( vs_ti   , qs_ti    , AcceleratorRead);

  Vector<mobj> my_Dq_spec{Du_spec};
  mobj * Dq_spec_p = &my_Dq_spec[0];

  if(op == "Q1"){
    accelerator_for(ss, grid->oSites(), grid->Nsimd(), {
      auto Dq_loop = vq_loop(ss);
      auto Dd_tf = vd_tf(ss);
      auto Ds_ti = vs_ti(ss);
      typedef decltype(coalescedRead(vcorr[0])) spinor;
      spinor result=Zero();
      SigmaToNucleonQ1EyeSite(Dq_loop,Dq_spec_p[0],Dd_tf,Ds_ti,Gamma_H,GammaB_sigma,GammaB_nucl,result);
      coalescedWrite(vcorr[ss],result);
    });//end loop over lattice sites
  } else if(op == "Q2"){
    accelerator_for(ss, grid->oSites(), grid->Nsimd(), {
      auto Dq_loop = vq_loop(ss);
      auto Dd_tf = vd_tf(ss);
      auto Ds_ti = vs_ti(ss);
      typedef decltype(coalescedRead(vcorr[0])) spinor;
      spinor result=Zero();
      SigmaToNucleonQ2EyeSite(Dq_loop,Dq_spec_p[0],Dd_tf,Ds_ti,Gamma_H,GammaB_sigma,GammaB_nucl,result);
      coalescedWrite(vcorr[ss],result);
    });//end loop over lattice sites
  } else {
    assert(0 && "Weak Operator not correctly specified");
  }
}

template<class FImpl>
template <class mobj>
void BaryonUtils<FImpl>::SigmaToNucleonNonEye(const PropagatorField &qq_ti,
             const PropagatorField &qq_tf,
             const mobj &Du_spec,
             const PropagatorField &qd_tf,
             const PropagatorField &qs_ti,
             const Gamma Gamma_H,
             const Gamma GammaB_sigma,
             const Gamma GammaB_nucl,
             const std::string op,
             SpinMatrixField &stn_corr)
{

  assert(Ns==4 && "Baryon code only implemented for N_spin = 4");
  assert(Nc==3 && "Baryon code only implemented for N_colour = 3");

  GridBase *grid = qs_ti.Grid();

  autoView( vcorr , stn_corr , AcceleratorWrite );
  autoView( vq_ti , qq_ti    , AcceleratorRead  );
  autoView( vq_tf , qq_tf    , AcceleratorRead  );
  autoView( vd_tf , qd_tf    , AcceleratorRead  );
  autoView( vs_ti , qs_ti    , AcceleratorRead  );
  
  Vector<mobj> my_Dq_spec{Du_spec};
  mobj * Dq_spec_p = &my_Dq_spec[0];

  if(op == "Q1"){
    accelerator_for(ss, grid->oSites(), grid->Nsimd(), {
      auto Dq_ti = vq_ti(ss);
      auto Dq_tf = vq_tf(ss);
      auto Dd_tf = vd_tf(ss);
      auto Ds_ti = vs_ti(ss);
      typedef decltype(coalescedRead(vcorr[0])) spinor;
      spinor result=Zero();
      SigmaToNucleonQ1NonEyeSite(Dq_ti,Dq_tf,Dq_spec_p[0],Dd_tf,Ds_ti,Gamma_H,GammaB_sigma,GammaB_nucl,result);
      coalescedWrite(vcorr[ss],result);
    });//end loop over lattice sites
  } else if(op == "Q2"){
    accelerator_for(ss, grid->oSites(), grid->Nsimd(), {
      auto Dq_ti = vq_ti(ss);
      auto Dq_tf = vq_tf(ss);
      auto Dd_tf = vd_tf(ss);
      auto Ds_ti = vs_ti(ss);
      typedef decltype(coalescedRead(vcorr[0])) spinor;
      spinor result=Zero();
      SigmaToNucleonQ2NonEyeSite(Dq_ti,Dq_tf,Dq_spec_p[0],Dd_tf,Ds_ti,Gamma_H,GammaB_sigma,GammaB_nucl,result);
      coalescedWrite(vcorr[ss],result);
    });//end loop over lattice sites
  } else {
    assert(0 && "Weak Operator not correctly specified");
  }
}


/***********************************************************************
 * The following code is for Xi -> Sigma rare hypeon decays            *
 **********************************************************************/

/* Dq_loop is a quark line from t_H to t_H
 * Dd_spec is a quark line from t_i to t_f
 * Ds_spec is a quark line from t_i to t_f
 * Dd_tf is a quark line from t_f to t_H
 * Ds_ti is a quark line from t_i to t_H */
template <class FImpl>
template <class mobj, class mobj2, class robj> accelerator_inline
void BaryonUtils<FImpl>::XiToSigmaQ1EyeSite(const mobj &Dq_loop,
						 const mobj2 &Dd_spec,
						 const mobj2 &Ds_spec,
						 const mobj &Dd_tf,
						 const mobj &Ds_ti,
				                 const Gamma Gamma_H,
				                 const Gamma GammaB_xi,
		                 		 const Gamma GammaB_sigma,
						 robj &result)
{

  Gamma g5(Gamma::Algebra::Gamma5); 

  auto DdG = Dd_spec * GammaB_sigma;
  auto GDs = GammaB_xi * Ds_spec;
  // Ds * \gamma_\mu^L * (\gamma_5 * Dd^\dagger * \gamma_5)
  auto DsGDd = Ds_ti * Gamma_H * g5 * adj(Dd_tf) * g5;
  // DsGDd * GammaB
  auto DsGDdG = DsGDd * GammaB_sigma;
  // GammaB * DsGDd
  auto GDsGDd = GammaB_xi * DsGDd;
  // GammaB * DsGDd * GammaB
  auto GDsGDdG = GDsGDd * GammaB_sigma;
  // \gamma_\mu^L * Dq_loop 
  auto trGDq = TensorRemove(trace(Gamma_H * Dq_loop)); 

  Real ee;

  for (int ie_s=0; ie_s < 6 ; ie_s++){
      int a_s = (ie_s < 3 ? ie_s       : (6-ie_s)%3 ); //epsilon[ie_s][0]; //a'
      int b_s = (ie_s < 3 ? (ie_s+1)%3 : (8-ie_s)%3 ); //epsilon[ie_s][1]; //b'
      int c_s = (ie_s < 3 ? (ie_s+2)%3 : (7-ie_s)%3 ); //epsilon[ie_s][2]; //c'
      int eSgn_s = (ie_s < 3 ? 1 : -1);
    for (int ie_x=0; ie_x < 6 ; ie_x++){
      int a_x = (ie_x < 3 ? ie_x       : (6-ie_x)%3 ); //epsilon[ie_x][0]; //a'
      int b_x = (ie_x < 3 ? (ie_x+1)%3 : (8-ie_x)%3 ); //epsilon[ie_x][1]; //b'
      int c_x = (ie_x < 3 ? (ie_x+2)%3 : (7-ie_x)%3 ); //epsilon[ie_x][2]; //c'
      int eSgn_x = (ie_x < 3 ? 1 : -1);
      ee = Real(eSgn_s * eSgn_x);
      auto ee_GD = ee * trGDq;
      for (int alpha_x=0; alpha_x<Ns; alpha_x++){
      for (int beta_s=0; beta_s<Ns; beta_s++){
        auto GDsGDdG_ab_ba = GDsGDd()(alpha_x,beta_s)(b_x,a_s);
        auto Ds_ab_ab = Ds_spec()(alpha_x,beta_s)(a_x,b_s);
        auto DsGDdG_ab_aa = DsGDd()(alpha_x,beta_s)(a_x,a_s);
        auto GDs_ab_bb = GDs()(alpha_x,beta_s)(b_x,b_s);
        for (int gamma_x=0; gamma_x<Ns; gamma_x++){
          auto DdG_cb_ca = DdG()(gamma_x,beta_s)(c_x,a_s);
          for (int gamma_s=0; gamma_s<Ns; gamma_s++){
            result()(gamma_x,gamma_s)() -= ee_GD * Dd_spec()(gamma_x,gamma_s)(c_x,c_s) * GDsGDdG_ab_ba * Ds_ab_ab;
            result()(gamma_x,gamma_s)() += ee_GD * DdG_cb_ca * GDsGDd()(alpha_x,gamma_s)(b_x,c_s) * Ds_ab_ab;
            result()(gamma_x,gamma_s)() += ee_GD * Dd_spec()(gamma_x,gamma_s)(c_x,c_s) * DsGDdG_ab_aa * GDs_ab_bb;
            result()(gamma_x,gamma_s)() -= ee_GD * DdG_cb_ca * DsGDd()(alpha_x,gamma_s)(a_x,c_s) * GDs_ab_bb;
          }
	}
      }}
    }
  }
}


//Equivalent to "One-trace"
/* Dq_loop is a quark line from t_H to t_H
 * Dd_spec is a quark line from t_i to t_f
 * Ds_spec is a quark line from t_i to t_f
 * Dd_tf is a quark line from t_f to t_H
 * Ds_ti is a quark line from t_i to t_H */
template <class FImpl>
template <class mobj, class mobj2, class robj> accelerator_inline
void BaryonUtils<FImpl>::XiToSigmaQ2EyeSite(const mobj &Dq_loop,
						 const mobj2 &Dd_spec,
						 const mobj2 &Ds_spec,
						 const mobj &Dd_tf,
						 const mobj &Ds_ti,
				                 const Gamma Gamma_H,
				                 const Gamma GammaB_xi,
		                 		 const Gamma GammaB_sigma,
						 robj &result)
{

  Gamma g5(Gamma::Algebra::Gamma5); 

  auto DdG = Dd_spec * GammaB_sigma;
  auto GDs = GammaB_xi * Ds_spec;
  // Ds * \gamma_\mu^L * Dq_loop * \gamma_\mu^L * (\gamma_5 * Dd^\dagger * \gamma_5)
  auto DsGDqGDd = Ds_ti * Gamma_H * Dq_loop * Gamma_H * g5 * adj(Dd_tf) * g5;
  // DsGDd * GammaB
  auto DsGDqGDdG = DsGDqGDd * GammaB_sigma;
  // GammaB * DsGDd
  auto GDsGDqGDd = GammaB_xi * DsGDqGDd;
  // GammaB * DsGDd * GammaB
  auto GDsGDqGDdG = GDsGDqGDd * GammaB_sigma;

  Real ee;

  for (int ie_s=0; ie_s < 6 ; ie_s++){
      int a_s = (ie_s < 3 ? ie_s       : (6-ie_s)%3 ); //epsilon[ie_s][0]; //a'
      int b_s = (ie_s < 3 ? (ie_s+1)%3 : (8-ie_s)%3 ); //epsilon[ie_s][1]; //b'
      int c_s = (ie_s < 3 ? (ie_s+2)%3 : (7-ie_s)%3 ); //epsilon[ie_s][2]; //c'
      int eSgn_s = (ie_s < 3 ? 1 : -1);
    for (int ie_x=0; ie_x < 6 ; ie_x++){
      int a_x = (ie_x < 3 ? ie_x       : (6-ie_x)%3 ); //epsilon[ie_x][0]; //a'
      int b_x = (ie_x < 3 ? (ie_x+1)%3 : (8-ie_x)%3 ); //epsilon[ie_x][1]; //b'
      int c_x = (ie_x < 3 ? (ie_x+2)%3 : (7-ie_x)%3 ); //epsilon[ie_x][2]; //c'
      int eSgn_x = (ie_x < 3 ? 1 : -1);
      ee = Real(eSgn_s * eSgn_x);
      for (int alpha_x=0; alpha_x<Ns; alpha_x++){
      for (int beta_s=0; beta_s<Ns; beta_s++){
        auto GDsGDqGDdG_ab_ba = GDsGDqGDdG()(alpha_x,beta_s)(b_x,a_s);
        auto Ds_ab_ab = Ds_spec()(alpha_x,beta_s)(a_x,b_s);
        auto DsGDqGDdG_ab_aa = GDsGDqGDdG()(alpha_x,beta_s)(a_x,a_s);
        auto GDs_ab_bb = GDs()(alpha_x,beta_s)(b_x,b_s);
          for (int gamma_x=0; gamma_x<Ns; gamma_x++){
            auto DdG_cb_ca = DdG()(gamma_x, beta_s)(c_x,a_s);
          for (int gamma_s=0; gamma_s<Ns; gamma_s++){
            result()(gamma_x,gamma_s)() -= ee * Dd_spec()(gamma_x,gamma_s)(c_x,c_s) * GDsGDqGDdG_ab_ba * Ds_ab_ab;
            result()(gamma_x,gamma_s)() += ee * DdG_cb_ca * GDsGDqGDd()(alpha_x,gamma_s)(b_x,c_s) * Ds_ab_ab;
            result()(gamma_x,gamma_s)() += ee * Dd_spec()(gamma_x,gamma_s)(c_x,c_s) * DsGDqGDdG_ab_aa * GDs_ab_bb;
            result()(gamma_x,gamma_s)() -= ee * DdG_cb_ca * DsGDqGDd()(alpha_x,gamma_s)(a_x,c_s) * GDs_ab_bb;
          }}
      }}
    }
  }
}

template<class FImpl>
template <class mobj>
void BaryonUtils<FImpl>::XiToSigmaEye(const PropagatorField &qq_loop,
						 const mobj &Dd_spec,
						 const mobj &Ds_spec,
						 const PropagatorField &qd_tf,
						 const PropagatorField &qs_ti,
				                 const Gamma Gamma_H,
				                 const Gamma GammaB_xi,
		                 		 const Gamma GammaB_sigma,
						 const std::string op,
						 SpinMatrixField &xts_corr)
{

  assert(Ns==4 && "Baryon code only implemented for N_spin = 4");
  assert(Nc==3 && "Baryon code only implemented for N_colour = 3");

  GridBase *grid = qs_ti.Grid();

  autoView( vcorr   , xts_corr , AcceleratorWrite);
  autoView( vq_loop , qq_loop  , AcceleratorRead);
  autoView( vd_tf   , qd_tf    , AcceleratorRead);
  autoView( vs_ti   , qs_ti    , AcceleratorRead);

  Vector<mobj> my_Dq_spec{Dd_spec,Ds_spec};
  mobj * Dq_spec_p = &my_Dq_spec[0];

  if(op == "Q1"){
    accelerator_for(ss, grid->oSites(), grid->Nsimd(), {
      auto Dq_loop = vq_loop(ss);
      auto Dd_tf   = vd_tf(ss);
      auto Ds_ti   = vs_ti(ss);
      typedef decltype(coalescedRead(vcorr[0])) spinor;
      spinor result=Zero();
      XiToSigmaQ1EyeSite(Dq_loop,Dq_spec_p[0],Dq_spec_p[1],Dd_tf,Ds_ti,Gamma_H,GammaB_xi,GammaB_sigma,result);
      coalescedWrite(vcorr[ss],result);
    }  );//end loop over lattice sites
  } else if(op == "Q2"){
    accelerator_for(ss, grid->oSites(), grid->Nsimd(), {
      auto Dq_loop = vq_loop(ss);
      auto Dd_tf   = vd_tf(ss);
      auto Ds_ti   = vs_ti(ss);
      typedef decltype(coalescedRead(vcorr[0])) spinor;
      spinor result=Zero();
      XiToSigmaQ2EyeSite(Dq_loop,Dq_spec_p[0],Dq_spec_p[1],Dd_tf,Ds_ti,Gamma_H,GammaB_xi,GammaB_sigma,result);
      coalescedWrite(vcorr[ss],result);
    }  );//end loop over lattice sites
  } else {
    assert(0 && "Weak Operator not correctly specified");
  }
}


NAMESPACE_END(Grid);
