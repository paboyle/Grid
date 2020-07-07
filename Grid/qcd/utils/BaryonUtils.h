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
//#include <Grid/Hadrons/Global.hpp>
#include <Grid/Eigen/unsupported/CXX11/Tensor>

NAMESPACE_BEGIN(Grid);

template <typename FImpl>
class BaryonUtils 
{
public:
  typedef typename FImpl::ComplexField ComplexField;
  typedef typename FImpl::FermionField FermionField;
  typedef typename FImpl::PropagatorField PropagatorField;

  typedef typename FImpl::SitePropagator pobj;
  typedef typename ComplexField::vector_object vobj;

  typedef Lattice<iSpinMatrix<typename FImpl::Simd>> SpinMatrixField;
  typedef typename SpinMatrixField::vector_object sobj;

  static const int epsilon[6][3] ;
  static const Real epsilon_sgn[6];

  private: 
  template <class mobj, class robj>
  static void baryon_site(const mobj &D1,
				 const mobj &D2,
				 const mobj &D3,
				 const Gamma GammaA_left,
				 const Gamma GammaB_left,
				 const Gamma GammaA_right,
				 const Gamma GammaB_right,
				 const int parity,
				 const bool * wick_contractions,
  				 robj &result);
  public:
  static void Wick_Contractions(std::string qi, 
                 std::string qf, 
                 bool* wick_contractions);
  static void ContractBaryons(const PropagatorField &q1_left,
				 const PropagatorField &q2_left,
				 const PropagatorField &q3_left,
				 const Gamma GammaA_left,
				 const Gamma GammaB_left,
				 const Gamma GammaA_right,
				 const Gamma GammaB_right,
				 const bool* wick_contractions,
				 const int parity,
				 ComplexField &baryon_corr);
  template <class mobj, class robj>
  static void ContractBaryons_Sliced(const mobj &D1,
				 const mobj &D2,
				 const mobj &D3,
				 const Gamma GammaA_left,
				 const Gamma GammaB_left,
				 const Gamma GammaA_right,
				 const Gamma GammaB_right,
				 const bool* wick_contractions,
				 const int parity,
				 const int nt,
				 robj &result);
  private:
  template <class mobj, class mobj2, class robj>
  static void Baryon_Gamma_3pt_Group1_Site(
           const mobj &Dq1_ti,
           const mobj2 &Dq2_spec,
           const mobj2 &Dq3_spec,
           const mobj &Dq4_tf,
                   const Gamma GammaJ,
                   const Gamma GammaBi,
                   const Gamma GammaBf,
           int wick_contraction,
           robj &result);

  template <class mobj, class mobj2, class robj>
  static void Baryon_Gamma_3pt_Group2_Site(
           const mobj2 &Dq1_spec,
           const mobj &Dq2_ti,
           const mobj2 &Dq3_spec,
           const mobj &Dq4_tf,
                   const Gamma GammaJ,
                   const Gamma GammaBi,
                   const Gamma GammaBf,
           int wick_contraction,
           robj &result);

  template <class mobj, class mobj2, class robj>
  static void Baryon_Gamma_3pt_Group3_Site(
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
  static void Baryon_Gamma_3pt(
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
  template <class mobj, class mobj2, class robj>
  static void Sigma_to_Nucleon_Q1_Eye_site(const mobj &Dq_loop,
						 const mobj2 &Du_spec,
						 const mobj &Dd_tf,
						 const mobj &Ds_ti,
				                 const Gamma Gamma_H,
				                 const Gamma GammaB_sigma,
		                 		 const Gamma GammaB_nucl,
						 robj &result);
  template <class mobj, class mobj2, class robj>
  static void Sigma_to_Nucleon_Q1_NonEye_site(const mobj &Du_ti,
						 const mobj &Du_tf,
						 const mobj2 &Du_spec,
						 const mobj &Dd_tf,
						 const mobj &Ds_ti,
				                 const Gamma Gamma_H,
				                 const Gamma GammaB_sigma,
		                 		 const Gamma GammaB_nucl,
						 robj &result);


  template <class mobj, class mobj2, class robj>
  static void Sigma_to_Nucleon_Q2_Eye_site(const mobj &Dq_loop,
						 const mobj2 &Du_spec,
						 const mobj &Dd_tf,
						 const mobj &Ds_ti,
				                 const Gamma Gamma_H,
				                 const Gamma GammaB_sigma,
		                 		 const Gamma GammaB_nucl,
						 robj &result);
  template <class mobj, class mobj2, class robj>
  static void Sigma_to_Nucleon_Q2_NonEye_site(const mobj &Du_ti,
						 const mobj &Du_tf,
						 const mobj2 &Du_spec,
						 const mobj &Dd_tf,
						 const mobj &Ds_ti,
				                 const Gamma Gamma_H,
				                 const Gamma GammaB_sigma,
		                 		 const Gamma GammaB_nucl,
						 robj &result);
  public:
  template <class mobj>
  static void Sigma_to_Nucleon_Eye(const PropagatorField &qq_loop,
				 const mobj &Du_spec,
				 const PropagatorField &qd_tf,
				 const PropagatorField &qs_ti,
				 const Gamma Gamma_H,
				 const Gamma GammaB_sigma,
				 const Gamma GammaB_nucl,
		                 const std::string op,
				 SpinMatrixField &stn_corr);
  template <class mobj>
  static void Sigma_to_Nucleon_NonEye(const PropagatorField &qq_ti,
				 const PropagatorField &qq_tf,
				 const mobj &Du_spec,
				 const PropagatorField &qd_tf,
				 const PropagatorField &qs_ti,
				 const Gamma Gamma_H,
				 const Gamma GammaB_sigma,
				 const Gamma GammaB_nucl,
		                 const std::string op,
				 SpinMatrixField &stn_corr);
};

template <class FImpl> 
const int BaryonUtils<FImpl>::epsilon[6][3] = {{0,1,2},{1,2,0},{2,0,1},{0,2,1},{2,1,0},{1,0,2}};
/*template <class FImpl> 
const Complex BaryonUtils<FImpl>::epsilon_sgn[6] = {Complex(1),
						    Complex(1),
						    Complex(1),
						    Complex(-1),
						    Complex(-1),
						    Complex(-1)};
*/
template <class FImpl> 
const Real BaryonUtils<FImpl>::epsilon_sgn[6] = {1.,1.,1.,-1.,-1.,-1.};

//This is the old version
template <class FImpl>
template <class mobj, class robj>
void BaryonUtils<FImpl>::baryon_site(const mobj &D1,
                const mobj &D2,
                const mobj &D3,
                         const Gamma GammaA_i,
                         const Gamma GammaB_i,
                         const Gamma GammaA_f,
                         const Gamma GammaB_f,
                const int parity,
                const bool * wick_contraction,
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

    for (int ie_f=0; ie_f < 6 ; ie_f++){
        int a_f = epsilon[ie_f][0]; //a
        int b_f = epsilon[ie_f][1]; //b
        int c_f = epsilon[ie_f][2]; //c
    for (int ie_i=0; ie_i < 6 ; ie_i++){
        int a_i = epsilon[ie_i][0]; //a'
        int b_i = epsilon[ie_i][1]; //b'
        int c_i = epsilon[ie_i][2]; //c'

        Real ee = epsilon_sgn[ie_f] * epsilon_sgn[ie_i];
        //This is the \delta_{456}^{123} part
        if (wick_contraction[0]){
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
        if (wick_contraction[1]){
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
        if (wick_contraction[2]){
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
        if (wick_contraction[3]){
            for (int rho=0; rho<Ns; rho++){
                auto GAf_D1_GAi_P_rr_cc = GAf_D1_GAi_P()(rho,rho)(c_f,c_i);
                for (int alpha_f=0; alpha_f<Ns; alpha_f++){
                for (int beta_i=0; beta_i<Ns; beta_i++){
                    result()()() -= ee  * GAf_D1_GAi_P_rr_cc
                                        * GBf_D2_GBi    ()(alpha_f,beta_i)(b_f,a_i)
                                        * D3            ()(alpha_f,beta_i)(a_f,b_i);
              }
        }}
        }   
        //This is the \delta_{456}^{321} part
        if (wick_contraction[4]){
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
        if (wick_contraction[5]){
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
    }}
}

/* Computes which wick contractions should be performed for a    *
 * baryon 2pt function given the initial and finals state quark  *
 * flavours.                                                     *
 * The array wick_contractions must be of length 6               */
template<class FImpl>
void BaryonUtils<FImpl>::Wick_Contractions(std::string qi, std::string qf, bool* wick_contractions) {
    const int epsilon[6][3] = {{0,1,2},{1,2,0},{2,0,1},{0,2,1},{2,1,0},{1,0,2}};
    for (int ie=0; ie < 6 ; ie++) {
        wick_contractions[ie] = (qi.size() == 3 && qf.size() == 3
                                    && qi[0] == qf[epsilon[ie][0]] 
                                    && qi[1] == qf[epsilon[ie][1]] 
                                    && qi[2] == qf[epsilon[ie][2]]);
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
						 const bool* wick_contractions,
						 const int parity,
						 ComplexField &baryon_corr)
{

  assert(Ns==4 && "Baryon code only implemented for N_spin = 4");
  assert(Nc==3 && "Baryon code only implemented for N_colour = 3");

  std::cout << "GammaA (left) " << (GammaA_left.g) <<  std::endl;
  std::cout << "GammaB (left) " << (GammaB_left.g) <<  std::endl;
  std::cout << "GammaA (right) " << (GammaA_right.g) <<  std::endl;
  std::cout << "GammaB (right) " << (GammaB_right.g) <<  std::endl;
 
  assert(parity==1 || parity == -1 && "Parity must be +1 or -1");

  GridBase *grid = q1_left.Grid();
  
  autoView(vbaryon_corr, baryon_corr,CpuWrite);
  autoView( v1 , q1_left, CpuRead);
  autoView( v2 , q2_left, CpuRead);
  autoView( v3 , q3_left, CpuRead);

  Real bytes =0.;
  bytes += grid->oSites() * (432.*sizeof(vComplex) + 126.*sizeof(int) + 36.*sizeof(Real));
  for (int ie=0; ie < 6 ; ie++){
    if(ie==0 or ie==3){
       bytes += grid->oSites() * (4.*sizeof(int) + 4752.*sizeof(vComplex)) * wick_contractions[ie];
    }
    else{
       bytes += grid->oSites() * (64.*sizeof(int) + 5184.*sizeof(vComplex)) * wick_contractions[ie];
    }
  }
  Real t=0.;
  t =-usecond();

  accelerator_for(ss, grid->oSites(), grid->Nsimd(), {
    auto D1 = v1[ss];
    auto D2 = v2[ss];
    auto D3 = v3[ss];
    vobj result=Zero();
    baryon_site(D1,D2,D3,GammaA_left,GammaB_left,GammaA_right,GammaB_right,parity,wick_contractions,result);
    vbaryon_corr[ss] = result; 
  }  );//end loop over lattice sites

  t += usecond();

  std::cout << std::setw(10) << bytes/t*1.0e6/1024/1024/1024 << " GB/s " << std::endl;

}

/* The array wick_contractions must be of length 6. The order     * 
 * corresponds to the to that shown in the Hadrons documentation  *
 * at https://aportelli.github.io/Hadrons-doc/#/mcontraction      *
 * This can also be computed from the quark flavours using the    *
 * Wick_Contractions function above                               */
template <class FImpl>
template <class mobj, class robj>
void BaryonUtils<FImpl>::ContractBaryons_Sliced(const mobj &D1,
						 const mobj &D2,
						 const mobj &D3,
				                 const Gamma GammaA_left,
				                 const Gamma GammaB_left,
				                 const Gamma GammaA_right,
		                 		 const Gamma GammaB_right,
						 const bool* wick_contractions,
						 const int parity,
						 const int nt,
						 robj &result)
{

  assert(Ns==4 && "Baryon code only implemented for N_spin = 4");
  assert(Nc==3 && "Baryon code only implemented for N_colour = 3");

  std::cout << "GammaA (left) " << (GammaA_left.g) <<  std::endl;
  std::cout << "GammaB (left) " << (GammaB_left.g) <<  std::endl;
  std::cout << "GammaA (right) " << (GammaA_right.g) <<  std::endl;
  std::cout << "GammaB (right) " << (GammaB_right.g) <<  std::endl;
 
  assert(parity==1 || parity == -1 && "Parity must be +1 or -1");

  for (int t=0; t<nt; t++) {
    baryon_site(D1[t],D2[t],D3[t],GammaA_left,GammaB_left,GammaA_right,GammaB_right,parity,wick_contractions,result[t]);
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
template <class mobj, class mobj2, class robj>
void BaryonUtils<FImpl>::Baryon_Gamma_3pt_Group1_Site(
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

    auto adjD4_g_D1     = g5 * adj(Dq4_tf) * g5 * GammaJ * Dq1_ti;
    auto Gf_adjD4_g_D1  = GammaBf * adjD4_g_D1;
    auto D2_Gi          = Dq2_spec * GammaBi;
    auto Gf_D2_Gi       = GammaBf * D2_Gi;
    auto Gf_D3          = GammaBf * Dq3_spec;

    int a_f, b_f, c_f;
    int a_i, b_i, c_i;

    Real ee;

    for (int ie_f=0; ie_f < 6 ; ie_f++){
        a_f = epsilon[ie_f][0]; //a
        b_f = epsilon[ie_f][1]; //b
        c_f = epsilon[ie_f][2]; //c
    for (int ie_i=0; ie_i < 6 ; ie_i++){
        a_i = epsilon[ie_i][0]; //a'
        b_i = epsilon[ie_i][1]; //b'
        c_i = epsilon[ie_i][2]; //c'

        ee = epsilon_sgn[ie_f] * epsilon_sgn[ie_i];

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
            }}
        }}
    }}
}

/* Dq1_spec is a quark line from t_i to t_f
 * Dq2_ti is a quark line from t_i to t_J
 * Dq3_spec is a quark line from t_i to t_f
 * Dq4_tf is a quark line from t_f to t_J */
template<class FImpl>
template <class mobj, class mobj2, class robj>
void BaryonUtils<FImpl>::Baryon_Gamma_3pt_Group2_Site(
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

    int a_f, b_f, c_f;
    int a_i, b_i, c_i;

    Real ee;

    for (int ie_f=0; ie_f < 6 ; ie_f++){
        a_f = epsilon[ie_f][0]; //a
        b_f = epsilon[ie_f][1]; //b
        c_f = epsilon[ie_f][2]; //c
    for (int ie_i=0; ie_i < 6 ; ie_i++){
        a_i = epsilon[ie_i][0]; //a'
        b_i = epsilon[ie_i][1]; //b'
        c_i = epsilon[ie_i][2]; //c'

        ee = epsilon_sgn[ie_f] * epsilon_sgn[ie_i];

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
            }}
        }}
    }}
}

/* Dq1_spec is a quark line from t_i to t_f
 * Dq2_spec is a quark line from t_i to t_f
 * Dq3_ti is a quark line from t_i to t_J
 * Dq4_tf is a quark line from t_f to t_J */
template<class FImpl>
template <class mobj, class mobj2, class robj>
void BaryonUtils<FImpl>::Baryon_Gamma_3pt_Group3_Site(
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

    int a_f, b_f, c_f;
    int a_i, b_i, c_i;

    Real ee;

    for (int ie_f=0; ie_f < 6 ; ie_f++){
        a_f = epsilon[ie_f][0]; //a
        b_f = epsilon[ie_f][1]; //b
        c_f = epsilon[ie_f][2]; //c
    for (int ie_i=0; ie_i < 6 ; ie_i++){
        a_i = epsilon[ie_i][0]; //a'
        b_i = epsilon[ie_i][1]; //b'
        c_i = epsilon[ie_i][2]; //c'

        ee = epsilon_sgn[ie_f] * epsilon_sgn[ie_i];

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
            }}
        }}
    }}
}

/* The group indicates which inital state quarks the current is  * 
 * connected to. It must be in the range 1-3.                    *
 * The wick_contraction must be in the range 1-6 correspond to   *
 * the contractions given in the Hadrons documentation at        *
 * https://aportelli.github.io/Hadrons-doc/#/mcontraction        */
template<class FImpl>
template <class mobj>
void BaryonUtils<FImpl>::Baryon_Gamma_3pt(
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
    GridBase *grid = q_tf.Grid();

    autoView( vcorr, stn_corr, CpuWrite);
    autoView( vq_ti , q_ti, CpuRead);
    autoView( vq_tf , q_tf, CpuRead);

    if (group == 1) {
        accelerator_for(ss, grid->oSites(), grid->Nsimd(), {
            auto Dq_ti = vq_ti[ss];
            auto Dq_tf = vq_tf[ss];
            sobj result=Zero();
            Baryon_Gamma_3pt_Group1_Site(Dq_ti,Dq_spec1,Dq_spec2,Dq_tf,GammaJ,GammaBi,GammaBf,wick_contraction,result);
            vcorr[ss] += result; 
        });//end loop over lattice sites
    } else if (group == 2) {
        accelerator_for(ss, grid->oSites(), grid->Nsimd(), {
            auto Dq_ti = vq_ti[ss];
            auto Dq_tf = vq_tf[ss];
            sobj result=Zero();
            Baryon_Gamma_3pt_Group2_Site(Dq_spec1,Dq_ti,Dq_spec2,Dq_tf,GammaJ,GammaBi,GammaBf,wick_contraction,result);
            vcorr[ss] += result; 
        });//end loop over lattice sites
    } else if (group == 3) {
        accelerator_for(ss, grid->oSites(), grid->Nsimd(), {
            auto Dq_ti = vq_ti[ss];
            auto Dq_tf = vq_tf[ss];
            sobj result=Zero();
            Baryon_Gamma_3pt_Group3_Site(Dq_spec1,Dq_spec2,Dq_ti,Dq_tf,GammaJ,GammaBi,GammaBf,wick_contraction,result);

            vcorr[ss] += result; 
        });//end loop over lattice sites
    }
}


/***********************************************************************
 * End of BaryonGamma3pt-function code.                                *
 *																	   *
 * The following code is for Sigma -> N rare hypeon decays             *
 **********************************************************************/

/* Dq_loop is a quark line from t_H to t_H
 * Du_spec is a quark line from t_i to t_f
 * Dd_tf is a quark line from t_f to t_H
 * Ds_ti is a quark line from t_i to t_H */
template <class FImpl>
template <class mobj, class mobj2, class robj>
void BaryonUtils<FImpl>::Sigma_to_Nucleon_Q1_Eye_site(const mobj &Dq_loop,
						 const mobj2 &Du_spec,
						 const mobj &Dd_tf,
						 const mobj &Ds_ti,
				                 const Gamma Gamma_H,
				                 const Gamma GammaB_sigma,
		                 		 const Gamma GammaB_nucl,
						 robj &result)
{

  Gamma g5(Gamma::Algebra::Gamma5); 

  auto DuG = Du_spec * GammaB_nucl;
  // Gamma^B * Ds * \gamma_\mu^L * (\gamma_5 * Dd^\dagger * \gamma_5)
  auto GDsGDd = GammaB_sigma * Ds_ti * Gamma_H * g5 * adj(Dd_tf) * g5;
  // Dq_loop * \gamma_\mu^L
  auto DqG = Dq_loop * Gamma_H;

  for (int ie_n=0; ie_n < 6 ; ie_n++){
    int a_n = epsilon[ie_n][0]; //a
    int b_n = epsilon[ie_n][1]; //b
    int c_n = epsilon[ie_n][2]; //c
    for (int ie_s=0; ie_s < 6 ; ie_s++){
      int a_s = epsilon[ie_s][0]; //a'
      int b_s = epsilon[ie_s][1]; //b'
      int c_s = epsilon[ie_s][2]; //c'
      for (int alpha_s=0; alpha_s<Ns; alpha_s++){
      for (int beta_n=0; beta_n<Ns; beta_n++){
        auto GDsGDd_ab_bb = GDsGDd()(alpha_s,beta_n)(b_s,b_n);
        for (int tau2=0; tau2<Ns; tau2++){
        for (int j=0; j<Nc; j++){
          auto DqG_tt_jj = DqG()(tau2,tau2)(j,j);
          auto ee_GDGDDG = epsilon_sgn[ie_n] * epsilon_sgn[ie_s] * GDsGDd_ab_bb * DqG_tt_jj;
          for (int gamma_s=0; gamma_s<Ns; gamma_s++){
          for (int gamma_n=0; gamma_n<Ns; gamma_n++){
            result()(gamma_s,gamma_n)() += ee_GDGDDG * DuG()(alpha_s, beta_n)(a_s,a_n) * Du_spec()(gamma_s,gamma_n)(c_s,c_n);
            result()(gamma_s,gamma_n)() -= ee_GDGDDG * DuG()(gamma_s, beta_n)(c_s,a_n) * Du_spec()(alpha_s,gamma_n)(a_s,c_n);
          }}
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
template <class mobj, class mobj2, class robj>
void BaryonUtils<FImpl>::Sigma_to_Nucleon_Q1_NonEye_site(const mobj &Du_ti,
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

  auto DuG = Du_spec * GammaB_nucl;
  auto adjDu = g5 * adj(Du_tf) * g5;
  auto adjDuG = adjDu * GammaB_nucl;
  // Gamma^B * Ds * \gamma_\mu^L * (\gamma_5 * Dd^\dagger * \gamma_5)
  auto GDsGDd = GammaB_sigma * Ds_ti * Gamma_H * g5 * adj(Dd_tf) * g5;
  // Dq_loop * \gamma_\mu^L
  auto DuGH = Du_ti * Gamma_H;

  for (int ie_n=0; ie_n < 6 ; ie_n++){
    int a_n = epsilon[ie_n][0]; //a
    int b_n = epsilon[ie_n][1]; //b
    int c_n = epsilon[ie_n][2]; //c
    for (int ie_s=0; ie_s < 6 ; ie_s++){
      int a_s = epsilon[ie_s][0]; //a'
      int b_s = epsilon[ie_s][1]; //b'
      int c_s = epsilon[ie_s][2]; //c'
      for (int alpha_s=0; alpha_s<Ns; alpha_s++){
      for (int beta_n=0; beta_n<Ns; beta_n++){
        auto GDsGDd_ab_bb = GDsGDd()(alpha_s,beta_n)(b_s,b_n);
        for (int tau2=0; tau2<Ns; tau2++){
        for (int j=0; j<Nc; j++){
          auto DuGH_at_aj = DuGH()(alpha_s,tau2)(a_s,j);
          auto ee_GDGDDG_a = epsilon_sgn[ie_n] * epsilon_sgn[ie_s] * GDsGDd_ab_bb * DuGH_at_aj;
          for (int gamma_s=0; gamma_s<Ns; gamma_s++){
            auto DuGH_gt_cj = DuGH()(gamma_s,tau2)(c_s,j);
            auto ee_GDGDDG_c = epsilon_sgn[ie_n] * epsilon_sgn[ie_s] * GDsGDd_ab_bb * DuGH_gt_cj;
            for (int gamma_n=0; gamma_n<Ns; gamma_n++){
              result()(gamma_s,gamma_n)() += ee_GDGDDG_a * DuG()(gamma_s, beta_n)(c_s,a_n) * adjDu()(tau2,gamma_n)(j,c_n);
              result()(gamma_s,gamma_n)() += ee_GDGDDG_c * adjDuG()(tau2, beta_n)(j,a_n) * Du_spec()(alpha_s,gamma_n)(a_s,c_n);
              result()(gamma_s,gamma_n)() -= ee_GDGDDG_a * adjDuG()(tau2, beta_n)(j,a_n) * Du_spec()(gamma_s,gamma_n)(c_s,c_n);
              result()(gamma_s,gamma_n)() -= ee_GDGDDG_c * DuG()(alpha_s, beta_n)(a_s,a_n) * adjDu()(tau2,gamma_n)(j,c_n);
            }
	  }
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
template <class mobj, class mobj2, class robj>
void BaryonUtils<FImpl>::Sigma_to_Nucleon_Q2_Eye_site(const mobj &Dq_loop,
						 const mobj2 &Du_spec,
						 const mobj &Dd_tf,
						 const mobj &Ds_ti,
				                 const Gamma Gamma_H,
				                 const Gamma GammaB_sigma,
		                 		 const Gamma GammaB_nucl,
						 robj &result)
{

  Gamma g5(Gamma::Algebra::Gamma5); 

  auto DuG = Du_spec * GammaB_nucl;
  // Gamma^B * Ds * \gamma_\mu^L
  auto GDsG = GammaB_sigma * Ds_ti * Gamma_H;
  // Dq_loop * \gamma_\mu^L * (\gamma_5 * Dd^\dagger * \gamma_5)
  auto DqGDd = Dq_loop * Gamma_H * g5 * adj(Dd_tf) * g5;

  for (int ie_n=0; ie_n < 6 ; ie_n++){
    int a_n = epsilon[ie_n][0]; //a
    int b_n = epsilon[ie_n][1]; //b
    int c_n = epsilon[ie_n][2]; //c
    for (int ie_s=0; ie_s < 6 ; ie_s++){
      int a_s = epsilon[ie_s][0]; //a'
      int b_s = epsilon[ie_s][1]; //b'
      int c_s = epsilon[ie_s][2]; //c'
      for (int alpha_s=0; alpha_s<Ns; alpha_s++){
      for (int tau=0; tau<Ns; tau++){
      for (int i=0; i<Nc; i++){
	auto GDsG_at_bi = GDsG()(alpha_s,tau)(b_s,i);
        for (int beta_n=0; beta_n<Ns; beta_n++){
          auto DqGDd_tb_ib = DqGDd()(tau,beta_n)(i,b_n);
	  auto ee_GDGDGD = epsilon_sgn[ie_n] * epsilon_sgn[ie_s] * GDsG_at_bi * DqGDd_tb_ib;
          for (int gamma_s=0; gamma_s<Ns; gamma_s++){
          for (int gamma_n=0; gamma_n<Ns; gamma_n++){
            result()(gamma_s,gamma_n)() -= ee_GDGDGD * DuG()(alpha_s, beta_n)(a_s,a_n) * Du_spec()(gamma_s,gamma_n)(c_s,c_n);
            result()(gamma_s,gamma_n)() += ee_GDGDGD * DuG()(gamma_s, beta_n)(c_s,a_n) * Du_spec()(alpha_s,gamma_n)(a_s,c_n);
          }}
	}
      }}}
    }
  }
}

/* Du_ti is a quark line from t_i to t_H
 * Du_tf is a quark line from t_f to t_H
 * Du_spec is a quark line from t_i to t_f
 * Dd_tf is a quark line from t_f to t_H
 * Ds_ti is a quark line from t_i to t_H */
template <class FImpl>
template <class mobj, class mobj2, class robj>
void BaryonUtils<FImpl>::Sigma_to_Nucleon_Q2_NonEye_site(const mobj &Du_ti,
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

  auto DuG = Du_spec * GammaB_nucl;
  auto adjDu = g5 * adj(Du_tf) * g5;
  auto adjDuG = adjDu * GammaB_nucl;
  // Gamma^B * Ds * \gamma_\mu^L
  auto GDsG = GammaB_sigma * Ds_ti * Gamma_H;
  // Du * \gamma_\mu^L * (\gamma_5 * Dd^\dagger * \gamma_5)
  auto DuGDd = Du_ti * Gamma_H * g5 * adj(Dd_tf) * g5;

  for (int ie_n=0; ie_n < 6 ; ie_n++){
    int a_n = epsilon[ie_n][0]; //a
    int b_n = epsilon[ie_n][1]; //b
    int c_n = epsilon[ie_n][2]; //c
    for (int ie_s=0; ie_s < 6 ; ie_s++){
      int a_s = epsilon[ie_s][0]; //a'
      int b_s = epsilon[ie_s][1]; //b'
      int c_s = epsilon[ie_s][2]; //c'
      for (int alpha_s=0; alpha_s<Ns; alpha_s++){
      for (int tau=0; tau<Ns; tau++){
      for (int i=0; i<Nc; i++){
	auto GDsG_at_bi = GDsG()(alpha_s,tau)(b_s,i);
        for (int beta_n=0; beta_n<Ns; beta_n++){
          auto DuGDd_ab_ab = DuGDd()(alpha_s,beta_n)(a_s,b_n);
	  auto ee_GDGDGD_a = epsilon_sgn[ie_n] * epsilon_sgn[ie_s] * GDsG_at_bi * DuGDd_ab_ab;
          for (int gamma_s=0; gamma_s<Ns; gamma_s++){
            auto DuGDd_gb_cb = DuGDd()(gamma_s,beta_n)(c_s,b_n);
	    auto ee_GDGDGD_c = epsilon_sgn[ie_n] * epsilon_sgn[ie_s] * GDsG_at_bi * DuGDd_gb_cb;
            for (int gamma_n=0; gamma_n<Ns; gamma_n++){
              result()(gamma_s,gamma_n)() -= ee_GDGDGD_a * DuG()(gamma_s, beta_n)(c_s,a_n) * adjDu()(tau,gamma_n)(i,c_n);
              result()(gamma_s,gamma_n)() -= ee_GDGDGD_c * adjDuG()(tau, beta_n)(i,a_n) * Du_spec()(alpha_s,gamma_n)(a_s,c_n);
              result()(gamma_s,gamma_n)() += ee_GDGDGD_a * adjDuG()(tau, beta_n)(i,a_n) * Du_spec()(gamma_s,gamma_n)(c_s,c_n);
              result()(gamma_s,gamma_n)() += ee_GDGDGD_c * DuG()(alpha_s, beta_n)(a_s,a_n) * adjDu()(tau,gamma_n)(i,c_n);
            }
	  }
	}
      }}}
    }
  }
}

template<class FImpl>
template <class mobj>
void BaryonUtils<FImpl>::Sigma_to_Nucleon_Eye(const PropagatorField &qq_loop,
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

  autoView( vcorr, stn_corr, CpuWrite);
  autoView( vq_loop , qq_loop, CpuRead);
  autoView( vd_tf , qd_tf, CpuRead);
  autoView( vs_ti , qs_ti, CpuRead);

  accelerator_for(ss, grid->oSites(), grid->Nsimd(), {
    auto Dq_loop = vq_loop[ss];
    auto Dd_tf = vd_tf[ss];
    auto Ds_ti = vs_ti[ss];
    sobj result=Zero();
    if(op == "Q1"){
      Sigma_to_Nucleon_Q1_Eye_site(Dq_loop,Du_spec,Dd_tf,Ds_ti,Gamma_H,GammaB_sigma,GammaB_nucl,result);
    } else if(op == "Q2"){
      Sigma_to_Nucleon_Q2_Eye_site(Dq_loop,Du_spec,Dd_tf,Ds_ti,Gamma_H,GammaB_sigma,GammaB_nucl,result);
    } else {
      assert(0 && "Weak Operator not correctly specified");
    }
      vcorr[ss] = result; 
  }  );//end loop over lattice sites
}

template<class FImpl>
template <class mobj>
void BaryonUtils<FImpl>::Sigma_to_Nucleon_NonEye(const PropagatorField &qq_ti,
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

  autoView( vcorr , stn_corr, CpuWrite);
  autoView( vq_ti , qq_ti, CpuRead);
  autoView( vq_tf , qq_tf, CpuRead);
  autoView( vd_tf , qd_tf, CpuRead);
  autoView( vs_ti , qs_ti, CpuRead);
 // accelerator_for(ss, grid->oSites(), grid->Nsimd(), {
  thread_for(ss,grid->oSites(),{
    auto Dq_ti = vq_ti[ss];
    auto Dq_tf = vq_tf[ss];
    auto Dd_tf = vd_tf[ss];
    auto Ds_ti = vs_ti[ss];
    sobj result=Zero();
    if(op == "Q1"){
      Sigma_to_Nucleon_Q1_NonEye_site(Dq_ti,Dq_tf,Du_spec,Dd_tf,Ds_ti,Gamma_H,GammaB_sigma,GammaB_nucl,result);
    } else if(op == "Q2"){
      Sigma_to_Nucleon_Q2_NonEye_site(Dq_ti,Dq_tf,Du_spec,Dd_tf,Ds_ti,Gamma_H,GammaB_sigma,GammaB_nucl,result);
    } else {
      assert(0 && "Weak Operator not correctly specified");
    }
      vcorr[ss] = result; 
  }  );//end loop over lattice sites
}

NAMESPACE_END(Grid);
