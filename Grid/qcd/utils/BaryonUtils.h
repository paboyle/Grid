/*************************************************************************************
 
 Grid physics library, www.github.com/paboyle/Grid
 
 Source file: ./lib/qcd/utils/BaryonUtils.h
 
 Copyright (C) 2019
 
 Author: Felix Erben <felix.erben@ed.ac.uk>

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
  static const Complex epsilon_sgn[6];

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
				 const int * wick_contractions,
  				 robj &result);
  public:
  static void ContractBaryons(const PropagatorField &q1_left,
				 const PropagatorField &q2_left,
				 const PropagatorField &q3_left,
				 const Gamma GammaA_left,
				 const Gamma GammaB_left,
				 const Gamma GammaA_right,
				 const Gamma GammaB_right,
				 const char * quarks_left,
				 const char * quarks_right,
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
				 const char * quarks_left,
				 const char * quarks_right,
				 const int parity,
				 robj &result);
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
  template <class mobj, class mobj2, class robj>
  static void Xi_to_Sigma_Q1_Eye_site(const mobj &Dq_loop,
						 const mobj2 &Dd_spec,
						 const mobj &Dd_tf,
						 const mobj2 &Ds_spec,
						 const mobj &Ds_ti,
				                 const Gamma Gamma_H,
				                 const Gamma GammaB_xi,
		                 		 const Gamma GammaB_sigma,
						 robj &result);
  template <class mobj, class mobj2, class robj>
  static void Xi_to_Sigma_Q2_Eye_site(const mobj &Dq_loop,
						 const mobj2 &Dd_spec,
						 const mobj &Dd_tf,
						 const mobj2 &Ds_spec,
						 const mobj &Ds_ti,
				                 const Gamma Gamma_H,
				                 const Gamma GammaB_xi,
		                 		 const Gamma GammaB_sigma,
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
  template <class mobj>
  static void Xi_to_Sigma_Eye(const PropagatorField &qq_loop,
				 const mobj &Dd_spec,
				 const PropagatorField &qd_tf,
				 const mobj &Ds_spec,
				 const PropagatorField &qs_ti,
				 const Gamma Gamma_H,
				 const Gamma GammaB_xi,
				 const Gamma GammaB_sigma,
		                 const std::string op,
				 SpinMatrixField &stn_corr);
};

template <class FImpl> 
const int BaryonUtils<FImpl>::epsilon[6][3] = {{0,1,2},{1,2,0},{2,0,1},{0,2,1},{2,1,0},{1,0,2}};
template <class FImpl> 
const Complex BaryonUtils<FImpl>::epsilon_sgn[6] = {Complex(1),
						    Complex(1),
						    Complex(1),
						    Complex(-1),
						    Complex(-1),
						    Complex(-1)};

template <class FImpl>
template <class mobj, class robj>
void BaryonUtils<FImpl>::baryon_site(const mobj &D1,
						 const mobj &D2,
						 const mobj &D3,
				                 const Gamma GammaA_left,
				                 const Gamma GammaB_left,
				                 const Gamma GammaA_right,
		                 		 const Gamma GammaB_right,
						 const int parity,
						 const int * wick_contraction,
						 robj &result)
{

  Gamma g4(Gamma::Algebra::GammaT); //needed for parity P_\pm = 0.5*(1 \pm \gamma_4)

    auto gD1a = GammaA_left * GammaA_right * D1;
    auto gD1b = GammaA_left * g4 * GammaA_right * D1;
    auto pD1 = 0.5* (gD1a + (double)parity * gD1b);
    auto gD3 = GammaB_right * D3;

    for (int ie_left=0; ie_left < 6 ; ie_left++){
      int a_left = epsilon[ie_left][0]; //a
      int b_left = epsilon[ie_left][1]; //b
      int c_left = epsilon[ie_left][2]; //c
      for (int ie_right=0; ie_right < 6 ; ie_right++){
        int a_right = epsilon[ie_right][0]; //a'
        int b_right = epsilon[ie_right][1]; //b'
        int c_right = epsilon[ie_right][2]; //c'
        //This is the \delta_{456}^{123} part
	if (wick_contraction[0]){
          auto D2g = D2 * GammaB_left;
	  for (int alpha_right=0; alpha_right<Ns; alpha_right++){
	  for (int beta_left=0; beta_left<Ns; beta_left++){
	  for (int gamma_left=0; gamma_left<Ns; gamma_left++){
	    result()()() += epsilon_sgn[ie_left] * epsilon_sgn[ie_right] * pD1()(gamma_left,gamma_left)(c_right,c_left)*D2g()(alpha_right,beta_left)(a_right,a_left)*gD3()(alpha_right,beta_left)(b_right,b_left);
          }}}
  	}	  
        //This is the \delta_{456}^{231} part
	if (wick_contraction[1]){
          auto pD1g = pD1 * GammaB_left;
	  for (int alpha_right=0; alpha_right<Ns; alpha_right++){
	  for (int beta_left=0; beta_left<Ns; beta_left++){
	  for (int gamma_left=0; gamma_left<Ns; gamma_left++){
	    result()()() += epsilon_sgn[ie_left] * epsilon_sgn[ie_right] * pD1g()(gamma_left,beta_left)(c_right,a_left)*D2()(alpha_right,beta_left)(a_right,b_left)*gD3()(alpha_right,gamma_left)(b_right,c_left);
          }}}
        }	  
        //This is the \delta_{456}^{312} part
	if (wick_contraction[2]){
          auto gD3g = gD3 * GammaB_left;
	  for (int alpha_right=0; alpha_right<Ns; alpha_right++){
	  for (int beta_left=0; beta_left<Ns; beta_left++){
	  for (int gamma_left=0; gamma_left<Ns; gamma_left++){
	    result()()() += epsilon_sgn[ie_left] * epsilon_sgn[ie_right] * pD1()(gamma_left,beta_left)(c_right,b_left)*D2()(alpha_right,gamma_left)(a_right,c_left)*gD3g()(alpha_right,beta_left)(b_right,a_left);
          }}}
        }	  
        //This is the \delta_{456}^{132} part
	if (wick_contraction[3]){
          auto gD3g = gD3 * GammaB_left;
	  for (int alpha_right=0; alpha_right<Ns; alpha_right++){
	  for (int beta_left=0; beta_left<Ns; beta_left++){
	  for (int gamma_left=0; gamma_left<Ns; gamma_left++){
	    result()()() -= epsilon_sgn[ie_left] * epsilon_sgn[ie_right] * pD1()(gamma_left,gamma_left)(c_right,c_left)*D2()(alpha_right,beta_left)(a_right,b_left)*gD3g()(alpha_right,beta_left)(b_right,a_left);
          }}}
        }	  
        //This is the \delta_{456}^{321} part
	if (wick_contraction[4]){
          auto D2g = D2 * GammaB_left;
	  for (int alpha_right=0; alpha_right<Ns; alpha_right++){
	  for (int beta_left=0; beta_left<Ns; beta_left++){
	  for (int gamma_left=0; gamma_left<Ns; gamma_left++){
	    result()()() -= epsilon_sgn[ie_left] * epsilon_sgn[ie_right] * pD1()(gamma_left,beta_left)(c_right,b_left)*D2g()(alpha_right,beta_left)(a_right,a_left)*gD3()(alpha_right,gamma_left)(b_right,c_left);
          }}}
        }	  
        //This is the \delta_{456}^{213} part
	if (wick_contraction[5]){
          auto pD1g = pD1 * GammaB_left;
	  for (int alpha_right=0; alpha_right<Ns; alpha_right++){
	  for (int beta_left=0; beta_left<Ns; beta_left++){
	  for (int gamma_left=0; gamma_left<Ns; gamma_left++){
	    result()()() -= epsilon_sgn[ie_left] * epsilon_sgn[ie_right] * pD1g()(gamma_left,beta_left)(c_right,a_left)*D2()(alpha_right,gamma_left)(a_right,c_left)*gD3()(alpha_right,beta_left)(b_right,b_left);
          }}}
        }	  
      }
    }
}

template<class FImpl>
void BaryonUtils<FImpl>::ContractBaryons(const PropagatorField &q1_left,
						 const PropagatorField &q2_left,
						 const PropagatorField &q3_left,
				                 const Gamma GammaA_left,
				                 const Gamma GammaB_left,
				                 const Gamma GammaA_right,
		                 		 const Gamma GammaB_right,
						 const char * quarks_left,
						 const char * quarks_right,
						 const int parity,
						 ComplexField &baryon_corr)
{
  std::cout << "Contraction <" << quarks_right[0] << quarks_right[1] << quarks_right[2] << "|" << quarks_left[0] << quarks_left[1] << quarks_left[2] << ">" << std::endl;
    std::cout << "GammaA (left) " << (GammaA_left.g) <<  std::endl;
    std::cout << "GammaB (left) " << (GammaB_left.g) <<  std::endl;
    std::cout << "GammaA (right) " << (GammaA_right.g) <<  std::endl;
    std::cout << "GammaB (right) " << (GammaB_right.g) <<  std::endl;
 
  assert(parity==1 || parity == -1 && "Parity must be +1 or -1");

  GridBase *grid = q1_left.Grid();

  int wick_contraction[6];
  for (int ie=0; ie < 6 ; ie++)
    wick_contraction[ie] = (quarks_left[0] == quarks_right[epsilon[ie][0]] && quarks_left[1] == quarks_right[epsilon[ie][1]] && quarks_left[2] == quarks_right[epsilon[ie][2]]) ? 1 : 0;

  auto vbaryon_corr= baryon_corr.View();
  auto v1 = q1_left.View();
  auto v2 = q2_left.View();
  auto v3 = q3_left.View();

 // accelerator_for(ss, grid->oSites(), grid->Nsimd(), {
  thread_for(ss,grid->oSites(),{
  //for(int ss=0; ss < grid->oSites(); ss++){

    auto D1 = v1[ss];
    auto D2 = v2[ss];
    auto D3 = v3[ss];

    vobj result=Zero();
    baryon_site(D1,D2,D3,GammaA_left,GammaB_left,GammaA_right,GammaB_right,parity,wick_contraction,result);
    vbaryon_corr[ss] = result; 
  }  );//end loop over lattice sites
}
template <class FImpl>
template <class mobj, class robj>
void BaryonUtils<FImpl>::ContractBaryons_Sliced(const mobj &D1,
						 const mobj &D2,
						 const mobj &D3,
				                 const Gamma GammaA_left,
				                 const Gamma GammaB_left,
				                 const Gamma GammaA_right,
		                 		 const Gamma GammaB_right,
						 const char * quarks_left,
						 const char * quarks_right,
						 const int parity,
						 robj &result)
{
  std::cout << "Contraction <" << quarks_right[0] << quarks_right[1] << quarks_right[2] << "|" << quarks_left[0] << quarks_left[1] << quarks_left[2] << ">" << std::endl;
    std::cout << "GammaA (left) " << (GammaA_left.g) <<  std::endl;
    std::cout << "GammaB (left) " << (GammaB_left.g) <<  std::endl;
    std::cout << "GammaA (right) " << (GammaA_right.g) <<  std::endl;
    std::cout << "GammaB (right) " << (GammaB_right.g) <<  std::endl;
 
  assert(parity==1 || parity == -1 && "Parity must be +1 or -1");

  int wick_contraction[6];
  for (int ie=0; ie < 6 ; ie++)
    wick_contraction[ie] = (quarks_left[0] == quarks_right[epsilon[ie][0]] && quarks_left[1] == quarks_right[epsilon[ie][1]] && quarks_left[2] == quarks_right[epsilon[ie][2]]) ? 1 : 0;

     result=Zero();
     baryon_site(D1,D2,D3,GammaA_left,GammaB_left,GammaA_right,GammaB_right,parity,wick_contraction,result);
}

/***********************************************************************
 * End of Baryon 2pt-function code.                                    *
 *                                                                     *
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
  GridBase *grid = qs_ti.Grid();

  auto vcorr= stn_corr.View();
  auto vq_loop = qq_loop.View();
  auto vd_tf = qd_tf.View();
  auto vs_ti = qs_ti.View();

 // accelerator_for(ss, grid->oSites(), grid->Nsimd(), {
  thread_for(ss,grid->oSites(),{
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
  GridBase *grid = qs_ti.Grid();

  auto vcorr= stn_corr.View();
  auto vq_ti = qq_ti.View();
  auto vq_tf = qq_tf.View();
  auto vd_tf = qd_tf.View();
  auto vs_ti = qs_ti.View();

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

/***********************************************************************
 * The following code is for Xi -> Sigma rare hypeon decays            *
 **********************************************************************/

/* Dq_loop is a quark line from t_H to t_H
 * Dd_spec is a quark line from t_i to t_f
 * Dd_tf is a quark line from t_f to t_H
 * Ds_spec is a quark line from t_i to t_f
 * Ds_ti is a quark line from t_i to t_H */
template <class FImpl>
template <class mobj, class mobj2, class robj>
void BaryonUtils<FImpl>::Xi_to_Sigma_Q1_Eye_site(const mobj &Dq_loop,
						 const mobj2 &Dd_spec,
						 const mobj &Dd_tf,
						 const mobj2 &Ds_spec,
						 const mobj &Ds_ti,
				                 const Gamma Gamma_H,
				                 const Gamma GammaB_xi,
		                 		 const Gamma GammaB_sigma,
						 robj &result)
{

  Gamma g5(Gamma::Algebra::Gamma5); 

  auto DdG = Dd_spec * GammaB_sigma;
  auto GDs = GammaB_xi * Ds_spec;
  //  Ds * \gamma_\mu^L * Dd
  auto DsGDd = Ds_ti * Gamma_H * Dd_tf;
  // Ds * \gamma_\mu^L * Dd * \Gamma^B
  auto GDsGDd = GammaB_sigma * DsGDd;
  // Ds * \gamma_\mu^L * Dd * \Gamma^B
  auto DsGDdG = DsGDd * GammaB_sigma;
  // Gamma^B * Ds * \gamma_\mu^L * Dd * \Gamma^B
  auto GDsGDdG = GammaB_xi * DsGDdG;
  // \gamma_\mu^L * Dq_loop 
  auto GDq = Gamma_H * Dq_loop;

  for (int ie_s=0; ie_s < 6 ; ie_s++){
    int a_s = epsilon[ie_s][0]; //a
    int b_s = epsilon[ie_s][1]; //b
    int c_s = epsilon[ie_s][2]; //c
    for (int ie_x=0; ie_x < 6 ; ie_x++){
      int a_x = epsilon[ie_x][0]; //a'
      int b_x = epsilon[ie_x][1]; //b'
      int c_x = epsilon[ie_x][2]; //c'
      for (int sigma2=0; sigma2<Ns; sigma2++){
      for (int j=0; j<Nc; j++){
        auto GDq_jj_ss=GDq()(sigma2,sigma2)(j,j);
        auto ee_GD = epsilon_sgn[ie_s] * epsilon_sgn[ie_x] * GDq_jj_ss;
        for (int alpha_x=0; alpha_x<Ns; alpha_x++){
        for (int beta_s=0; beta_s<Ns; beta_s++){
          auto GDsGDd_ab_bb = GDsGDd()(alpha_x,beta_s)(b_x,b_s);
          for (int gamma_x=0; gamma_x<Ns; gamma_x++){
          for (int gamma_s=0; gamma_s<Ns; gamma_s++){
            result()(gamma_x,gamma_s)() -= ee_GD * Dd_spec()(gamma_x, gamma_s)(c_x,c_s) * GDsGDdG()(alpha_x, beta_s)(b_x,a_s) * Ds_spec()(alpha_x, beta_s)(a_x, b_s);
            result()(gamma_x,gamma_s)() += ee_GD * DdG()(gamma_x, beta_s)(c_x,a_s) * GDsGDd()(alpha_x, gamma_s)(b_x,c_s) * Ds_spec()(alpha_x, beta_s)(a_x, b_s);
            result()(gamma_x,gamma_s)() += ee_GD * Dd_spec()(gamma_x, gamma_s)(c_x,c_s) * DsGDdG()(alpha_x, beta_s)(a_x,a_s) * GDs()(alpha_x, beta_s)(b_x, b_s);
            result()(gamma_x,gamma_s)() -= ee_GD * DdG()(gamma_x, beta_s)(c_x,a_s) * DsGDd()(alpha_x, gamma_s)(a_x,c_s) * GDs()(alpha_x, beta_s)(b_x, b_s);
          }}
	}}
      }}
    }
  }
}

/* Dq_loop is a quark line from t_H to t_H
 * Dd_spec is a quark line from t_i to t_f
 * Dd_tf is a quark line from t_f to t_H
 * Ds_spec is a quark line from t_i to t_f
 * Ds_ti is a quark line from t_i to t_H */
template <class FImpl>
template <class mobj, class mobj2, class robj>
void BaryonUtils<FImpl>::Xi_to_Sigma_Q2_Eye_site(const mobj &Dq_loop,
						 const mobj2 &Dd_spec,
						 const mobj &Dd_tf,
						 const mobj2 &Ds_spec,
						 const mobj &Ds_ti,
				                 const Gamma Gamma_H,
				                 const Gamma GammaB_xi,
		                 		 const Gamma GammaB_sigma,
						 robj &result)
{

  Gamma g5(Gamma::Algebra::Gamma5); 

  auto DdG = Dd_spec * GammaB_sigma;
  auto GDs = GammaB_xi * Ds_spec;
  //  Ds * gamma_mu^L * Dq * \gamma_\mu^L * Dd
  auto DsGDqGDd = Ds_ti * Gamma_H * Dq_loop * Gamma_H * Dd_tf;
  //  Ds * gamma_mu^L * Dq * \gamma_\mu^L * Dd * GammaB
  auto DsGDqGDdG = DsGDqGDd * GammaB_sigma;
  //  GammaB * Ds * gamma_mu^L * Dq * \gamma_\mu^L * Dd
  auto GDsGDqGDd = GammaB_xi * DsGDqGDd;
  //  GammaB * Ds * gamma_mu^L * Dq * \gamma_\mu^L * Dd * GammaB
  auto GDsGDqGDdG = GammaB_xi * DsGDqGDdG; 

  for (int ie_s=0; ie_s < 6 ; ie_s++){
    int a_s = epsilon[ie_s][0]; //a
    int b_s = epsilon[ie_s][1]; //b
    int c_s = epsilon[ie_s][2]; //c
    for (int ie_x=0; ie_x < 6 ; ie_x++){
      int a_x = epsilon[ie_x][0]; //a'
      int b_x = epsilon[ie_x][1]; //b'
      int c_x = epsilon[ie_x][2]; //c'
      auto ee = epsilon_sgn[ie_s] * epsilon_sgn[ie_x];
      for (int alpha_x=0; alpha_x<Ns; alpha_x++){
      for (int beta_s=0; beta_s<Ns; beta_s++){
        auto Ds_ab_ab = Ds_spec()(alpha_x, beta_s)(a_x,b_s);
        auto GDs_ab_bb = GDs()(alpha_x, beta_s)(b_x,b_s);
        auto DsGDqGDdG_ab_aa = DsGDqGDdG()(alpha_x, beta_s)(a_x,a_s);
        auto GDsGDqGDdG_ab_ba = GDsGDqGDdG()(alpha_x, beta_s)(a_x,b_s);
        for (int gamma_x=0; gamma_x<Ns; gamma_x++){
        for (int gamma_s=0; gamma_s<Ns; gamma_s++){
          result()(gamma_x,gamma_s)() -= ee * Dd_spec()(gamma_x, gamma_s)(c_x,c_s) * GDsGDqGDdG_ab_ba * Ds_ab_ab;
          result()(gamma_x,gamma_s)() += ee * DdG()(gamma_x, beta_s)(c_x,a_s) * GDsGDqGDd()(alpha_x,gamma_s)(b_x,c_s) * Ds_ab_ab;
          result()(gamma_x,gamma_s)() += ee * Dd_spec()(gamma_x, gamma_s)(c_x,c_s) * DsGDqGDdG_ab_aa * GDs_ab_bb;
          result()(gamma_x,gamma_s)() -= ee * DdG()(gamma_x, beta_s)(c_x,a_s) * DsGDqGDd()(alpha_x,gamma_s)(a_x,c_s) * GDs_ab_bb;
        }}
      }}
    }
  }
}



template<class FImpl>
template <class mobj>
void BaryonUtils<FImpl>::Xi_to_Sigma_Eye(const PropagatorField &qq_loop,
						 const mobj &Dd_spec,
						 const PropagatorField &qd_tf,
						 const mobj &Ds_spec,
						 const PropagatorField &qs_ti,
				                 const Gamma Gamma_H,
				                 const Gamma GammaB_xi,
		                 		 const Gamma GammaB_sigma,
						 const std::string op,
						 SpinMatrixField &stn_corr)
{
  GridBase *grid = qs_ti.Grid();

  auto vcorr= stn_corr.View();
  auto vq_loop = qq_loop.View();
  auto vd_tf = qd_tf.View();
  auto vs_ti = qs_ti.View();

 // accelerator_for(ss, grid->oSites(), grid->Nsimd(), {
  thread_for(ss,grid->oSites(),{
    auto Dq_loop = vq_loop[ss];
    auto Dd_tf = vd_tf[ss];
    auto Ds_ti = vs_ti[ss];
    sobj result=Zero();
    if(op == "Q1"){
      Xi_to_Sigma_Q1_Eye_site(Dq_loop,Dd_spec,Dd_tf,Ds_spec,Ds_ti,Gamma_H,GammaB_xi,GammaB_sigma,result);
    } else if(op == "Q2"){
      Xi_to_Sigma_Q2_Eye_site(Dq_loop,Dd_spec,Dd_tf,Ds_spec,Ds_ti,Gamma_H,GammaB_xi,GammaB_sigma,result);
    } else {
      assert(0 && "Weak Operator not correctly specified");
    }
      vcorr[ss] = result; 
  }  );//end loop over lattice sites
}


NAMESPACE_END(Grid);
