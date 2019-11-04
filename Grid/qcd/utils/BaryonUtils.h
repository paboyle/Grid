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
  
  static constexpr int epsilon[6][3] = {{0,1,2},{1,2,0},{2,0,1},{0,2,1},{2,1,0},{1,0,2}};
  static constexpr int epsilon_sgn[6]= {1,1,1,-1,-1,-1};

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
};

template <class FImpl>
constexpr int BaryonUtils<FImpl>::epsilon[6][3];
template <class FImpl>
constexpr int BaryonUtils<FImpl>::epsilon_sgn[6];

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
	    result()()() += static_cast<Complex>(epsilon_sgn[ie_left] * epsilon_sgn[ie_right]) * pD1()(gamma_left,gamma_left)(c_right,c_left)*D2g()(alpha_right,beta_left)(a_right,a_left)*gD3()(alpha_right,beta_left)(b_right,b_left);
          }}}
  	}	  
        //This is the \delta_{456}^{231} part
	if (wick_contraction[1]){
          auto pD1g = pD1 * GammaB_left;
	  for (int alpha_right=0; alpha_right<Ns; alpha_right++){
	  for (int beta_left=0; beta_left<Ns; beta_left++){
	  for (int gamma_left=0; gamma_left<Ns; gamma_left++){
	    result()()() += static_cast<Complex>(epsilon_sgn[ie_left] * epsilon_sgn[ie_right]) * pD1g()(gamma_left,beta_left)(c_right,a_left)*D2()(alpha_right,beta_left)(a_right,b_left)*gD3()(alpha_right,gamma_left)(b_right,c_left);
          }}}
        }	  
        //This is the \delta_{456}^{312} part
	if (wick_contraction[2]){
          auto gD3g = gD3 * GammaB_left;
	  for (int alpha_right=0; alpha_right<Ns; alpha_right++){
	  for (int beta_left=0; beta_left<Ns; beta_left++){
	  for (int gamma_left=0; gamma_left<Ns; gamma_left++){
	    result()()() += static_cast<Complex>(epsilon_sgn[ie_left] * epsilon_sgn[ie_right]) * pD1()(gamma_left,beta_left)(c_right,b_left)*D2()(alpha_right,gamma_left)(a_right,c_left)*gD3g()(alpha_right,beta_left)(b_right,a_left);
          }}}
        }	  
        //This is the \delta_{456}^{132} part
	if (wick_contraction[3]){
          auto gD3g = gD3 * GammaB_left;
	  for (int alpha_right=0; alpha_right<Ns; alpha_right++){
	  for (int beta_left=0; beta_left<Ns; beta_left++){
	  for (int gamma_left=0; gamma_left<Ns; gamma_left++){
	    result()()() -= static_cast<Complex>(epsilon_sgn[ie_left] * epsilon_sgn[ie_right]) * pD1()(gamma_left,gamma_left)(c_right,c_left)*D2()(alpha_right,beta_left)(a_right,b_left)*gD3g()(alpha_right,beta_left)(b_right,a_left);
          }}}
        }	  
        //This is the \delta_{456}^{321} part
	if (wick_contraction[4]){
          auto D2g = D2 * GammaB_left;
	  for (int alpha_right=0; alpha_right<Ns; alpha_right++){
	  for (int beta_left=0; beta_left<Ns; beta_left++){
	  for (int gamma_left=0; gamma_left<Ns; gamma_left++){
	    result()()() -= static_cast<Complex>(epsilon_sgn[ie_left] * epsilon_sgn[ie_right]) * pD1()(gamma_left,beta_left)(c_right,b_left)*D2g()(alpha_right,beta_left)(a_right,a_left)*gD3()(alpha_right,gamma_left)(b_right,c_left);
          }}}
        }	  
        //This is the \delta_{456}^{213} part
	if (wick_contraction[5]){
          auto pD1g = pD1 * GammaB_left;
	  for (int alpha_right=0; alpha_right<Ns; alpha_right++){
	  for (int beta_left=0; beta_left<Ns; beta_left++){
	  for (int gamma_left=0; gamma_left<Ns; gamma_left++){
	    result()()() -= static_cast<Complex>(epsilon_sgn[ie_left] * epsilon_sgn[ie_right]) * pD1g()(gamma_left,beta_left)(c_right,a_left)*D2()(alpha_right,gamma_left)(a_right,c_left)*gD3()(alpha_right,beta_left)(b_right,b_left);
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
NAMESPACE_END(Grid);
