/*************************************************************************************
 
 Grid physics library, www.github.com/paboyle/Grid
 
 Source file: ./lib/qcd/utils/BaryonUtils.h
 
 Copyright (C) 2019
 
 Author: Felix Erben <felix.erben@ed.ac.uk>
 Author: Michael Marshall <Michael.Marshall@ed.ac.uk>

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

#undef DELTA_F_EQ_2

template <typename FImpl>
class BaryonUtils 
{
public:
  typedef typename FImpl::ComplexField ComplexField;
  typedef typename FImpl::FermionField FermionField;
  typedef typename FImpl::PropagatorField PropagatorField;

  typedef typename FImpl::SitePropagator pobj;
  typedef typename FImpl::SiteSpinor vobj;
  typedef typename vobj::scalar_object sobj;
  typedef typename vobj::scalar_type scalar_type;
  typedef typename vobj::vector_type vector_type;

  static void ContractBaryons(const PropagatorField &q1_src,
					 const PropagatorField &q2_src,
					 const PropagatorField &q3_src,
					 const Gamma GammaA,
					 const Gamma GammaB,
					 const char quarks_snk[],
					 const char quarks_src[],
					 const int parity,
					 ComplexField &baryon_corr);
 
};


template<class FImpl>
void BaryonUtils<FImpl>::ContractBaryons(const PropagatorField &q1_src,
						 const PropagatorField &q2_src,
						 const PropagatorField &q3_src,
						 const Gamma GammaA,
						 const Gamma GammaB,
						 const char quarks_snk[],
						 const char quarks_src[],
						 const int parity,
						 ComplexField &baryon_corr)
{

  assert(parity==1 || parity == -1 && "Parity must be +1 or -1");

  GridBase *grid = q1_src.Grid();

  Gamma g4(Gamma::Algebra::GammaT); //needed for parity P_\pm = 0.5*(1 \pm \gamma_4)

  std::vector<std::vector<int>> epsilon = {{0,1,2},{1,2,0},{2,0,1},{0,2,1},{2,1,0},{1,0,2}};
  std::vector<int> epsilon_sgn = {1,1,1,-1,-1,-1};
  std::vector<int> wick_contraction = {0,0,0,0,0,0};

  for (int ie=0; ie < 6 ; ie++)
    if (quarks_src[0] == quarks_snk[epsilon[ie][0]] && quarks_src[1] == quarks_snk[epsilon[ie][1]] && quarks_src[2] == quarks_snk[epsilon[ie][2]])
      wick_contraction[ie]=1;

  typedef typename ComplexField::vector_object vobj;
  LatticeView<vobj> vbaryon_corr{ baryon_corr };
  accelerator_for(ss, grid->oSites(), grid->Nsimd(), {

    LatticeView<pobj> v1(q1_src);
    LatticeView<pobj> v2(q2_src);
    LatticeView<pobj> v3(q3_src);
    auto D1 = v1[ss];
    auto D2 = v2[ss];
    auto D3 = v3[ss];

    auto gD1a = GammaA * GammaA * D1;
    auto gD1b = GammaA * g4 * GammaA * D1;
    auto pD1 = 0.5* (gD1a + (double)parity * gD1b);
    auto gD3 = GammaB * D3;

    vobj result{ 0 };
    
    for (int ie_src=0; ie_src < 6 ; ie_src++){
      int a_src = epsilon[ie_src][0]; //a
      int b_src = epsilon[ie_src][1]; //b
      int c_src = epsilon[ie_src][2]; //c
      for (int ie_snk=0; ie_snk < 6 ; ie_snk++){
        int a_snk = epsilon[ie_snk][0]; //a'
        int b_snk = epsilon[ie_snk][1]; //b'
        int c_snk = epsilon[ie_snk][2]; //c'
        //This is the \delta_{456}^{123} part
	if (wick_contraction[0]){
          auto D2g = D2 * GammaB;
	  for (int alpha_snk=0; alpha_snk<Ns; alpha_snk++){
	  for (int beta_src=0; beta_src<Ns; beta_src++){
	  for (int gamma_src=0; gamma_src<Ns; gamma_src++){
	    result()()() += epsilon_sgn[ie_src] * epsilon_sgn[ie_snk] * pD1()(gamma_src,gamma_src)(c_snk,c_src)*D2g()(alpha_snk,beta_src)(a_snk,a_src)*gD3()(alpha_snk,beta_src)(b_snk,b_src);
          }}}
  	}	  
        //This is the \delta_{456}^{231} part
	if (wick_contraction[1]){
          auto pD1g = pD1 * GammaB;
	  for (int alpha_snk=0; alpha_snk<Ns; alpha_snk++){
	  for (int beta_src=0; beta_src<Ns; beta_src++){
	  for (int gamma_src=0; gamma_src<Ns; gamma_src++){
	    result()()() += epsilon_sgn[ie_src] * epsilon_sgn[ie_snk] * pD1g()(gamma_src,beta_src)(c_snk,a_src)*D2()(alpha_snk,beta_src)(a_snk,b_src)*gD3()(alpha_snk,gamma_src)(b_snk,c_src);
          }}}
        }	  
        //This is the \delta_{456}^{312} part
	if (wick_contraction[2]){
          auto gD3g = gD3 * GammaB;
	  for (int alpha_snk=0; alpha_snk<Ns; alpha_snk++){
	  for (int beta_src=0; beta_src<Ns; beta_src++){
	  for (int gamma_src=0; gamma_src<Ns; gamma_src++){
	    result()()() += epsilon_sgn[ie_src] * epsilon_sgn[ie_snk] * pD1()(gamma_src,beta_src)(c_snk,b_src)*D2()(alpha_snk,gamma_src)(a_snk,c_src)*gD3g()(alpha_snk,beta_src)(b_snk,a_src);
          }}}
        }	  
        //This is the \delta_{456}^{132} part
	if (wick_contraction[3]){
          auto gD3g = gD3 * GammaB;
	  for (int alpha_snk=0; alpha_snk<Ns; alpha_snk++){
	  for (int beta_src=0; beta_src<Ns; beta_src++){
	  for (int gamma_src=0; gamma_src<Ns; gamma_src++){
	    result()()() -= epsilon_sgn[ie_src] * epsilon_sgn[ie_snk] * pD1()(gamma_src,gamma_src)(c_snk,c_src)*D2()(alpha_snk,beta_src)(a_snk,b_src)*gD3g()(alpha_snk,beta_src)(b_snk,a_src);
          }}}
        }	  
        //This is the \delta_{456}^{321} part
	if (wick_contraction[4]){
          auto D2g = D2 * GammaB;
	  for (int alpha_snk=0; alpha_snk<Ns; alpha_snk++){
	  for (int beta_src=0; beta_src<Ns; beta_src++){
	  for (int gamma_src=0; gamma_src<Ns; gamma_src++){
	    result()()() -= epsilon_sgn[ie_src] * epsilon_sgn[ie_snk] * pD1()(gamma_src,beta_src)(c_snk,b_src)*D2g()(alpha_snk,beta_src)(a_snk,a_src)*gD3()(alpha_snk,gamma_src)(b_snk,c_src);
          }}}
        }	  
        //This is the \delta_{456}^{213} part
	if (wick_contraction[5]){
          auto pD1g = pD1 * GammaB;
	  for (int alpha_snk=0; alpha_snk<Ns; alpha_snk++){
	  for (int beta_src=0; beta_src<Ns; beta_src++){
	  for (int gamma_src=0; gamma_src<Ns; gamma_src++){
	    result()()() -= epsilon_sgn[ie_src] * epsilon_sgn[ie_snk] * pD1g()(gamma_src,beta_src)(c_snk,a_src)*D2()(alpha_snk,gamma_src)(a_snk,c_src)*gD3()(alpha_snk,beta_src)(b_snk,b_src);
          }}}
        }	  
      }
    }
    vbaryon_corr[ss] = result; 

  } );//end loop over lattice sites
}
NAMESPACE_END(Grid);
