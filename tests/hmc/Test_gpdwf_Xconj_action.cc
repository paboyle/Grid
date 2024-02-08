    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/hmc/Test_gpdwf_Xconj_action.cc

    Copyright (C) 2015

Author: Christopher Kelly <ckelly@bnl.gov>
Author: Dennis Bollweg <db3516@columbia.edu>
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

/**
 * Tests for the action of X-conjugate Dirac operator
 *
 */

#include <Grid/Grid.h>

using namespace std;
using namespace Grid;

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  const int Ls=4;

  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  auto const & latt_size = GridDefaultLatt();
  size_t vol4d = latt_size[0]*latt_size[1]*latt_size[2]*latt_size[3];

  std::vector<int> seeds4({1,2,3,4});
  std::vector<int> seeds5({5,6,7,8});
  GridParallelRNG          RNG5(FGrid);  RNG5.SeedFixedIntegers(seeds5);
  GridParallelRNG          RNG4(UGrid);  RNG4.SeedFixedIntegers(seeds4);
  GridParallelRNG          RNG5rb(FrbGrid);  RNG5.SeedFixedIntegers(seeds5);

  LatticeGaugeField Umu(UGrid); 
  SU<Nc>::HotConfiguration(RNG4, Umu);

  Gamma C = Gamma(Gamma::Algebra::MinusGammaY) * Gamma(Gamma::Algebra::GammaT);
  Gamma g5 = Gamma(Gamma::Algebra::Gamma5);
  Gamma X = C*g5;
  
  //Set up a regular MDWF action instance as well as X-conj and Xbar-conj versions
  RealD m1=0.1;
  RealD m2=1.0;

  RealD M5=1.8;
  RealD mob_b=1.5;

  typedef typename GparityMobiusFermionD::FermionField FermionField2f;
  typedef typename XconjugateMobiusFermionD::FermionField FermionField1f;

  GparityMobiusFermionD ::ImplParams params;
  std::vector<int> twists({1,1,0,1}); //G-parity in x,y  periodic in z, antiperiodic in t
  params.twists = twists;    
  GparityMobiusFermionD reg_action_m1(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,m1,M5,mob_b,mob_b-1.,params);
  GparityMobiusFermionD reg_action_m2(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,m2,M5,mob_b,mob_b-1.,params);

  XconjugateMobiusFermionD::ImplParams xparams;
  xparams.twists = twists;
  xparams.boundary_phase = 1.0;
  XconjugateMobiusFermionD xconj_action_m1(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,m1,M5,mob_b,mob_b-1.,xparams);
  XconjugateMobiusFermionD xconj_action_m2(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,m2,M5,mob_b,mob_b-1.,xparams);

  xparams.boundary_phase = -1.0;
  XconjugateMobiusFermionD xbarconj_action_m1(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,m1,M5,mob_b,mob_b-1.,xparams);
  XconjugateMobiusFermionD xbarconj_action_m2(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,m2,M5,mob_b,mob_b-1.,xparams);

  //Gauge BCs
  typedef ConjugateGimplD GimplD;
  std::vector<int> gauge_twists = twists;
  gauge_twists[3] = 0; //periodic in time
  GimplD::setDirections(twists);

  typedef typename GparityMobiusFermionD::Impl_t FImplD;
  typedef typename XconjugateMobiusFermionD::Impl_t XFImplD;

  //First demonstrate that the X-conjugate action is computed correctly
  //It should be the same as the G-parity Dirac operator applied to an X-conjugate vector
  {   
     typedef TwoFlavourEvenOddRatioPseudoFermionAction<XFImplD> XconjTwoFlavorHMC;
     typedef TwoFlavourEvenOddRatioPseudoFermionAction<FImplD> GparityTwoFlavorHMC;
     
     ConjugateGradient<FermionField2f>  CG2f(1e-10,10000);
     ConjugateGradient<FermionField1f>  CG1f(1e-10,10000);

     XconjTwoFlavorHMC hmc_xconj(xconj_action_m2, xconj_action_m1, CG1f, CG1f);
     GparityTwoFlavorHMC hmc_gp(reg_action_m2, reg_action_m1, CG2f, CG2f);

     RealD scale = std::sqrt(0.5);
     FermionField1f eta1f(FGrid);
     gaussian(RNG5,eta1f);
     
     FermionField2f eta2f(FGrid);

     FermionField1f tmp = -(X*conjugate(eta1f));
     PokeIndex<GparityFlavourIndex>(eta2f, eta1f, 0);
     PokeIndex<GparityFlavourIndex>(eta2f, tmp, 1);
     
     hmc_xconj.setPhi(eta1f);
     RealD S_xconj = 2. * hmc_xconj.S(Umu);

     hmc_gp.setPhi(eta2f);
     RealD S_gp = hmc_gp.S(Umu);
     
     std::cout << "Test  2 rho^dag [M_X M_X^dag]^-1 rho  = \\hat rho^dag [M_GP M_GP^dag] \\hat rho (expect 0): " << S_gp - S_xconj << std::endl;
     assert(fabs(S_gp - S_xconj) < 1e-8);
  }
  //Next demonstrate the G-parity Dirac operator acting on a general 2f vector can be computed in terms of the X-conjugate action
  {   
     typedef TwoFlavourEvenOddRatioPseudoFermionAction<XFImplD> XconjTwoFlavorHMC;
     typedef TwoFlavourEvenOddRatioPseudoFermionAction<FImplD> GparityTwoFlavorHMC;
     
     ConjugateGradient<FermionField2f>  CG2f(1e-10,10000);
     ConjugateGradient<FermionField1f>  CG1f(1e-10,10000);

     XconjTwoFlavorHMC hmc_xconj(xconj_action_m2, xconj_action_m1, CG1f, CG1f);
     XconjTwoFlavorHMC hmc_xbarconj(xbarconj_action_m2, xbarconj_action_m1, CG1f, CG1f);
     GparityTwoFlavorHMC hmc_gp(reg_action_m2, reg_action_m1, CG2f, CG2f);

     RealD scale = std::sqrt(0.5);
     FermionField2f eta2f(FGrid);
     gaussian(RNG5,eta2f);

     FermionField1f eta2f_0 = PeekIndex<GparityFlavourIndex>(eta2f,0);
     FermionField1f eta2f_1 = PeekIndex<GparityFlavourIndex>(eta2f,1);
     
     FermionField1f rho1f = 0.5*( eta2f_0 + X*conjugate(eta2f_1) );
     FermionField1f tau1f = Complex(0,0.5)*( eta2f_0 - X*conjugate(eta2f_1) );
     
     hmc_gp.setPhi(eta2f);
     RealD S_gp = hmc_gp.S(Umu);

     hmc_xconj.setPhi(rho1f);   
     RealD S_xconj = 2.*hmc_xconj.S(Umu);
     hmc_xconj.setPhi(tau1f);
     RealD S_xconj2 = 2.*hmc_xconj.S(Umu);

     RealD S_xtotal = S_xconj + S_xconj2;

     std::cout << "Test eta^dag [M_GP M_GP^dag]^-1 eta  = 2 rho^dag [M_X M_X^dag]^-1 rho + 2 tau^dag [M_X M_X^dag]^-1 tau (expect 0): " << S_xtotal-S_gp  << std::endl;
     assert(fabs(S_xtotal-S_gp) < 1e-8);

     //The above can also be computed with the Xbar-conjugate action after transforming the source
     FermionField1f xi1f = Complex(0,-1)*tau1f;
     hmc_xbarconj.setPhi(xi1f);
     RealD S_xbarconj = 2.*hmc_xbarconj.S(Umu);

     S_xtotal = S_xconj + S_xbarconj;

     std::cout << "Test eta^dag [M_GP M_GP^dag]^-1 eta  = 2 rho^dag [M_X M_X^dag]^-1 rho + 2 xi^dag [M_Xbar M_Xbar^dag]^-1 xi (expect 0): " << S_xtotal-S_gp  << std::endl;
     assert(fabs(S_xtotal-S_gp) < 1e-8);     
  }
  //Do the same but with a rational approximation to (M M^dag)^-1/2
  {   
    typedef GeneralEvenOddRatioRationalPseudoFermionAction<FImplD> GparityRHMC;
    typedef GeneralEvenOddRatioRationalPseudoFermionAction<XFImplD> XconjRHMC;

    RationalActionParams rat_act_params;
    rat_act_params.inv_pow  = 2; // (M^dag M)^{1/2}
    rat_act_params.precision= 60;
    rat_act_params.MaxIter  = 50000;
    rat_act_params.action_degree = 15;
    rat_act_params.action_tolerance = 1.0e-10;
    rat_act_params.md_degree = 12;
    rat_act_params.md_tolerance = 1.0e-5;
    rat_act_params.lo = 0.01;
    rat_act_params.hi = 80.0;
    
    GparityRHMC rhmc_gp(reg_action_m2, reg_action_m1, rat_act_params);
    XconjRHMC rhmc_xconj(xconj_action_m2, xconj_action_m1, rat_act_params);
    XconjRHMC rhmc_xbarconj(xbarconj_action_m2, xbarconj_action_m1, rat_act_params);

    RealD scale = std::sqrt(0.5);
    FermionField2f eta2f(FGrid);
    gaussian(RNG5,eta2f);
    
    FermionField1f eta2f_0 = PeekIndex<GparityFlavourIndex>(eta2f,0);
    FermionField1f eta2f_1 = PeekIndex<GparityFlavourIndex>(eta2f,1);
    
    FermionField1f rho1f = 0.5*( eta2f_0 + X*conjugate(eta2f_1) );
    FermionField1f tau1f = Complex(0,0.5)*( eta2f_0 - X*conjugate(eta2f_1) );
    
    rhmc_gp.setPhi(eta2f);
    RealD S_gp = rhmc_gp.S(Umu);

    rhmc_xconj.setPhi(rho1f);   
    RealD S_xconj = 2.*rhmc_xconj.S(Umu);
    rhmc_xconj.setPhi(tau1f);
    RealD S_xconj2 = 2.*rhmc_xconj.S(Umu);

    RealD S_xtotal = S_xconj + S_xconj2;
    
    std::cout << "Test eta^dag [M_GP M_GP^dag]^-1/2 eta  = 2 rho^dag [M_X M_X^dag]^-1/2 rho + 2 tau^dag [M_X M_X^dag]^-1/2 tau using RHMC (expect 0): " << S_xtotal-S_gp  << std::endl;
    assert(fabs(S_xtotal-S_gp) < 1e-8);

    //The above can also be computed with the Xbar-conjugate action after transforming the source
    FermionField1f xi1f = Complex(0,-1)*tau1f;
    rhmc_xbarconj.setPhi(xi1f);
    RealD S_xbarconj = 2.*rhmc_xbarconj.S(Umu);
    
    S_xtotal = S_xconj + S_xbarconj;
    
    std::cout << "Test eta^dag [M_GP M_GP^dag]^-1/2 eta  = 2 rho^dag [M_X M_X^dag]^-1/2 rho + 2 xi^dag [M_Xbar M_Xbar^dag]^-1/2 xi using RHMC (expect 0): " << S_xtotal-S_gp  << std::endl;
    assert(fabs(S_xtotal-S_gp) < 1e-8);     
  }

  //Examine the EOFA Dirac operator
  {   
    typedef XconjugateMobiusEOFAFermionD XconjEOFAaction;
    typedef ExactOneFlavourRatioPseudoFermionAction<XFImplD> XconjEOFA;
    typedef GparityMobiusEOFAFermionD GparityEOFAaction;
    typedef ExactOneFlavourRatioPseudoFermionAction<FImplD> GparityEOFA;

    ConjugateGradient<FermionField2f>  CG2f(1e-10,10000);
    ConjugateGradient<FermionField1f>  CG1f(1e-10,10000);
     
    RealD iml = m1;
    RealD ipv = m2;

    OneFlavourRationalParams eofa_params; //not needed because we are not doing the refresh
    eofa_params.lo = 1.0;
    eofa_params.hi = 2.0;
    eofa_params.degree = 12;
    eofa_params.tolerance = 1e-4; //this is just for doing the Remez

    xparams.boundary_phase = 1.0;
    XconjEOFAaction Lop_xconj(Umu, *FGrid, *FrbGrid, *UGrid, *UrbGrid, iml, iml, ipv, 0.0, -1, M5, mob_b, mob_b-1.,xparams);
    XconjEOFAaction Rop_xconj(Umu, *FGrid, *FrbGrid, *UGrid, *UrbGrid, ipv, iml, ipv, -1.0, 1, M5, mob_b, mob_b-1.,xparams);
    XconjEOFA eofa_xconj(Lop_xconj, Rop_xconj, CG1f, eofa_params);

    xparams.boundary_phase = -1.0;
    XconjEOFAaction Lop_xbarconj(Umu, *FGrid, *FrbGrid, *UGrid, *UrbGrid, iml, iml, ipv, 0.0, -1, M5, mob_b, mob_b-1.,xparams);
    XconjEOFAaction Rop_xbarconj(Umu, *FGrid, *FrbGrid, *UGrid, *UrbGrid, ipv, iml, ipv, -1.0, 1, M5, mob_b, mob_b-1.,xparams);
    XconjEOFA eofa_xbarconj(Lop_xbarconj, Rop_xbarconj, CG1f, eofa_params);

    GparityEOFAaction Lop_gp(Umu, *FGrid, *FrbGrid, *UGrid, *UrbGrid, iml, iml, ipv, 0.0, -1, M5, mob_b, mob_b-1.,params);
    GparityEOFAaction Rop_gp(Umu, *FGrid, *FrbGrid, *UGrid, *UrbGrid, ipv, iml, ipv, -1.0, 1, M5, mob_b, mob_b-1.,params);
    GparityEOFA eofa_gp(Lop_gp, Rop_gp, CG2f, eofa_params);

    RealD scale = std::sqrt(0.5);
    FermionField2f eta2f(FGrid);
    gaussian(RNG5,eta2f);
    
    FermionField1f eta2f_0 = PeekIndex<GparityFlavourIndex>(eta2f,0);
    FermionField1f eta2f_1 = PeekIndex<GparityFlavourIndex>(eta2f,1);
    
    FermionField1f rho1f = 0.5*( eta2f_0 + X*conjugate(eta2f_1) );
    FermionField1f tau1f = Complex(0,0.5)*( eta2f_0 - X*conjugate(eta2f_1) );
    
    eofa_gp.setPhi(eta2f);
    RealD S_gp = eofa_gp.S(Umu);

    eofa_xconj.setPhi(rho1f);   
    RealD S_xconj = 2.*eofa_xconj.S(Umu);
    eofa_xconj.setPhi(tau1f);
    RealD S_xconj2 = 2.*eofa_xconj.S(Umu);

    RealD S_xtotal = S_xconj + S_xconj2;
    
    std::cout << "Test eta^dag [M_GP M_GP^dag]^-1/2 eta  = 2 rho^dag [M_X M_X^dag]^-1/2 rho + 2 tau^dag [M_X M_X^dag]^-1/2 tau using EOFA (expect 0): " << S_xtotal-S_gp  << std::endl;
    assert(fabs(S_xtotal-S_gp) < 1e-8);

    //The above can also be computed with the Xbar-conjugate action after transforming the source
    FermionField1f xi1f = Complex(0,-1)*tau1f;
    eofa_xbarconj.setPhi(xi1f);
    RealD S_xbarconj = 2.*eofa_xbarconj.S(Umu);
    
    S_xtotal = S_xconj + S_xbarconj;
    
    std::cout << "Test eta^dag [M_GP M_GP^dag]^-1/2 eta  = 2 rho^dag [M_X M_X^dag]^-1/2 rho + 2 xi^dag [M_Xbar M_Xbar^dag]^-1/2 xi using EOFA (expect 0): " << S_xtotal-S_gp  << std::endl;
    assert(fabs(S_xtotal-S_gp) < 1e-8);     
  }




  Grid_finalize();
}
