/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/smearing/HISQSmearing.h

Copyright (C) 2023

Author: D. A. Clarke <clarke.davida@gmail.com> 

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

See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
/*
    @file HISQSmearing.h
    @brief Declares classes related to the HISQ action 
*/


#pragma once
#include <Grid/Grid.h>
#include <Grid/lattice/PaddedCell.h>
#include <Grid/stencil/GeneralLocalStencil.h>
NAMESPACE_BEGIN(Grid);


/*!  @brief structure holding the link treatment for a given smear */
template<class floatT>
struct SmearingParameters{
    floatT c_1;               // 1 link
    floatT c_naik;            // Naik term
    floatT c_3;               // 3 link
    floatT c_5;               // 5 link
    floatT c_7;               // 7 link
    floatT c_lp;              // Lepage
    SmearingParameters(floatT c1, floatT cnaik, floatT c3, floatT c5, floatT c7, floatT clp) 
        : c_1(c1),
          c_naik(cnaik),
          c_3(c3),
          c_5(c5),
          c_7(c7),
          c_lp(clp){}
};


// There are 6 quarks in nature, and 3 never need a Naik epsilon
int const GRID_MAX_NAIK = 3;


/*!  @brief structure holding all input parameters related to the HISQ action */
template<class floatT>
struct HISQParameters{
    // Structure from QOP/QDP 
    int n_naiks;
    std::array<floatT,GRID_MAX_NAIK> eps_naiks;
    floatT fat7_c1  ; floatT fat7_c3  ; floatT fat7_c5  ; floatT fat7_c7  ; floatT fat7_clp;
    floatT asqtad_c1; floatT asqtad_c3; floatT asqtad_c5; floatT asqtad_c7; floatT asqtad_clp; floatT asqtad_cnaik;
    floatT diff_c1  ; floatT diff_cnaik;
    HISQParameters(int n_naiks_in, std::array<floatT,GRID_MAX_NAIK> eps_naiks_in, 
                   floatT fat7_one_link  , floatT fat7_three_staple  , floatT fat7_five_staple  , floatT fat7_seven_staple  , floatT fat7_lepage, 
                   floatT asqtad_one_link, floatT asqtad_three_staple, floatT asqtad_five_staple, floatT asqtad_seven_staple, floatT asqtad_lepage,
                   floatT asqtad_naik    , floatT difference_one_link, floatT difference_naik) 
        : n_naiks(n_naiks_in),
          eps_naiks(eps_naiks_in), 
          fat7_c1(fat7_one_link),
          fat7_c3(fat7_three_staple),
          fat7_c5(fat7_five_staple),
          fat7_c7(fat7_seven_staple),
          fat7_clp(fat7_lepage),
          asqtad_c1(asqtad_one_link),
          asqtad_c3(asqtad_three_staple),
          asqtad_c5(asqtad_five_staple),
          asqtad_c7(asqtad_seven_staple),
          asqtad_clp(asqtad_lepage),
          asqtad_cnaik(asqtad_naik),
          diff_c1(difference_one_link),
          diff_cnaik(difference_naik) {}
};


/*!  @brief Sometimes in the U(3) projection we use SVD cuts; here we collect related parameters */
template<class floatT>
struct HISQReunitSVDParameters{
    // Structure from QOP/QDP 
    bool   allow_svd;
    bool   svd_only;
    floatT svd_rel_error;
    floatT svd_abs_error;
    floatT force_filter;
    HISQReunitSVDParameters(bool allow, bool only, floatT rel, floatT abs, floatT filter)
        : allow_svd(allow),
          svd_only(only),
          svd_rel_error(rel),
          svd_abs_error(abs),
          force_filter(filter) {}
};


/*!  @brief Get the link U_mu(x).  */
template<class link> accelerator_inline
auto getLink(const link& __restrict__ U, GeneralStencilEntry* x, int mu) {
    return coalescedReadGeneralPermute(U[x->_offset](mu), x->_permute, Nd); 
}
#define setLink coalescedWrite 


/*! @brief Figure out the stencil index from mu and nu. */
accelerator_inline 
int HISQStencilIndex(int mu, int nu, int rho=0, int sig=0, std::string kind="3STAPLE") {
    int res;
    if (kind=="3STAPLE") 
        res = 5*(nu + Nd*mu);
    else if (kind=="5STAPLE")
        res = 17*(rho + Nd*nu + Nd*Nd*mu);
    else if (kind=="7STAPLE")
        // seems correct
        res = 46*(sig + Nd*rho + Nd*Nd*nu + Nd*Nd*Nd*mu);
    else 
        Grid_error("Unknown staple kind",kind);
    return res;
}


/*! @brief Create various stencils needed for HISQ calculations */ 
inline
std::vector<Coordinate> createHISQStencil(std::string kind="3STAPLE") {
    std::vector<Coordinate> shifts;
    // We allow nu=mu and rho=nu, rho=mu to make indexing easier, but these
    // entries will not be used.
    if (kind=="3STAPLE") {
        for(int mu=0;mu<Nd;mu++)
        for(int nu=0;nu<Nd;nu++) {
            appendShift<Nd>(shifts,mu);
            appendShift<Nd>(shifts,nu);
            appendShift<Nd>(shifts,shiftSignal::NO_SHIFT);
            appendShift<Nd>(shifts,mu,Back(nu));
            appendShift<Nd>(shifts,Back(nu));
        }
    } else if (kind=="5STAPLE") {
        for(int mu =0;mu <Nd;mu++)
        for(int nu =0;nu <Nd;nu++) 
        for(int rho=0;rho<Nd;rho++) {
            appendShift<Nd>(shifts,nu,Back(rho));
            appendShift<Nd>(shifts,nu);
            appendShift<Nd>(shifts,Back(rho));
            appendShift<Nd>(shifts,shiftSignal::NO_SHIFT);
            appendShift<Nd>(shifts,rho);
            appendShift<Nd>(shifts,Back(nu),Back(rho));
            appendShift<Nd>(shifts,Back(nu));
            appendShift<Nd>(shifts,Back(nu),rho);
            appendShift<Nd>(shifts,mu,nu,Back(rho));
            appendShift<Nd>(shifts,mu,nu);
            appendShift<Nd>(shifts,mu,Back(rho));
            appendShift<Nd>(shifts,mu);
            appendShift<Nd>(shifts,mu,rho);
            appendShift<Nd>(shifts,mu,Back(nu),Back(rho));
            appendShift<Nd>(shifts,mu,Back(nu));
            appendShift<Nd>(shifts,mu,Back(nu),rho);
            appendShift<Nd>(shifts,nu,rho);
        }
    } else if (kind=="7STAPLE") {
        for(int mu =0;mu <Nd;mu++)
        for(int nu =0;nu <Nd;nu++) 
        for(int rho=0;rho<Nd;rho++)
        for(int sig=0;sig<Nd;sig++) {
            appendShift<Nd>(shifts,shiftSignal::NO_SHIFT);
            appendShift<Nd>(shifts,mu);
            appendShift<Nd>(shifts,mu,nu);
            appendShift<Nd>(shifts,mu,nu,rho);
            appendShift<Nd>(shifts,mu,nu,rho,Back(sig));
            appendShift<Nd>(shifts,mu,nu,Back(rho));
            appendShift<Nd>(shifts,mu,nu,Back(rho),Back(sig));
            appendShift<Nd>(shifts,mu,nu,Back(sig));
            appendShift<Nd>(shifts,mu,Back(nu));
            appendShift<Nd>(shifts,mu,Back(nu),rho);
            appendShift<Nd>(shifts,mu,Back(nu),rho,sig);
            appendShift<Nd>(shifts,mu,Back(nu),rho,Back(sig));
            appendShift<Nd>(shifts,mu,Back(nu),Back(rho));
            appendShift<Nd>(shifts,mu,Back(nu),Back(rho),sig);
            appendShift<Nd>(shifts,mu,Back(nu),Back(rho),Back(sig));
            appendShift<Nd>(shifts,mu,Back(nu),sig);
            appendShift<Nd>(shifts,mu,Back(nu),Back(sig));
            appendShift<Nd>(shifts,mu,rho);
            appendShift<Nd>(shifts,mu,rho,sig);
            appendShift<Nd>(shifts,mu,rho,Back(sig));
            appendShift<Nd>(shifts,mu,Back(rho));
            appendShift<Nd>(shifts,mu,Back(rho),sig);
            appendShift<Nd>(shifts,mu,Back(rho),Back(sig));
            appendShift<Nd>(shifts,mu,sig);
            appendShift<Nd>(shifts,mu,Back(sig));
            appendShift<Nd>(shifts,nu);
            appendShift<Nd>(shifts,nu,rho);
            appendShift<Nd>(shifts,nu,rho,sig);
            appendShift<Nd>(shifts,nu,rho,Back(sig));
            appendShift<Nd>(shifts,nu,Back(rho));
            appendShift<Nd>(shifts,nu,Back(rho),sig);
            appendShift<Nd>(shifts,nu,Back(rho),Back(sig));
            appendShift<Nd>(shifts,rho);
            appendShift<Nd>(shifts,rho,Back(nu));
            appendShift<Nd>(shifts,rho,sig);
            appendShift<Nd>(shifts,rho,Back(sig));
            appendShift<Nd>(shifts,Back(nu));
            appendShift<Nd>(shifts,Back(nu),rho);
            appendShift<Nd>(shifts,Back(nu),rho,sig);
            appendShift<Nd>(shifts,Back(nu),rho,Back(sig));
            appendShift<Nd>(shifts,Back(nu),Back(rho));
            appendShift<Nd>(shifts,Back(nu),Back(rho),sig);
            appendShift<Nd>(shifts,Back(nu),Back(rho),Back(sig));
            appendShift<Nd>(shifts,Back(rho));
            appendShift<Nd>(shifts,Back(rho),sig);
            appendShift<Nd>(shifts,Back(rho),Back(sig));
        }
    } else { 
        Grid_error("Unknown staple kind",kind);
    }
    return shifts;
}


/*! @brief Retrieve 3-link stencil entries. */
template<class acc> accelerator_inline
std::tuple<GeneralStencilEntry*,GeneralStencilEntry*,GeneralStencilEntry*,GeneralStencilEntry*,GeneralStencilEntry*> 
get3StaplePoints(acc sView, int sIndex, int site) {
    GeneralStencilEntry* x_p_mu      = sView.GetEntry(sIndex+0,site);
    GeneralStencilEntry* x_p_nu      = sView.GetEntry(sIndex+1,site);
    GeneralStencilEntry* x           = sView.GetEntry(sIndex+2,site);
    GeneralStencilEntry* x_p_mu_m_nu = sView.GetEntry(sIndex+3,site);
    GeneralStencilEntry* x_m_nu      = sView.GetEntry(sIndex+4,site);
    return { x_p_mu, x_p_nu, x, x_p_mu_m_nu, x_m_nu };
}


/*! @brief Retrieve 5-link stencil entries. */
template<class acc> accelerator_inline
std::tuple<GeneralStencilEntry*,GeneralStencilEntry*,GeneralStencilEntry*,
           GeneralStencilEntry*,GeneralStencilEntry*,GeneralStencilEntry*,
           GeneralStencilEntry*,GeneralStencilEntry*,GeneralStencilEntry*,
           GeneralStencilEntry*,GeneralStencilEntry*,GeneralStencilEntry*,
           GeneralStencilEntry*,GeneralStencilEntry*,GeneralStencilEntry*,
           GeneralStencilEntry*,GeneralStencilEntry*>
get5StaplePoints(acc sView, int sIndex, int site) {
    GeneralStencilEntry* x_p_nu_m_rho      = sView.GetEntry(sIndex+0 ,site);
    GeneralStencilEntry* x_p_nu            = sView.GetEntry(sIndex+1 ,site);
    GeneralStencilEntry* x_m_rho           = sView.GetEntry(sIndex+2 ,site);
    GeneralStencilEntry* x                 = sView.GetEntry(sIndex+3 ,site);
    GeneralStencilEntry* x_p_rho           = sView.GetEntry(sIndex+4 ,site);
    GeneralStencilEntry* x_m_nu_m_rho      = sView.GetEntry(sIndex+5 ,site);
    GeneralStencilEntry* x_m_nu            = sView.GetEntry(sIndex+6 ,site);
    GeneralStencilEntry* x_m_nu_p_rho      = sView.GetEntry(sIndex+7 ,site);
    GeneralStencilEntry* x_p_mu_p_nu_m_rho = sView.GetEntry(sIndex+8 ,site);
    GeneralStencilEntry* x_p_mu_p_nu       = sView.GetEntry(sIndex+9 ,site);
    GeneralStencilEntry* x_p_mu_m_rho      = sView.GetEntry(sIndex+10,site);
    GeneralStencilEntry* x_p_mu            = sView.GetEntry(sIndex+11,site);
    GeneralStencilEntry* x_p_mu_p_rho      = sView.GetEntry(sIndex+12,site);
    GeneralStencilEntry* x_p_mu_m_nu_m_rho = sView.GetEntry(sIndex+13,site);
    GeneralStencilEntry* x_p_mu_m_nu       = sView.GetEntry(sIndex+14,site);
    GeneralStencilEntry* x_p_mu_m_nu_p_rho = sView.GetEntry(sIndex+15,site);
    GeneralStencilEntry* x_p_nu_p_rho      = sView.GetEntry(sIndex+16,site);
    return {x_p_nu_m_rho     , x_p_nu           , x_m_rho          , 
            x                , x_p_rho          , x_m_nu_m_rho     , 
            x_m_nu           , x_m_nu_p_rho     , x_p_mu_p_nu_m_rho, 
            x_p_mu_p_nu      , x_p_mu_m_rho     , x_p_mu           ,
            x_p_mu_p_rho     , x_p_mu_m_nu_m_rho, x_p_mu_m_nu      , 
            x_p_mu_m_nu_p_rho, x_p_nu_p_rho}; 
}


/*! @brief Retrieve 7-link stencil entries. */
template<class acc> accelerator_inline
std::tuple<GeneralStencilEntry*,GeneralStencilEntry*,GeneralStencilEntry*,GeneralStencilEntry*,
           GeneralStencilEntry*,GeneralStencilEntry*,GeneralStencilEntry*,GeneralStencilEntry*,
           GeneralStencilEntry*,GeneralStencilEntry*,GeneralStencilEntry*,GeneralStencilEntry*,
           GeneralStencilEntry*,GeneralStencilEntry*,GeneralStencilEntry*,GeneralStencilEntry*,
           GeneralStencilEntry*,GeneralStencilEntry*,GeneralStencilEntry*,GeneralStencilEntry*,
           GeneralStencilEntry*,GeneralStencilEntry*,GeneralStencilEntry*,GeneralStencilEntry*,
           GeneralStencilEntry*,GeneralStencilEntry*,GeneralStencilEntry*,GeneralStencilEntry*,
           GeneralStencilEntry*,GeneralStencilEntry*,GeneralStencilEntry*,GeneralStencilEntry*,
           GeneralStencilEntry*,GeneralStencilEntry*,GeneralStencilEntry*,GeneralStencilEntry*,
           GeneralStencilEntry*,GeneralStencilEntry*,GeneralStencilEntry*,GeneralStencilEntry*,
           GeneralStencilEntry*,GeneralStencilEntry*,GeneralStencilEntry*,GeneralStencilEntry*,
           GeneralStencilEntry*,GeneralStencilEntry*>
get7StaplePoints(acc sView, int sIndex, int site) {
    GeneralStencilEntry* x                       = sView.GetEntry(sIndex+0 ,site); 
    GeneralStencilEntry* x_p_mu                  = sView.GetEntry(sIndex+1 ,site);
    GeneralStencilEntry* x_p_mu_p_nu             = sView.GetEntry(sIndex+2 ,site);
    GeneralStencilEntry* x_p_mu_p_nu_p_rho       = sView.GetEntry(sIndex+3 ,site);
    GeneralStencilEntry* x_p_mu_p_nu_p_rho_m_sig = sView.GetEntry(sIndex+4 ,site);
    GeneralStencilEntry* x_p_mu_p_nu_m_rho       = sView.GetEntry(sIndex+5 ,site);
    GeneralStencilEntry* x_p_mu_p_nu_m_rho_m_sig = sView.GetEntry(sIndex+6 ,site);
    GeneralStencilEntry* x_p_mu_p_nu_m_sig       = sView.GetEntry(sIndex+7 ,site);
    GeneralStencilEntry* x_p_mu_m_nu             = sView.GetEntry(sIndex+8 ,site);
    GeneralStencilEntry* x_p_mu_m_nu_p_rho       = sView.GetEntry(sIndex+9 ,site);
    GeneralStencilEntry* x_p_mu_m_nu_p_rho_p_sig = sView.GetEntry(sIndex+10,site);
    GeneralStencilEntry* x_p_mu_m_nu_p_rho_m_sig = sView.GetEntry(sIndex+11,site);
    GeneralStencilEntry* x_p_mu_m_nu_m_rho       = sView.GetEntry(sIndex+12,site);
    GeneralStencilEntry* x_p_mu_m_nu_m_rho_p_sig = sView.GetEntry(sIndex+13,site);
    GeneralStencilEntry* x_p_mu_m_nu_m_rho_m_sig = sView.GetEntry(sIndex+14,site);
    GeneralStencilEntry* x_p_mu_m_nu_p_sig       = sView.GetEntry(sIndex+15,site);
    GeneralStencilEntry* x_p_mu_m_nu_m_sig       = sView.GetEntry(sIndex+16,site);
    GeneralStencilEntry* x_p_mu_p_rho            = sView.GetEntry(sIndex+17,site);
    GeneralStencilEntry* x_p_mu_p_rho_p_sig      = sView.GetEntry(sIndex+18,site);
    GeneralStencilEntry* x_p_mu_p_rho_m_sig      = sView.GetEntry(sIndex+19,site);
    GeneralStencilEntry* x_p_mu_m_rho            = sView.GetEntry(sIndex+20,site);
    GeneralStencilEntry* x_p_mu_m_rho_p_sig      = sView.GetEntry(sIndex+21,site);
    GeneralStencilEntry* x_p_mu_m_rho_m_sig      = sView.GetEntry(sIndex+22,site);
    GeneralStencilEntry* x_p_mu_p_sig            = sView.GetEntry(sIndex+23,site); 
    GeneralStencilEntry* x_p_mu_m_sig            = sView.GetEntry(sIndex+24,site);
    GeneralStencilEntry* x_p_nu                  = sView.GetEntry(sIndex+25,site);
    GeneralStencilEntry* x_p_nu_p_rho            = sView.GetEntry(sIndex+26,site);
    GeneralStencilEntry* x_p_nu_p_rho_p_sig      = sView.GetEntry(sIndex+27,site);
    GeneralStencilEntry* x_p_nu_p_rho_m_sig      = sView.GetEntry(sIndex+28,site);
    GeneralStencilEntry* x_p_nu_m_rho            = sView.GetEntry(sIndex+29,site);
    GeneralStencilEntry* x_p_nu_m_rho_p_sig      = sView.GetEntry(sIndex+30,site);
    GeneralStencilEntry* x_p_nu_m_rho_m_sig      = sView.GetEntry(sIndex+31,site);
    GeneralStencilEntry* x_p_rho                 = sView.GetEntry(sIndex+32,site);
    GeneralStencilEntry* x_p_rho_m_nu            = sView.GetEntry(sIndex+33,site);
    GeneralStencilEntry* x_p_rho_p_sig           = sView.GetEntry(sIndex+34,site);
    GeneralStencilEntry* x_p_rho_m_sig           = sView.GetEntry(sIndex+35,site);
    GeneralStencilEntry* x_m_nu                  = sView.GetEntry(sIndex+36,site);
    GeneralStencilEntry* x_m_nu_p_rho            = sView.GetEntry(sIndex+37,site);
    GeneralStencilEntry* x_m_nu_p_rho_p_sig      = sView.GetEntry(sIndex+38,site);
    GeneralStencilEntry* x_m_nu_p_rho_m_sig      = sView.GetEntry(sIndex+39,site);
    GeneralStencilEntry* x_m_nu_m_rho            = sView.GetEntry(sIndex+40,site);
    GeneralStencilEntry* x_m_nu_m_rho_p_sig      = sView.GetEntry(sIndex+41,site);
    GeneralStencilEntry* x_m_nu_m_rho_m_sig      = sView.GetEntry(sIndex+42,site);
    GeneralStencilEntry* x_m_rho                 = sView.GetEntry(sIndex+43,site);
    GeneralStencilEntry* x_m_rho_p_sig           = sView.GetEntry(sIndex+44,site);
    GeneralStencilEntry* x_m_rho_m_sig           = sView.GetEntry(sIndex+45,site);
    return {x                      , x_p_mu                 , x_p_mu_p_nu            , x_p_mu_p_nu_p_rho      ,
            x_p_mu_p_nu_p_rho_m_sig, x_p_mu_p_nu_m_rho      , x_p_mu_p_nu_m_rho_m_sig, x_p_mu_p_nu_m_sig      ,
            x_p_mu_m_nu            , x_p_mu_m_nu_p_rho      , x_p_mu_m_nu_p_rho_p_sig, x_p_mu_m_nu_p_rho_m_sig,
            x_p_mu_m_nu_m_rho      , x_p_mu_m_nu_m_rho_p_sig, x_p_mu_m_nu_m_rho_m_sig, x_p_mu_m_nu_p_sig      ,
            x_p_mu_m_nu_m_sig      , x_p_mu_p_rho           , x_p_mu_p_rho_p_sig     , x_p_mu_p_rho_m_sig     , 
            x_p_mu_m_rho           , x_p_mu_m_rho_p_sig     , x_p_mu_m_rho_m_sig     , x_p_mu_p_sig           , 
            x_p_mu_m_sig           , x_p_nu                 , x_p_nu_p_rho           , x_p_nu_p_rho_p_sig     , 
            x_p_nu_p_rho_m_sig     , x_p_nu_m_rho           , x_p_nu_m_rho_p_sig     , x_p_nu_m_rho_m_sig     , 
            x_p_rho                , x_p_rho_m_nu           , x_p_rho_p_sig          , x_p_rho_m_sig          , 
            x_m_nu                 , x_m_nu_p_rho           , x_m_nu_p_rho_p_sig     , x_m_nu_p_rho_m_sig     , 
            x_m_nu_m_rho           , x_m_nu_m_rho_p_sig     , x_m_nu_m_rho_m_sig     , x_m_rho                , 
            x_m_rho_p_sig          , x_m_rho_m_sig};
}


/*!  @brief Allows for ASQTAD-like smearings. */
template<class Gimpl> 
class Smear_HISQ : public Gimpl {
public:

    GridCartesian* const _grid;

    // Sort out the Gimpl. This handles BCs and part of the precision. 
    INHERIT_GIMPL_TYPES(Gimpl);
    typedef typename Gimpl::GaugeField      GF;
    typedef typename Gimpl::GaugeLinkField  LF;
    typedef typename Gimpl::ComplexField    CF;
    typedef typename Gimpl::Scalar          ComplexScalar;
    typedef decltype(real(ComplexScalar())) RealScalar;
    typedef iColourMatrix<ComplexScalar>    ComplexColourMatrix;

    RealScalar _Scut; // Cutoff for U(3) projection eigenvalues, set at initialization
    int _HaloDepth=1; 

    SmearingParameters<RealScalar> _linkTreatment;

    void initialize() {
        if (sizeof(RealScalar)==4) {
            _Scut=1e-5; // Maybe should be higher? e.g. 1e-4
        } else if (sizeof(RealScalar)==8) { 
            _Scut=1e-8;
        } else {
            Grid_error("HISQ smearing only implemented for single and double");
        }
        assert(Nc == 3 && "HISQ smearing currently implemented only for Nc==3");
        assert(Nd == 4 && "HISQ smearing only defined for Nd==4");
    }

    Smear_HISQ(GridCartesian* grid, RealScalar c1, RealScalar cnaik, RealScalar c3, RealScalar c5, RealScalar c7, RealScalar clp) 
        : _grid(grid), 
          _linkTreatment(c1,cnaik,c3,c5,c7,clp) {
        initialize();
    }

    // Allow to pass a pointer to a C-style array for MILC convenience
    Smear_HISQ(GridCartesian* grid, double* coeff) 
        : _grid(grid), 
          _linkTreatment(coeff[0],coeff[1],coeff[2],coeff[3],coeff[4],coeff[5]) {
        initialize();
    }
    Smear_HISQ(GridCartesian* grid, float* coeff) 
        : _grid(grid), 
          _linkTreatment(coeff[0],coeff[1],coeff[2],coeff[3],coeff[4],coeff[5]) {
        initialize();
    }

    ~Smear_HISQ() {}


    // Intent: OUT--U_3link (sum of left and right 3-staples attached to U)
    //              U_fat (accmulates the fat smearing)
    //          IN--U (thin links)
    //              gStencil (3-link stencil)
    //              mu
    void threeLinkStaple(GF& U_fat, GF& U_3link, GF& U, GeneralLocalStencil gStencil, int mu) const {

        SmearingParameters<RealScalar> lt = this->_linkTreatment;
        autoView(U_v       , U       , AcceleratorRead);
        autoView(U_fat_v   , U_fat   , AcceleratorWrite);
        autoView(U_3link_v , U_3link , AcceleratorWrite);
        auto gStencil_v = gStencil.View(AcceleratorRead); 
        typedef decltype(getLink(U_v,gStencil_v.GetEntry(0,0),0)) U3matrix;
        int Nsites = U_v.size();

        accelerator_for(site,Nsites,Simd::Nsimd(),{
            U3matrix res;
            for(int nu=0;nu<Nd;nu++) {
                if(nu==mu) continue;
                int s = HISQStencilIndex(mu,nu);

                auto [x_p_mu, x_p_nu, x, x_p_mu_m_nu, x_m_nu] = get3StaplePoints(gStencil_v,s,site);      

                // When you're deciding whether to take an adjoint, the question is: how is the
                // stored link oriented compared to the one you want? If I imagine myself traveling
                // with the to-be-updated link, I have two possible, alternative 3-link paths I can
                // take, one starting by going to the left, the other starting by going to the right.

                //  "left" + "right"
                res =       getLink(U_v,x,nu)       * getLink(U_v,x_p_nu,mu) * adj(getLink(U_v,x_p_mu,nu))
                      + adj(getLink(U_v,x_m_nu,nu)) * getLink(U_v,x_m_nu,mu) *     getLink(U_v,x_p_mu_m_nu,nu);

                // Save 3-link construct for later 
                setLink(U_3link_v[x->_offset](nu), res);

                // The index operator (x) returns the coalesced read on GPU. The view [] index returns 
                // a reference to the vector object. The [x](mu) returns a reference to the densely 
                // packed (contiguous in memory) mu-th element of the vector object. 
                setLink(U_fat_v[x->_offset](mu), U_fat_v(x->_offset)(mu) + lt.c_3*res);
            }
        })
        return;
    }


    // Intent: OUT--U_5link (sum of 5-staples attached to U) 
    //              U_fat (accmulates the fat smearing)
    //          IN--U (thin links)
    //              U_3link (sum of left and right 3-staples attached to U)
    //              gStencil (3-link stencil)
    //              mu
    //              updateFatLinks (in the force, you only want U_5link_v)
    void fiveLinkStaple(GF& U_fat, GF& U_5linkA, GF& U_5linkB, GF& U_3link, 
                        GF& U, GeneralLocalStencil gStencil, int mu) const {

        SmearingParameters<RealScalar> lt = this->_linkTreatment;
        autoView(U_v       , U       , AcceleratorRead);
        autoView(U_fat_v   , U_fat   , AcceleratorWrite);
        autoView(U_3link_v , U_3link , AcceleratorWrite);
        autoView(U_5linkA_v, U_5linkA, AcceleratorWrite);
        autoView(U_5linkB_v, U_5linkB, AcceleratorWrite);
        auto gStencil_v = gStencil.View(AcceleratorRead); 
        typedef decltype(getLink(U_v,gStencil_v.GetEntry(0,0),0)) U3matrix;
        int Nsites = U_v.size();

        accelerator_for(site,Nsites,Simd::Nsimd(),{
            U3matrix res;
            int sigmaIndex = 0;
            for(int nu=0;nu<Nd;nu++) {
                if(nu==mu) continue;
                int s = HISQStencilIndex(mu,nu);
                for(int rho=0;rho<Nd;rho++) {
                    if (rho == mu || rho == nu) continue;

                    auto [x_p_mu, x_p_nu, x, x_p_mu_m_nu, x_m_nu] = get3StaplePoints(gStencil_v,s,site);      

                    res =       getLink(U_v,x,nu )      * getLink(U_3link_v,x_p_nu,rho) * adj(getLink(U_v,x_p_mu,nu))
                          + adj(getLink(U_v,x_m_nu,nu)) * getLink(U_3link_v,x_m_nu,rho) *     getLink(U_v,x_p_mu_m_nu,nu);

                    // Counting 3-link staples: there are three planes attached to the to-be-updated link,
                    // which corresponds to three (forward+backward) staples. For the 5-link staples, for
                    // each plane, there are two remaining directions, so that there are six 5-link staples
                    // altogether. That will not fit in a single GaugeField object, so we use two. You can
                    // think of sigmaIndex and rho together as being the labels that pick out a particular
                    // 5-link staple. They therefore should not be interpreted as directions.
                    if(sigmaIndex<3)
                        setLink(U_5linkA_v[x->_offset](rho), res);
                    else
                        setLink(U_5linkB_v[x->_offset](rho), res);

                    setLink(U_fat_v[x->_offset](mu), U_fat_v(x->_offset)(mu) + lt.c_5*res);
                    sigmaIndex++;
                }
            }
        })
        return;
    }


    // Intent: OUT--U_fat (accmulates the fat smearing)
    //          IN--U (thin links)
    //              U_5link (sum of 5-staples attached to U)
    //              gStencil (3-link stencil)
    //              mu
    void sevenLinkStaple(GF& U_fat, GF& U_5linkA, GF& U_5linkB, GF& U, GeneralLocalStencil gStencil, int mu) const {

        SmearingParameters<RealScalar> lt = this->_linkTreatment;
        autoView(U_v       , U       , AcceleratorRead);
        autoView(U_fat_v   , U_fat   , AcceleratorWrite);
        autoView(U_5linkA_v, U_5linkA, AcceleratorWrite);
        autoView(U_5linkB_v, U_5linkB, AcceleratorWrite);
        auto gStencil_v = gStencil.View(AcceleratorRead); 
        typedef decltype(getLink(U_v,gStencil_v.GetEntry(0,0),0)) U3matrix;
        int Nsites = U_v.size();

        accelerator_for(site,Nsites,Simd::Nsimd(),{ 
            U3matrix res;
            int sigmaIndex = 0;
            for(int nu=0;nu<Nd;nu++) {
                if(nu==mu) continue;
                int s = HISQStencilIndex(mu,nu);
                for(int rho=0;rho<Nd;rho++) {
                    if (rho == mu || rho == nu) continue;
    
                    auto [x_p_mu, x_p_nu, x, x_p_mu_m_nu, x_m_nu] = get3StaplePoints(gStencil_v,s,site);      
    
                    if(sigmaIndex<3)
                        res =       getLink(U_v,x,nu)       * getLink(U_5linkB_v,x_p_nu,rho) * adj(getLink(U_v,x_p_mu,nu))
                              + adj(getLink(U_v,x_m_nu,nu)) * getLink(U_5linkB_v,x_m_nu,rho) *     getLink(U_v,x_p_mu_m_nu,nu);
                    else
                        res =       getLink(U_v,x,nu)       * getLink(U_5linkA_v,x_p_nu,rho) * adj(getLink(U_v,x_p_mu,nu))
                              + adj(getLink(U_v,x_m_nu,nu)) * getLink(U_5linkA_v,x_m_nu,rho) *     getLink(U_v,x_p_mu_m_nu,nu);
                    
                    setLink(U_fat_v[x->_offset](mu), U_fat_v(x->_offset)(mu) + lt.c_7*res);
                    sigmaIndex++;
                }
            }
        })
        return;
    }


    // Intent: OUT--u_smr (smeared links) 
    //              u_naik (Naik links)
    //          IN--u_thin (thin links)
    void smear(GF& u_smr, GF& u_naik, GF& u_thin) const {

        SmearingParameters<RealScalar> lt = this->_linkTreatment;
        auto grid = this->_grid;

        // Create a padded cell of extra padding depth=1 and fill the padding.
        PaddedCell Ghost(_HaloDepth,grid);
        GF Ughost = Ghost.Exchange(u_thin);

        // This is where auxiliary N-link fields and the final smear will be stored. As
        // implemented, this uses about 25% more memory than necessary. 
        GF U_fat(Ughost.Grid());
        GF U_3link(Ughost.Grid());
        GF U_5linkA(Ughost.Grid());
        GF U_5linkB(Ughost.Grid());

        // mu-nu plane stencil.
        std::vector<Coordinate> shifts = createHISQStencil("3STAPLE");

        // A GeneralLocalStencil has two indices: a site and stencil index 
        GeneralLocalStencil gStencil(Ughost.Grid(),shifts);

        // Store sum of 3-, 5-, 7-link contributions in U_fat 
        U_fat=Zero();
        for(int mu=0;mu<Nd;mu++) {
            U_3link=Zero(); U_5linkA=Zero(); U_5linkB=Zero();
            if((lt.c_7!=0) || (lt.c_5!=0) || (lt.c_3!=0))
                threeLinkStaple(U_fat,                     U_3link, Ughost, gStencil, mu);
            if((lt.c_7!=0) || (lt.c_5!=0)) 
                fiveLinkStaple( U_fat, U_5linkA, U_5linkB, U_3link, Ughost, gStencil, mu);
            if( lt.c_7!=0 ) 
                sevenLinkStaple(U_fat, U_5linkA, U_5linkB,          Ughost, gStencil, mu);
        }

        // Add 1-link contribution 
        u_smr = Ghost.Extract(U_fat) + lt.c_1*u_thin;

        // Load up U and V std::vectors to access thin and smeared links.
        std::vector<LF> U(Nd, grid);
        std::vector<LF> V(Nd, grid);
        std::vector<LF> Vnaik(Nd, grid);
        for (int mu = 0; mu < Nd; mu++) {
            U[mu]     = PeekIndex<LorentzIndex>(u_thin, mu);
            V[mu]     = PeekIndex<LorentzIndex>(u_smr, mu);
            Vnaik[mu] = Zero();
        }

        for(int mu=0;mu<Nd;mu++) {

            if(lt.c_naik!=0) { // Naik
                Vnaik[mu] = lt.c_naik*Gimpl::CovShiftForward(U[mu],mu,
                                        Gimpl::CovShiftForward(U[mu],mu,
                                          Gimpl::CovShiftIdentityForward(U[mu],mu)));
            }

            if(lt.c_lp!=0) { // LePage
                for (int nu=0;nu<Nd;nu++) {
                    if(mu==nu) continue; 
                    V[mu] = V[mu] + lt.c_lp*Gimpl::CovShiftForward(U[nu],nu,
                                              Gimpl::CovShiftForward(U[nu],nu,
                                                Gimpl::CovShiftForward(U[mu],mu,
                                                  Gimpl::CovShiftBackward(U[nu],nu,
                                                    Gimpl::CovShiftIdentityBackward(U[nu],nu)))))
                                  + lt.c_lp*Gimpl::CovShiftBackward(U[nu],nu,
                                              Gimpl::CovShiftBackward(U[nu],nu,
                                                Gimpl::CovShiftForward(U[mu],mu,
                                                  Gimpl::CovShiftForward(U[nu],nu,
                                                    Gimpl::CovShiftIdentityForward(U[nu],nu)))));
                }
            }
        }

        // Put V back into u_smr.
        for (int mu = 0; mu < Nd; mu++) {
            PokeIndex<LorentzIndex>(u_smr , V[mu]    , mu);
            PokeIndex<LorentzIndex>(u_naik, Vnaik[mu], mu);
        }
    };


    // Intent: OUT--u_proj (U3-projected links)
    //          IN--u_mu (to-be-projected links)
    void projectU3(GF& u_proj, GF& u_mu) const {

        // Open up the views
        autoView(uproj_v, u_proj, AcceleratorWrite);
        autoView(umu_v  , u_mu  , AcceleratorRead);

        // Make sure everyone is using the same Grid
        conformable(u_proj,u_mu);

        // Follow MILC 10.1103/PhysRevD.82.074501, eqs (B2-B3) and (C1-C8)
        accelerator_for(ss,umu_v.size(),Simd::Nsimd(),{
#ifdef GRID_SIMT
            { int blane=acceleratorSIMTlane(Simd::Nsimd());//
#else
            for(int blane=0;blane<Simd::Nsimd();blane++) {
#endif
                RealScalar g1, g2, g0;
                ComplexColourMatrix V;
                auto Vmu = extractLane(blane,umu_v[ss]);
                for (int mu = 0; mu < Nd; mu++) {
                    V()     = Vmu(mu);
                    auto Q  = adj(V)*V;
                    RealScalar c0 =        real(trace(Q))()()();
                    RealScalar c1 = (1/2.)*real(trace(Q*Q))()()();
                    RealScalar c2 = (1/3.)*real(trace(Q*Q*Q))()()();
                    RealScalar S  = (1/3.)*c1-(1/18.)*c0*c0;
                    if (abs(S)<_Scut) {
                        g0 = (1/3.)*c0; 
                        g1 = g0; 
                        g2 = g1;
                    } else {
                        RealScalar R     = (1/2.)*c2-(1/3. )*c0*c1+(1/27.)*c0*c0*c0;
                        RealScalar theta = acos(R*pow(S,-1.5));
                        g0 = (1/3.)*c0+2.*sqrt(S)*cos((1/3.)*theta-2*M_PI/3.);
                        g1 = (1/3.)*c0+2.*sqrt(S)*cos((1/3.)*theta          );
                        g2 = (1/3.)*c0+2.*sqrt(S)*cos((1/3.)*theta+2*M_PI/3.);
                    }
                    RealScalar u   = sqrt(g0) + sqrt(g1) + sqrt(g2);
                    RealScalar v   = sqrt(g0*g1) + sqrt(g0*g2) + sqrt(g1*g2);
                    RealScalar w   = sqrt(g0*g1*g2);
                    RealScalar den = w*(u*v-w);
                    RealScalar f0  = (-w*(u*u+v)+u*v*v)/den;
                    RealScalar f1  = (-w-u*u*u+2.*u*v)/den;
                    RealScalar f2  = u/den;

                    auto res = V*(f0 + f1*Q + f2*Q*Q);

                    insertLane(blane,uproj_v[ss](mu),res());
                }
            } 
        });
    };
};



/*!  @brief compute force from link variables */
template<class Gimpl> 
class Force_HISQ : public Gimpl {
public:

    GridCartesian* const   _grid;
    GridRedBlackCartesian* _gridRB;

    // Sort out the Gimpl. This handles BCs and part of the precision. 
    INHERIT_GIMPL_TYPES(Gimpl);
    typedef typename Gimpl::FermionField   FF;
    typedef typename Gimpl::GaugeField     GF;
    typedef typename Gimpl::GaugeLinkField LF;
    typedef typename Gimpl::ComplexField   CF;
    typedef typename Gimpl::Scalar ComplexScalar;
    typedef decltype(real(ComplexScalar())) RealScalar;
    typedef iColourMatrix<ComplexScalar> ComplexColourMatrix;

    RealScalar _Scut=-1; // Cutoff for U(3) projection eigenvalues, set at initialization
    int _HaloDepth=1;    // Depth of padded cell 

    HISQParameters<Real> _linkParams;
    HISQReunitSVDParameters<Real> _reunitParams;
    GF _Umu, _Vmu, _Wmu;

    void initialize() {
        if (sizeof(RealScalar)==4) {
            _Scut=1e-5; // Maybe should be higher? e.g. 1e-4
        } else if (sizeof(RealScalar)==8) { 
            _Scut=1e-8;
        } else {
            Grid_error("HISQ force only implemented for single and double");
        }
        _gridRB = SpaceTimeGrid::makeFourDimRedBlackGrid(_grid);
        assert(Nc == 3 && "HISQ force currently implemented only for Nc==3");
        assert(Nd == 4 && "HISQ force only defined for Nd==4");
    }

    Force_HISQ(GridCartesian* grid, HISQParameters<Real> linkParams, GF Wmu, GF Vmu, GF Umu, 
               HISQReunitSVDParameters<Real> reunitParams) 
        : _grid(grid), 
          _linkParams(linkParams),
          _Wmu(Wmu),
          _Vmu(Vmu),
          _Umu(Umu),
          _reunitParams(reunitParams) {
        initialize();
    }

    ~Force_HISQ() {}


    // Intent: OUT--u_deriv (dW/dV slotted into force)
    //          IN--u_mu (fat links), 
    //              u_force (slot derivative into this force), 
    //              delta (force cutoff)
    // Follow MILC 10.1103/PhysRevD.82.074501
    void projU3Deriv(GF& u_deriv, GF& u_mu, GF& u_force, RealScalar const delta=5e-5) {

        conformable(u_force,u_mu);
        conformable(u_deriv,u_mu);

        autoView(uderiv_v, u_deriv, AcceleratorWrite);
        autoView(umu_v   , u_mu   , AcceleratorRead);
        autoView(uforce_v, u_force, AcceleratorRead);

        // Follow MILC 10.1103/PhysRevD.82.074501, eqs (B2-B3) and (C1-C8)
        accelerator_for(ss,umu_v.size(),Simd::Nsimd(),{
#ifdef GRID_SIMT
            { int blane=acceleratorSIMTlane(Simd::Nsimd());//
#else
            for(int blane=0;blane<Simd::Nsimd();blane++) {
#endif
                RealScalar g1, g2, g0;
                ComplexColourMatrix V, force;
                auto Vmu     = extractLane(blane,umu_v[ss]);
                auto forcemu = extractLane(blane,uforce_v[ss]);
                for (int mu = 0; mu < Nd; mu++) {
                    V()     = Vmu(mu);
                    auto Q  = adj(V)*V;
                    RealScalar c0 =        real(trace(Q))()()();
                    RealScalar c1 = (1/2.)*real(trace(Q*Q))()()();
                    RealScalar c2 = (1/3.)*real(trace(Q*Q*Q))()()();
                    RealScalar S  = (1/3.)*c1-(1/18.)*c0*c0;
                    if (abs(S)<_Scut) {
                        g0 = (1/3.)*c0; 
                        g1 = g0; 
                        g2 = g1;
                    } else {
                        RealScalar R     = (1/2.)*c2-(1/3. )*c0*c1+(1/27.)*c0*c0*c0;
                        RealScalar theta = acos(R*pow(S,-1.5));
                        g0 = (1/3.)*c0+2.*sqrt(S)*cos((1/3.)*theta-2*M_PI/3.);
                        g1 = (1/3.)*c0+2.*sqrt(S)*cos((1/3.)*theta          );
                        g2 = (1/3.)*c0+2.*sqrt(S)*cos((1/3.)*theta+2*M_PI/3.);
                    }

                    if (g0 < delta || g1 < delta || g2 < delta) {
                        // force filter eq (C23)
                        g0 += delta;
                        g1 += delta;
                        g2 += delta;
                        Q   = Q + delta;
                    }
//                if (fabs(Q.determinant()/(g0*g1*g2)-1.0) > 1e-5) { SVD }

                    RealScalar u   = sqrt(g0) + sqrt(g1) + sqrt(g2);
                    RealScalar v   = sqrt(g0*g1) + sqrt(g0*g2) + sqrt(g1*g2);
                    RealScalar w   = sqrt(g0*g1*g2);
                    RealScalar den = w*(u*v-w);
                    RealScalar f0  = (-w*(u*u+v)+u*v*v)/den;
                    RealScalar f1  = (-w-u*u*u+2.*u*v)/den;
                    RealScalar f2  = u/den;

                    auto Qinvsq = f0 + f1*Q + f2*Q*Q;

                    force() = forcemu(mu);
                    auto forcedag = adj(force);

                    RealScalar u2, u3, u4, u5, u6, u7, u8, v2 ,v3, v4, v5, v6, w2, w3, w4, w5;

                    u2 = u*u; u3 = u2*u; u4 = u3*u; u5 = u4*u; u6 = u5*u; u7 = u6*u; u8 = u7*u;
                    v2 = v*v; v3 = v2*v; v4 = v3*v; v5 = v4*v; v6 = v5*v;
                    w2 = w*w; w3 = w2*w; w4 = w3*w; w5 = w4*w;
        
                    // eq (C10)
                    auto d = 2*w3*(u*v-w)*(u*v-w)*(u*v-w);
        
                    // eq (C11)
                    auto C00 = ( -w3*u6 + 3*v*w3*u4 + 3*v4*w*u4 - v6*u3 - 4*w4*u3 - 12*v3*w2*u3 + 16*v2*w3*u2 
                                 + 3*v5*w*u2 - 8*v*w4*u - 3*v4*w2*u + w5 + v3*w3 )/d;
                    auto C01 = ( -w2*u7 - v2*w*u6 + v4*u5 + 6*v*w2*u5 - 5*w3*u4 - v3*w*u4 - 2*v5*u3 - 6*v2*w2*u3 
                                 + 10*v*w3*u2 + 6*v4*w*u2 - 3*w4*u - 6*v3*w2*u + 2*v2*w3 )/d;
                    auto C02 = ( w2*u5 + v2*w*u4 - v4*u3 - 4*v*w2*u3 + 4*w3*u2 +3*v3*w*u2 - 3*v2*w2*u + v*w3 )/d;
                    auto C11 = ( -w*u8 - v2*u7 + 7*v*w*u6 + 4*v3*u5 - 5*w2*u5 - 16*v2*w*u4 - 4*v4*u3 + 16*v*w2*u3 
                                 - 3*w3*u2 + 12*v3*w*u2 - 12*v2*w2*u + 3*v*w3 )/d;
                    auto C12 = ( w*u6 + v2*u5 - 5*v*w*u4 - 2*v3*u3 + 4*w2*u3 + 6*v2*w*u2 - 6*v*w2*u + w3 )/d;
                    auto C22 = ( -w*u4 - v2*u3 + 3*v*w*u2 - 3*w2*u )/d;
        
                    // These are all used in the loop over color entries, and we want to avoid recomputing
                    // these products, which should be broadcast to all sites, 3*3*3*3=81 times. 
                    auto Vdag   = adj(V);
                    auto VVdag  = V*Vdag;
                    auto VQ     = V*Q;   
                    auto VQ2    = VQ*Q;
                    auto VQVdag = VQ*Vdag;
                    auto QVdag  = Q*Vdag;
                    auto Q2Vdag = Q*QVdag;
        
                    // eqs (C17-C19)
                    auto PVdag  = ( C00 + C01*Q + C02*Q*Q )*Vdag;
                    auto RVdag  = ( C01 + C11*Q + C12*Q*Q )*Vdag;
                    auto SVdag  = ( C02 + C12*Q + C22*Q*Q )*Vdag;
        
                    // eqs (C20) and (C21)
                    ComplexColourMatrix res = Zero();
                    for (int k = 0; k < 3; k++)
                    for (int l = 0; l < 3; l++)
                    for (int i = 0; i < 3; i++)
                    for (int j = 0; j < 3; j++) {
        
                        ComplexScalar deriv = 0.; // dWij/dVkl
                
                        if (k == i) deriv += Qinvsq()()(l,j);
                        if (l == j) deriv += f1*VVdag()()(i,k)+f2*VQVdag()()(i,k);

                        deriv += f2*VVdag()()(i,k)*Q()()(l,j) + V()()(i,j)*PVdag()()(l,k) 
                                 + VQ()()(i,j)*RVdag()()(l,k) + VQ2()()(i,j)*SVdag()()(l,k);
        
                        res()()(l,k) = res()()(l,k) + deriv*force()()(j,i);
                
                        // dWij^+/dVkl
                        deriv = (f1*Vdag()()(i,k)+f2*QVdag()()(i,k))*Vdag()()(l,j) 
                                + f2*Vdag()()(i,k)*QVdag()()(l,j) + Vdag()()(i,j)*PVdag()()(l,k) 
                                + QVdag()()(i,j)*RVdag()()(l,k)+Q2Vdag()()(i,j)*SVdag()()(l,k);
                
                        res()()(l,k) = res()()(l,k) + deriv*forcedag()()(j,i);
                  	}
        
                    insertLane(blane,uderiv_v[ss](mu),res());
                }
            } 
        });
    };


    // vecx (contains |X> and |Y>)
    // l (rat approx and Naik index)
    // sep (separation between |X> and |Y>)
    GF outerProductHISQ(std::vector<FF>& vecx, std::vector<Real> vecdt, std::vector<int> n_orders_naik, int n_naiks, int sep) {
        
        auto grid   = this->_grid;
        auto gridRB = this->_gridRB;

        GF XY(grid), XY_l(grid);
        FF X(grid), Y(grid), RB(gridRB);
        LF XYnu(grid), YXnu(grid);

        XY = Zero();

        // These four lines control the loop over rational approximation contributions. As explained in force(), 
        // l indexes over both Naik epsilon and rational approximation order.
        int l = 0;
        for (int inaik = 0; inaik < n_naiks; inaik++) {
            int rat_order = n_orders_naik[inaik];
            for (int i=0; i<rat_order; i++) {

                X = Zero(); Y = Zero(); 
                
                RB = Zero(); 
                pickCheckerboard(Even,RB,vecx[l]);
                setCheckerboard(X,RB);
                RB = Zero();
                pickCheckerboard(Odd ,RB,vecx[l]);
                setCheckerboard(Y,RB);
        
                XY_l = Zero(); XYnu = Zero(); YXnu = Zero(); 
                for (int nu = 0; nu < Nd; nu++) {
                    YXnu = outerProduct( Cshift(Y,nu,sep) ,X);
                    XYnu = outerProduct( Cshift(X,nu,sep) ,Y);
                    PokeIndex<LorentzIndex>(XY_l,(YXnu-XYnu),nu);
                }
                XY += vecdt[l]*XY_l; 
                l++;
            }
        }   
        return XY;
    }


    // Intent: OUT--Fghost (accumulates 3-link derivative contribution)
    //          IN--Ughost (thin links)
    //              XYCghost (outer product)
    //              gStencil3 (3-link stencil)
    //              c3
    //              mu
    void threeLinkDeriv(GF& Fghost, GF& Ughost, GF& XYCghost, GeneralLocalStencil gStencil3, Real c3, int mu) const {
        
        autoView(U_v , Ughost , AcceleratorRead);
        autoView(XY_v, XYCghost, AcceleratorRead);
        autoView(F_v , Fghost , AcceleratorWrite);
        int Nsites = U_v.size();
        auto gStencil3_v = gStencil3.View(AcceleratorRead);
        typedef decltype(getLink(U_v,gStencil3.GetEntry(0,0),0)) U3matrix;

        accelerator_for(site,Nsites,Simd::Nsimd(),{
            U3matrix res;
            for(int nu=0;nu<Nd;nu++) {
                if(nu==mu) continue;
                int s = HISQStencilIndex(mu,nu);

                auto [x_p_mu, x_p_nu, x, x_p_mu_m_nu, x_m_nu] = get3StaplePoints(gStencil3_v,s,site);      

                res =   adj(getLink(XY_v,x,nu))*    getLink(U_v ,x_p_nu,mu) *adj(getLink(U_v ,x_p_mu,nu))
                      +     getLink(U_v ,x,nu) *adj(getLink(XY_v,x_p_nu,mu))*adj(getLink(U_v ,x_p_mu,nu))
                      +     getLink(U_v ,x,nu) *    getLink(U_v ,x_p_nu,mu) *    getLink(XY_v,x_p_mu,nu)

                      +     getLink(XY_v,x_m_nu,nu) *    getLink(U_v ,x_m_nu,mu) *    getLink(U_v ,x_p_mu_m_nu,nu)
                      + adj(getLink(U_v ,x_m_nu,nu))*adj(getLink(XY_v,x_m_nu,mu))*    getLink(U_v ,x_p_mu_m_nu,nu)
                      + adj(getLink(U_v ,x_m_nu,nu))*    getLink(U_v ,x_m_nu,mu) *adj(getLink(XY_v,x_p_mu_m_nu,nu));

                setLink(F_v[x->_offset](mu), F_v(x->_offset)(mu) + c3*adj(res));
            }              
        })
    }


    // Intent: OUT--Fghost (accumulates 5-link derivative contribution)
    //          IN--Ughost (thin links)
    //              XYCghost (outer product)
    //              gStencil5 (5-link stencil)
    //              c5
    //              mu
    template<int term>
    void fiveLinkDeriv(GF& Fghost, GF& Ughost, GF& XYCghost, GeneralLocalStencil gStencil5, Real c5, int mu) const {
        
        autoView(U_v , Ughost , AcceleratorRead);
        autoView(XY_v, XYCghost, AcceleratorRead);
        autoView(F_v , Fghost , AcceleratorWrite);
        int Nsites = U_v.size();
        auto gStencil5_v = gStencil5.View(AcceleratorRead);
        typedef decltype(getLink(U_v,gStencil5.GetEntry(0,0),0)) U3matrix;

        accelerator_for(site,Nsites,Simd::Nsimd(),{ 
            U3matrix res; 
            for(int nu=0;nu<Nd;nu++) {
                if(nu==mu) continue;
                for(int rho=0;rho<Nd;rho++) {
                    if (rho == mu || rho == nu) continue;
                    int s = HISQStencilIndex(mu,nu,rho,0,"5STAPLE");

                    auto [x_p_nu_m_rho     , x_p_nu           , x_m_rho          , 
                          x                , x_p_rho          , x_m_nu_m_rho     , 
                          x_m_nu           , x_m_nu_p_rho     , x_p_mu_p_nu_m_rho, 
                          x_p_mu_p_nu      , x_p_mu_m_rho     , x_p_mu           ,
                          x_p_mu_p_rho     , x_p_mu_m_nu_m_rho, x_p_mu_m_nu      , 
                          x_p_mu_m_nu_p_rho, x_p_nu_p_rho     ] = get5StaplePoints(gStencil5_v,s,site);

                    res = Zero();

                    // The idea behind the constexpr syntax is to reduce compile times. These seem to grow
                    // with increasing kernel size. The template parameter term lets the user choose which
                    // part of the kernel to compile, and hence constexpr is evaluated at compile time.
if constexpr(term==0) {
                    res += (      getLink(U_v ,x_p_mu           ,rho)
                            *     getLink(U_v ,x_p_mu_p_rho     ,nu ) 
                            * adj(getLink(U_v ,x_p_mu_p_nu      ,rho))

                            + adj(getLink(U_v ,x_p_mu_m_rho     ,rho))
                            *     getLink(U_v ,x_p_mu_m_rho     ,nu )
                            *     getLink(U_v ,x_p_mu_p_nu_m_rho,rho)
                    )*adj(getLink(U_v,x_p_nu,mu))*getLink(XY_v,x,nu);

                    res += (      getLink(U_v ,x_p_mu           ,rho) 
                            * adj(getLink(U_v ,x_p_rho          ,mu )) 
                            * adj(getLink(U_v ,x_m_nu_p_rho     ,nu ))
                            *     getLink(XY_v,x_m_nu           ,rho)

                            + adj(getLink(U_v ,x_p_mu_m_rho     ,rho))
                            * adj(getLink(U_v ,x_m_rho          ,mu ))
                            * adj(getLink(U_v ,x_m_nu_m_rho     ,nu ))
                            * adj(getLink(XY_v,x_m_nu_m_rho     ,rho))
                    )*getLink(U_v,x_m_nu,nu);

                    res += (      getLink(U_v ,x_p_mu           ,rho) 
                            * adj(getLink(U_v ,x_p_mu_m_nu_p_rho,nu ))
                            * adj(getLink(U_v ,x_p_mu_m_nu      ,rho))
                            + adj(getLink(U_v ,x_p_mu_m_rho     ,rho))
                            * adj(getLink(U_v ,x_p_mu_m_nu_m_rho,nu ))
                            *     getLink(U_v ,x_p_mu_m_nu_m_rho,rho)
                    )*adj(getLink(U_v,x_m_nu,mu))*adj(getLink(XY_v,x_m_nu,nu));

                    res += (      getLink(U_v ,x_p_mu           ,rho) 
                            * adj(getLink(U_v ,x_p_rho          ,mu ))
                            *     getLink(U_v ,x_p_rho          ,nu )
                            *     getLink(XY_v,x_p_nu           ,rho) 
                            + adj(getLink(U_v ,x_p_mu_m_rho     ,rho))
                            * adj(getLink(U_v ,x_m_rho          ,mu ))
                            *     getLink(U_v ,x_m_rho          ,nu )
                            * adj(getLink(XY_v,x_p_nu_m_rho     ,rho))
                    )*adj(getLink(U_v,x,nu));

                    res += (      getLink(XY_v,x_p_mu_m_nu      ,nu ) 
                            * adj(getLink(U_v ,x_m_nu           ,mu ))
                            *     getLink(U_v ,x_m_nu           ,rho)
                            +     getLink(U_v ,x_p_mu           ,rho)
                            * adj(getLink(U_v ,x_p_mu_m_nu_p_rho,nu ))
                            *     getLink(XY_v,x_m_nu_p_rho     ,mu )
                    )*getLink(U_v,x_m_nu_p_rho,nu)*adj(getLink(U_v,x,rho));
}
if constexpr(term==1) {
                    res += (      getLink(U_v ,x_p_mu           ,nu ) 
                            * adj(getLink(XY_v,x_p_mu_p_nu      ,rho))
                            * adj(getLink(U_v ,x_p_mu_p_rho     ,nu ))
                            + adj(getLink(U_v ,x_p_mu_m_nu      ,nu ))
                            * adj(getLink(XY_v,x_p_mu_m_nu      ,rho))
                            *     getLink(U_v ,x_p_mu_m_nu_p_rho,nu )
                    )*adj(getLink(U_v,x_p_rho,mu))*adj(getLink(U_v,x,rho));

                    res += (      getLink(U_v ,x_p_mu           ,rho) 
                            *     getLink(U_v ,x_p_mu_p_rho     ,nu )
                            *     getLink(XY_v,x_p_nu_p_rho     ,mu )
                            + adj(getLink(XY_v,x_p_mu           ,nu ))
                            * adj(getLink(U_v ,x_p_nu           ,mu ))
                            *     getLink(U_v ,x_p_nu           ,rho)
                    )*adj(getLink(U_v,x_p_rho,nu))*adj(getLink(U_v,x,rho));

                    res += (  adj(getLink(XY_v,x_p_mu           ,nu ))
                            * adj(getLink(U_v ,x_p_nu           ,mu ))
                            * adj(getLink(U_v ,x_p_nu_m_rho     ,rho))
                            + adj(getLink(U_v ,x_p_mu_m_rho     ,rho))
                            *     getLink(U_v ,x_p_mu_m_rho     ,nu )
                            *     getLink(XY_v,x_p_nu_m_rho     ,mu )
                    )*adj(getLink(U_v,x_m_rho,nu))*getLink(U_v,x_m_rho,rho);

                    res += (      getLink(U_v ,x_p_mu           ,nu )
                            *     getLink(XY_v,x_p_mu_p_nu_m_rho,rho) 
                            * adj(getLink(U_v ,x_p_mu_m_rho     ,nu ))
                            + adj(getLink(U_v ,x_p_mu_m_nu      ,nu ))
                            *     getLink(XY_v,x_p_mu_m_nu_m_rho,rho)
                            *     getLink(U_v ,x_p_mu_m_nu_m_rho,nu )
                    )*adj(getLink(U_v,x_m_rho,mu))*getLink(U_v,x_m_rho,rho);

                    res += (      getLink(XY_v,x_p_mu_m_nu      ,nu )
                            * adj(getLink(U_v ,x_m_nu           ,mu )) 
                            * adj(getLink(U_v ,x_m_nu_m_rho     ,rho))
                            + adj(getLink(U_v ,x_p_mu_m_rho     ,rho))
                            * adj(getLink(U_v ,x_p_mu_m_nu_m_rho,nu ))
                            *     getLink(XY_v,x_m_nu_m_rho     ,mu )
                    )*getLink(U_v,x_m_nu_m_rho,nu)*getLink(U_v,x_m_rho,rho);
}
                    setLink(F_v[x->_offset](mu), F_v(x->_offset)(mu) + c5*res);
                }
            }
        })
    }


    template<int term>
    void sevenLinkDeriv(GF& Fghost, GF& Ughost, GF& XYCghost, GeneralLocalStencil gStencil7, Real c7, int mu) const {
        
        autoView(U_v , Ughost , AcceleratorRead);
        autoView(XY_v, XYCghost, AcceleratorRead);
        autoView(F_v , Fghost , AcceleratorWrite);
        int Nsites = U_v.size();
        auto gStencil7_v = gStencil7.View(AcceleratorRead);
        typedef decltype(getLink(U_v,gStencil7.GetEntry(0,0),0)) U3matrix;

        accelerator_for(site,Nsites,Simd::Nsimd(),{
            U3matrix res, U1; 
            for(int nu=0;nu<Nd;nu++) {
                if(nu==mu) continue;
                for(int rho=0;rho<Nd;rho++) {
                    if (rho == mu || rho == nu) continue;
                    for(int sig=0;sig<Nd;sig++) {
                        if (sig == mu || sig == nu || sig == rho) continue;
                        int s = HISQStencilIndex(mu,nu,rho,sig,"7STAPLE");

                        auto [x                      , x_p_mu                 , x_p_mu_p_nu            , x_p_mu_p_nu_p_rho      ,
                              x_p_mu_p_nu_p_rho_m_sig, x_p_mu_p_nu_m_rho      , x_p_mu_p_nu_m_rho_m_sig, x_p_mu_p_nu_m_sig      ,
                              x_p_mu_m_nu            , x_p_mu_m_nu_p_rho      , x_p_mu_m_nu_p_rho_p_sig, x_p_mu_m_nu_p_rho_m_sig,
                              x_p_mu_m_nu_m_rho      , x_p_mu_m_nu_m_rho_p_sig, x_p_mu_m_nu_m_rho_m_sig, x_p_mu_m_nu_p_sig      ,
                              x_p_mu_m_nu_m_sig      , x_p_mu_p_rho           , x_p_mu_p_rho_p_sig     , x_p_mu_p_rho_m_sig     , 
                              x_p_mu_m_rho           , x_p_mu_m_rho_p_sig     , x_p_mu_m_rho_m_sig     , x_p_mu_p_sig           , 
                              x_p_mu_m_sig           , x_p_nu                 , x_p_nu_p_rho           , x_p_nu_p_rho_p_sig     , 
                              x_p_nu_p_rho_m_sig     , x_p_nu_m_rho           , x_p_nu_m_rho_p_sig     , x_p_nu_m_rho_m_sig     , 
                              x_p_rho                , x_p_rho_m_nu           , x_p_rho_p_sig          , x_p_rho_m_sig          , 
                              x_m_nu                 , x_m_nu_p_rho           , x_m_nu_p_rho_p_sig     , x_m_nu_p_rho_m_sig     , 
                              x_m_nu_m_rho           , x_m_nu_m_rho_p_sig     , x_m_nu_m_rho_m_sig     , x_m_rho                , 
                              x_m_rho_p_sig          , x_m_rho_m_sig
                        ] = get7StaplePoints(gStencil7_v,s,site);

if constexpr(term==0) {
                        res = adj(getLink(XY_v,x_p_mu,nu))*adj(getLink(U_v,x_p_nu,mu))
                              *(getLink(U_v,x_p_nu,rho)
                                *(   getLink(U_v,x_p_nu_p_rho      ,sig)
                                *adj(getLink(U_v,x_p_rho_p_sig     ,nu ))
                                *adj(getLink(U_v,x_p_rho           ,sig))
                                +adj(getLink(U_v,x_p_nu_p_rho_m_sig,sig))
                                *adj(getLink(U_v,x_p_rho_m_sig     ,nu ))
                                *    getLink(U_v,x_p_rho_m_sig     ,sig)
                                 )*adj(getLink(U_v,x,rho))

                                +adj(getLink(U_v,x_p_nu_m_rho,rho))
                                *(   getLink(U_v,x_p_nu_m_rho      ,sig)
                                *adj(getLink(U_v,x_m_rho_p_sig     ,nu ))
                                *adj(getLink(U_v,x_m_rho           ,sig))
                                +adj(getLink(U_v,x_p_nu_m_rho_m_sig,sig))
                                *adj(getLink(U_v,x_m_rho_m_sig     ,nu ))
                                *    getLink(U_v,x_m_rho_m_sig     ,sig)
                                 )*getLink(U_v,x_m_rho,rho)
                               ); 
} 
if constexpr(term==1) {
                        res = getLink(XY_v,x_p_mu_m_nu,nu)*adj(getLink(U_v,x_m_nu,mu))
                              *(getLink(U_v,x_m_nu,rho)
                                *(   getLink(U_v,x_m_nu_p_rho      ,sig)
                                *    getLink(U_v,x_m_nu_p_rho_p_sig,nu )
                                *adj(getLink(U_v,x_p_rho           ,sig))
                                +adj(getLink(U_v,x_m_nu_p_rho_m_sig,sig))
                                *    getLink(U_v,x_m_nu_p_rho_m_sig,nu ) 
                                *    getLink(U_v,x_p_rho_m_sig     ,sig)
                                 )*adj(getLink(U_v,x,rho))

                                +adj(getLink(U_v,x_m_nu_m_rho,rho))
                                *(   getLink(U_v,x_m_nu_m_rho      ,sig)
                                *    getLink(U_v,x_m_nu_m_rho_p_sig,nu ) 
                                *adj(getLink(U_v,x_m_rho           ,sig))
                                +adj(getLink(U_v,x_m_nu_m_rho_m_sig,sig))
                                *    getLink(U_v,x_m_nu_m_rho_m_sig,nu ) 
                                *    getLink(U_v,x_m_rho_m_sig     ,sig)
                                 )*getLink(U_v,x_m_rho,rho)
                               );
}
if constexpr(term==2) {
                        U1 = adj(getLink(U_v,x_p_nu,mu));

                        res = (     getLink(U_v ,x_p_mu           ,sig)
                               *adj(getLink(XY_v,x_p_mu_p_sig     ,nu ))
                               *adj(getLink(U_v ,x_p_mu_p_nu      ,sig))
                               +adj(getLink(U_v ,x_p_mu_m_sig     ,sig))
                               *adj(getLink(XY_v,x_p_mu_m_sig     ,nu ))
                               *    getLink(U_v ,x_p_mu_p_nu_m_sig,sig)
                              )*U1*getLink(U_v,x_p_nu,rho)*adj(getLink(U_v,x_p_rho,nu))*adj(getLink(U_v,x,rho))

                            + (     getLink(U_v ,x_p_mu           ,sig) 
                               *adj(getLink(XY_v,x_p_mu_p_sig     ,nu ))
                               *adj(getLink(U_v ,x_p_mu_p_nu      ,sig))
                               +adj(getLink(U_v ,x_p_mu_m_sig     ,sig))
                               *adj(getLink(XY_v,x_p_mu_m_sig     ,nu ))
                               *    getLink(U_v ,x_p_mu_p_nu_m_sig,sig)
                              )*U1*adj(getLink(U_v,x_p_nu_m_rho,rho))*adj(getLink(U_v,x_m_rho,nu))*getLink(U_v,x_m_rho,rho); 
}
if constexpr(term==3) {
                        U1 = adj(getLink(U_v,x_m_nu,mu));
                        
                        res = (     getLink(U_v ,x_p_mu            ,sig) 
                               *    getLink(XY_v,x_p_mu_m_nu_p_sig ,nu )
                               *adj(getLink(U_v ,x_p_mu_m_nu       ,sig))
                               +adj(getLink(U_v ,x_p_mu_m_sig      ,sig))
                               *    getLink(XY_v,x_p_mu_m_nu_m_sig ,nu )
                               *    getLink(U_v ,x_p_mu_m_nu_m_sig ,sig)
                              )*U1*getLink(U_v,x_m_nu,rho)*getLink(U_v,x_m_nu_p_rho,nu)*adj(getLink(U_v,x,rho))

                            + (     getLink(U_v ,x_p_mu            ,sig) 
                               *    getLink(XY_v,x_p_mu_m_nu_p_sig ,nu )
                               *adj(getLink(U_v ,x_p_mu_m_nu       ,sig))
                               +adj(getLink(U_v ,x_p_mu_m_sig      ,sig))
                               *    getLink(XY_v,x_p_mu_m_nu_m_sig ,nu )
                               *    getLink(U_v ,x_p_mu_m_nu_m_sig ,sig)
                              )*U1*adj(getLink(U_v,x_m_nu_m_rho,rho))*getLink(U_v,x_m_nu_m_rho,nu)*getLink(U_v,x_m_rho,rho);
}
if constexpr(term==4) {
                        res = (getLink(U_v,x_p_mu,rho)
                                *(     getLink(U_v ,x_p_mu_p_rho           ,sig)
                                  *adj(getLink(XY_v,x_p_mu_p_rho_p_sig     ,nu ))
                                  *adj(getLink(U_v ,x_p_mu_p_nu_p_rho      ,sig))
                                  +adj(getLink(U_v ,x_p_mu_p_rho_m_sig     ,sig))
                                  *adj(getLink(XY_v,x_p_mu_p_rho_m_sig     ,nu ))
                                  *    getLink(U_v ,x_p_mu_p_nu_p_rho_m_sig,sig) 
                                 )*adj(getLink(U_v,x_p_mu_p_nu,rho))

                                +adj(getLink(U_v,x_p_mu_m_rho,rho))
                                *(     getLink(U_v ,x_p_mu_m_rho           ,sig)
                                  *adj(getLink(XY_v,x_p_mu_m_rho_p_sig     ,nu ))
                                  *adj(getLink(U_v ,x_p_mu_p_nu_m_rho      ,sig))
                                  +adj(getLink(U_v ,x_p_mu_m_rho_m_sig     ,sig))
                                  *adj(getLink(XY_v,x_p_mu_m_rho_m_sig     ,nu ))
                                  *    getLink(U_v ,x_p_mu_p_nu_m_rho_m_sig,sig)
                                 )*getLink(U_v,x_p_mu_p_nu_m_rho,rho)
                              )*adj(getLink(U_v,x_p_nu,mu))*adj(getLink(U_v,x,nu));
}
if constexpr(term==5) {
                        res = (getLink(U_v,x_p_mu,rho)
                               *(     getLink(U_v ,x_p_mu_p_rho           ,sig)
                                 *    getLink(XY_v,x_p_mu_m_nu_p_rho_p_sig,nu )
                                 *adj(getLink(U_v ,x_p_mu_m_nu_p_rho      ,sig))
                                 +adj(getLink(U_v ,x_p_mu_p_rho_m_sig     ,sig))
                                 *    getLink(XY_v,x_p_mu_m_nu_p_rho_m_sig,nu )
                                 *    getLink(U_v ,x_p_mu_m_nu_p_rho_m_sig,sig) 
                                )*adj(getLink(U_v,x_p_mu_m_nu,rho))

                              + adj(getLink(U_v,x_p_mu_m_rho,rho))
                               *(     getLink(U_v ,x_p_mu_m_rho           ,sig)
                                 *    getLink(XY_v,x_p_mu_m_nu_m_rho_p_sig,nu )
                                 *adj(getLink(U_v ,x_p_mu_m_nu_m_rho      ,sig))
                                 +adj(getLink(U_v ,x_p_mu_m_rho_m_sig     ,sig))
                                 *    getLink(XY_v,x_p_mu_m_nu_m_rho_m_sig,nu )
                                 *    getLink(U_v ,x_p_mu_m_nu_m_rho_m_sig,sig)
                                )* getLink(U_v,x_p_mu_m_nu_m_rho,rho)
                              )*adj(getLink(U_v,x_m_nu,mu))*getLink(U_v,x_m_nu,nu);
}
if constexpr(term==6) {
                        res = getLink(U_v,x_p_mu,nu)
                              *(getLink(U_v,x_p_mu_p_nu,rho)
                                *(    getLink(U_v ,x_p_mu_p_nu_p_rho      ,sig)
                                 *    getLink(XY_v,x_p_nu_p_rho_p_sig     ,mu )
                                 *adj(getLink(U_v ,x_p_nu_p_rho           ,sig))
                                 +adj(getLink(U_v ,x_p_mu_p_nu_p_rho_m_sig,sig))
                                 *    getLink(XY_v,x_p_nu_p_rho_m_sig     ,mu )
                                 *    getLink(U_v ,x_p_nu_p_rho_m_sig     ,sig)
                                )*adj(getLink(U_v ,x_p_nu                 ,rho))

                                +adj(getLink(U_v,x_p_mu_p_nu_m_rho,rho))
                                *(adj(getLink(U_v ,x_p_mu_p_nu_m_rho_m_sig,sig))
                                 *    getLink(XY_v,x_p_nu_m_rho_m_sig     ,mu )
                                 *    getLink(U_v ,x_p_nu_m_rho_m_sig     ,sig)
                                 +    getLink(U_v ,x_p_mu_p_nu_m_rho      ,sig)
                                 *    getLink(XY_v,x_p_nu_m_rho_p_sig     ,mu )
                                 *adj(getLink(U_v ,x_p_nu_m_rho           ,sig))
                                )*getLink(U_v,x_p_nu_m_rho,rho)
                              )*adj(getLink(U_v,x,nu));
}
if constexpr(term==7) {
                        res = adj(getLink(U_v,x_p_mu_m_nu,nu))
                              *(getLink(U_v,x_p_mu_m_nu,rho)
                                *(    getLink(U_v ,x_p_mu_m_nu_p_rho      ,sig)
                                 *    getLink(XY_v,x_m_nu_p_rho_p_sig     ,mu )
                                 *adj(getLink(U_v ,x_m_nu_p_rho           ,sig))
                                 +adj(getLink(U_v ,x_p_mu_m_nu_p_rho_m_sig,sig))
                                 *    getLink(XY_v,x_m_nu_p_rho_m_sig     ,mu )
                                 *    getLink(U_v ,x_m_nu_p_rho_m_sig     ,sig)
                                )*adj(getLink(U_v,x_m_nu,rho))

                                +adj(getLink(U_v,x_p_mu_m_nu_m_rho,rho))
                                *(    getLink(U_v ,x_p_mu_m_nu_m_rho      ,sig)
                                 *    getLink(XY_v,x_m_nu_m_rho_p_sig     ,mu )
                                 *adj(getLink(U_v ,x_m_nu_m_rho           ,sig))
                                 +adj(getLink(U_v ,x_p_mu_m_nu_m_rho_m_sig,sig))
                                 *    getLink(XY_v,x_m_nu_m_rho_m_sig     ,mu )
                                 *    getLink(U_v ,x_m_nu_m_rho_m_sig     ,sig)
                                )*getLink(U_v,x_m_nu_m_rho,rho)
                              )*getLink(U_v,x_m_nu,nu);
}
if constexpr(term==8) {
                        res = getLink(U_v,x_p_mu,nu)*adj(getLink(U_v,x_p_nu,mu))
                              *(getLink(U_v,x_p_nu,rho)
                                *(    getLink(U_v ,x_p_nu_p_rho           ,sig)
                                 *    getLink(XY_v,x_p_rho_p_sig          ,nu )
                                 *adj(getLink(U_v ,x_p_rho                ,sig))
                                 +adj(getLink(U_v ,x_p_nu_p_rho_m_sig     ,sig))
                                 *    getLink(XY_v,x_p_rho_m_sig          ,nu )
                                 *    getLink(U_v ,x_p_rho_m_sig          ,sig)
                                )*adj(getLink(U_v,x,rho))

                                +adj(getLink(U_v,x_p_nu_m_rho,rho))
                                *(    getLink(U_v ,x_p_nu_m_rho           ,sig)
                                 *    getLink(XY_v,x_m_rho_p_sig          ,nu )
                                 *adj(getLink(U_v ,x_m_rho                ,sig))
                                 +adj(getLink(U_v ,x_p_nu_m_rho_m_sig     ,sig))
                                 *    getLink(XY_v,x_m_rho_m_sig          ,nu )
                                 *    getLink(U_v ,x_m_rho_m_sig          ,sig)
                                )*getLink(U_v,x_m_rho,rho)
                              );
}
if constexpr(term==9) {
                        res = adj(getLink(U_v,x_p_mu_m_nu,nu))*adj(getLink(U_v,x_m_nu,mu))
                              *(getLink(U_v,x_m_nu,rho)
                                *(    getLink(U_v ,x_m_nu_p_rho           ,sig)
                                 *adj(getLink(XY_v,x_m_nu_p_rho_p_sig     ,nu ))
                                 *adj(getLink(U_v ,x_p_rho                ,sig))
                                 +adj(getLink(U_v ,x_m_nu_p_rho_m_sig     ,sig))
                                 *adj(getLink(XY_v,x_m_nu_p_rho_m_sig     ,nu ))
                                 *    getLink(U_v ,x_p_rho_m_sig          ,sig)
                                )*adj(getLink(U_v,x,rho))

                                +adj(getLink(U_v,x_m_nu_m_rho,rho))
                                *(    getLink(U_v ,x_m_nu_m_rho           ,sig)
                                 *adj(getLink(XY_v,x_m_nu_m_rho_p_sig     ,nu ))
                                 *adj(getLink(U_v ,x_m_rho                ,sig))
                                 +adj(getLink(U_v ,x_m_nu_m_rho_m_sig     ,sig))
                                 *adj(getLink(XY_v,x_m_nu_m_rho_m_sig     ,nu ))
                                 *    getLink(U_v ,x_m_rho_m_sig          ,sig)
                                )*getLink(U_v,x_m_rho,rho)
                              );
}
if constexpr(term==10) {
                        U1 = adj(getLink(U_v,x_p_nu,mu));
                        
                        res = (     getLink(U_v,x_p_mu            ,sig) 
                               *    getLink(U_v,x_p_mu_p_sig      ,nu )
                               *adj(getLink(U_v,x_p_mu_p_nu       ,sig))
                               +adj(getLink(U_v,x_p_mu_m_sig      ,sig)) 
                               *    getLink(U_v,x_p_mu_m_sig      ,nu )
                               *    getLink(U_v,x_p_mu_p_nu_m_sig ,sig)
                              )*U1*getLink(U_v,x_p_nu,rho)*getLink(XY_v,x_p_rho,nu)*adj(getLink(U_v,x,rho))

                              +(adj(getLink(U_v,x_p_mu_m_sig      ,sig)) 
                               *    getLink(U_v,x_p_mu_m_sig      ,nu )
                               *    getLink(U_v,x_p_mu_p_nu_m_sig ,sig)
                               +    getLink(U_v,x_p_mu            ,sig)
                               *    getLink(U_v,x_p_mu_p_sig      ,nu )
                               *adj(getLink(U_v,x_p_mu_p_nu       ,sig))
                              )*U1*adj(getLink(U_v,x_p_nu_m_rho,rho))*getLink(XY_v,x_m_rho,nu)*getLink(U_v,x_m_rho,rho);
}
if constexpr(term==11) {
                        U1 = adj(getLink(U_v,x_m_nu,mu));
                        
                        res = (     getLink(U_v,x_p_mu            ,sig) 
                               *adj(getLink(U_v,x_p_mu_m_nu_p_sig ,nu ))
                               *adj(getLink(U_v,x_p_mu_m_nu       ,sig))
                               +adj(getLink(U_v,x_p_mu_m_sig      ,sig))
                               *adj(getLink(U_v,x_p_mu_m_nu_m_sig ,nu ))
                               *    getLink(U_v,x_p_mu_m_nu_m_sig ,sig)
                              )*U1*getLink(U_v,x_m_nu,rho)*adj(getLink(XY_v,x_m_nu_p_rho,nu))*adj(getLink(U_v,x,rho))

                             +(     getLink(U_v,x_p_mu            ,sig)
                               *adj(getLink(U_v,x_p_mu_m_nu_p_sig ,nu ))
                               *adj(getLink(U_v,x_p_mu_m_nu       ,sig))
                               +adj(getLink(U_v,x_p_mu_m_sig      ,sig))
                               *adj(getLink(U_v,x_p_mu_m_nu_m_sig ,nu ))
                               *    getLink(U_v,x_p_mu_m_nu_m_sig ,sig)
                              )*U1*adj(getLink(U_v,x_m_nu_m_rho,rho))*adj(getLink(XY_v,x_m_nu_m_rho,nu))*getLink(U_v,x_m_rho,rho);
}
if constexpr(term==12) {
                        res = (getLink(U_v,x_p_mu,rho)
                               *(    getLink(U_v,x_p_mu_p_rho           ,sig)
                                *    getLink(U_v,x_p_mu_p_rho_p_sig     ,nu ) 
                                *adj(getLink(U_v,x_p_mu_p_nu_p_rho      ,sig))
                                +adj(getLink(U_v,x_p_mu_p_rho_m_sig     ,sig))
                                *    getLink(U_v,x_p_mu_p_rho_m_sig     ,nu )
                                *    getLink(U_v,x_p_mu_p_nu_p_rho_m_sig,sig) 
                               )*adj(getLink(U_v,x_p_mu_p_nu,rho))

                               +adj(getLink(U_v,x_p_mu_m_rho,rho))
                               *(    getLink(U_v,x_p_mu_m_rho           ,sig)
                                *    getLink(U_v,x_p_mu_m_rho_p_sig     ,nu )
                                *adj(getLink(U_v,x_p_mu_p_nu_m_rho      ,sig))
                                +adj(getLink(U_v,x_p_mu_m_rho_m_sig     ,sig))
                                *    getLink(U_v,x_p_mu_m_rho_m_sig     ,nu )
                                *    getLink(U_v,x_p_mu_p_nu_m_rho_m_sig,sig)
                               )*getLink(U_v,x_p_mu_p_nu_m_rho,rho)
                             )*adj(getLink(U_v,x_p_nu,mu))*getLink(XY_v,x,nu);
}
if constexpr(term==13) {
                        res = (adj(getLink(U_v,x_p_mu_m_rho,rho))
                               *(    getLink(U_v,x_p_mu_m_rho           ,sig)
                                *adj(getLink(U_v,x_p_mu_m_nu_m_rho_p_sig,nu ))
                                *adj(getLink(U_v,x_p_mu_m_nu_m_rho      ,sig))
                                +adj(getLink(U_v,x_p_mu_m_rho_m_sig     ,sig))
                                *adj(getLink(U_v,x_p_mu_m_nu_m_rho_m_sig,nu ))
                                *    getLink(U_v,x_p_mu_m_nu_m_rho_m_sig,sig) 
                               )*getLink(U_v,x_p_mu_m_nu_m_rho,rho)

                               +getLink(U_v,x_p_mu,rho)
                               *(adj(getLink(U_v,x_p_mu_p_rho_m_sig     ,sig))
                                *adj(getLink(U_v,x_p_mu_m_nu_p_rho_m_sig,nu ))
                                *    getLink(U_v,x_p_mu_m_nu_p_rho_m_sig,sig)
                                +    getLink(U_v,x_p_mu_p_rho           ,sig)
                                *adj(getLink(U_v,x_p_mu_m_nu_p_rho_p_sig,nu ))
                                *adj(getLink(U_v,x_p_mu_m_nu_p_rho      ,sig))
                               )*adj(getLink(U_v,x_p_mu_m_nu,rho))
                             )*adj(getLink(U_v,x_m_nu,mu))*adj(getLink(XY_v,x_m_nu,nu));
} 
                        setLink(F_v[x->_offset](mu), F_v(x->_offset)(mu) + c7*res);
                    }
                }
            }
        })
    
    }


    GF naikLinkDeriv(std::vector<Real> vecdt, std::vector<FF>& vecx, std::vector<int> n_orders_naik, int n_naiks, Real cnaik) {
        auto grid = this->_grid;
        GF temp(grid);
        std::vector<LF> Wv(Nd, grid);
        std::vector<LF> XYdag(Nd, grid);
        std::vector<LF> ddW(Nd, grid);
        temp = outerProductHISQ(vecx, vecdt, n_orders_naik, n_naiks, 3);
        for (int mu = 0; mu < Nd; mu++) {
            Wv[mu]    = PeekIndex<LorentzIndex>(_Wmu, mu);
            XYdag[mu] = PeekIndex<LorentzIndex>(temp, mu);
            ddW[mu]   = Zero();
        }
        for (int mu = 0; mu < Nd; mu++) {
            ddW[mu] =    Cshift(   Wv[mu],mu, 1)*Cshift(   Wv[mu],mu, 2)*       XYdag[mu]
                       + Cshift(   Wv[mu],mu, 1)*Cshift(XYdag[mu],mu,-1)*Cshift(   Wv[mu],mu,-1)
                       + Cshift(XYdag[mu],mu,-2)*Cshift(   Wv[mu],mu,-2)*Cshift(   Wv[mu],mu,-1); 
        }
        for (int mu = 0; mu < Nd; mu++) {
            PokeIndex<LorentzIndex>(temp, ddW[mu], mu);
        }
        return cnaik*temp;
    } 

    
    GF lepageLinkDeriv(GF& XY, Real clp) {

        auto grid = this->_grid;
        GF temp(grid);
        std::vector<LF> Wv(Nd, grid);
        std::vector<LF> XYdag(Nd, grid);
        std::vector<LF> ddW(Nd, grid);

        for (int mu = 0; mu < Nd; mu++) {
            Wv[mu]    = PeekIndex<LorentzIndex>(_Wmu, mu);
            ddW[mu]   = Zero();
            XYdag[mu] = adj(PeekIndex<LorentzIndex>(XY, mu));
        }

        for (int mu = 0; mu < Nd; mu++) 
        for (int nu = 0; nu < Nd; nu++) {
            if(mu==nu) continue;

            // (forward)
            ddW[mu] = ddW[mu] + adj(Gimpl::CovShiftBackward(Wv[mu],mu,
                                            Gimpl::CovShiftForward(Wv[nu],nu,
                                              Gimpl::CovShiftForward(Wv[mu],mu,
                                                Gimpl::CovShiftForward(Wv[mu],mu,
                                                  Gimpl::CovShiftIdentityBackward(XYdag[nu],nu))))));
            // (backward)
            ddW[mu] = ddW[mu] + adj(Gimpl::CovShiftBackward(Wv[mu],mu,
                                            Gimpl::CovShiftBackward(Wv[nu],nu,
                                              Gimpl::CovShiftForward(Wv[mu],mu,
                                                Gimpl::CovShiftForward(Wv[mu],mu,
                                                  Gimpl::CovShiftIdentityForward(XYdag[nu],nu))))));
            // (forward)
            ddW[mu] = ddW[mu] + adj(Gimpl::CovShiftForward(Wv[nu],nu,
                                            Gimpl::CovShiftForward(Wv[mu],mu,
                                              Gimpl::CovShiftForward(Wv[mu],mu,
                                                Gimpl::CovShiftBackward(XYdag[nu],nu,
                                                  Gimpl::CovShiftIdentityBackward(Wv[mu],mu))))));
            // (backward)
            ddW[mu] = ddW[mu] + adj(Gimpl::CovShiftBackward(Wv[nu],nu,
                                            Gimpl::CovShiftForward(Wv[mu],mu,
                                              Gimpl::CovShiftForward(Wv[mu],mu,
                                                Gimpl::CovShiftForward(XYdag[nu],nu,
                                                  Gimpl::CovShiftIdentityBackward(Wv[mu],mu))))));
            // (forward)
            ddW[mu] = ddW[mu] + adj(Gimpl::CovShiftForward(Wv[nu],nu,
                                            Gimpl::CovShiftForward(Wv[nu],nu,
                                              Gimpl::CovShiftForward(XYdag[mu],mu,
                                                Gimpl::CovShiftBackward(Wv[nu],nu,
                                                  Gimpl::CovShiftIdentityBackward(Wv[nu],nu))))));
            // (backward)
            ddW[mu] = ddW[mu] + adj(Gimpl::CovShiftBackward(Wv[nu],nu,
                                            Gimpl::CovShiftBackward(Wv[nu],nu,
                                              Gimpl::CovShiftForward(XYdag[mu],mu,
                                                Gimpl::CovShiftForward(Wv[nu],nu,
                                                  Gimpl::CovShiftIdentityForward(Wv[nu],nu))))));
            // (forward)
            ddW[mu] = ddW[mu] + adj(Gimpl::CovShiftBackward(Wv[mu],mu,
                                            Gimpl::CovShiftForward(XYdag[nu],nu,
                                              Gimpl::CovShiftForward(Wv[mu],mu,
                                                Gimpl::CovShiftForward(Wv[mu],mu,
                                                  Gimpl::CovShiftIdentityBackward(Wv[nu],nu))))));
            // (backward)
            ddW[mu] = ddW[mu] + adj(Gimpl::CovShiftBackward(Wv[mu],mu,
                                            Gimpl::CovShiftBackward(XYdag[nu],nu,
                                              Gimpl::CovShiftForward(Wv[mu],mu,
                                                Gimpl::CovShiftForward(Wv[mu],mu,
                                                  Gimpl::CovShiftIdentityForward(Wv[nu],nu))))));
            // (forward)
            ddW[mu] = ddW[mu] + adj(Gimpl::CovShiftForward(XYdag[nu],nu,
                                            Gimpl::CovShiftForward(Wv[mu],mu,
                                              Gimpl::CovShiftForward(Wv[mu],mu,
                                                Gimpl::CovShiftBackward(Wv[nu],nu,
                                                  Gimpl::CovShiftIdentityBackward(Wv[mu],mu))))));
            // (backward)
            ddW[mu] = ddW[mu] + adj(Gimpl::CovShiftBackward(XYdag[nu],nu,
                                            Gimpl::CovShiftForward(Wv[mu],mu,
                                              Gimpl::CovShiftForward(Wv[mu],mu,
                                                Gimpl::CovShiftForward(Wv[nu],nu,
                                                  Gimpl::CovShiftIdentityBackward(Wv[mu],mu))))));
        }

        for (int mu = 0; mu < Nd; mu++) {
            PokeIndex<LorentzIndex>(temp, ddW[mu], mu);
        }

        return clp*temp;
    }


    // We are calculating the force using the rational approximation. The goal is that we can approximate
    // (Mdag M)^(-nf/4) = alpha_0 + sum_l alpha_l/(M^dag M + beta_l). Hence the index l runs over the
    // introduce a different "Naik epsilon" for each M. Hence, we can think of the total application
    // of this operator as having an index inaik, running over the different Naik epsilons; for each inaik
    // there is a possibly different order_inaik, then the operator has an index l running up to order_inaik.
    // All terms with inaik=0 correspond to epsilon_Naik = 0.
    //
    // Intent: OUT--u_force
    //          IN--vecdt: Monte Carlo separation vector times alpha_{inaik,0}. 
    //              vecx: A vector of fermion fields coming from the MILC code. It is organized so that 
    //                    |X_l> = (Mdag M + beta_l)^-1 |Phi> is on even sites, |Y_l>=D|X_l> is on odd sites.
    //                    All the |X_l> for i=0 come first in memory, followed by all the |X_l> with
    //                    i=1 in memory, and so on.
    //              n_orders_naik: Indexed by unique naik epsilon.
    void force(GF& u_force, std::vector<Real> vecdt, std::vector<FF>& vecx, std::vector<int> n_orders_naik) {

        HISQParameters<Real> hp = this->_linkParams;
        auto grid = this->_grid;

        GF XY(grid);    // outer product field
        GF temp(grid);  // used to accumulate N-link force contributions and projU3Deriv 

        u_force = Zero();

        if(hp.asqtad_cnaik!=0) u_force += naikLinkDeriv(vecdt, vecx, n_orders_naik, hp.n_naiks, hp.asqtad_cnaik); 

        XY = outerProductHISQ(vecx, vecdt, n_orders_naik, hp.n_naiks, 1);   
        u_force += hp.asqtad_c1*XY;

        if(hp.asqtad_clp!=0) u_force += lepageLinkDeriv(XY, hp.asqtad_clp); 


        // ---------------------------------- N-LINK DERIVATIVES (ASQTAD)

        PaddedCell Ghost(_HaloDepth,grid);

        temp = Zero();

        GF UWghost  = Ghost.Exchange(_Wmu);    // Plays role of U or W 
        GF XYCghost = Ghost.Exchange(XY);      // Plays role of XY or chain rule = dW/dV*dX/dW
        GF Fghost   = Ghost.Exchange(temp);

        std::vector<Coordinate> shifts3 = createHISQStencil("3STAPLE");
        std::vector<Coordinate> shifts5 = createHISQStencil("5STAPLE");
        std::vector<Coordinate> shifts7 = createHISQStencil("7STAPLE");
        GeneralLocalStencil gStencil3(UWghost.Grid(),shifts3);
        GeneralLocalStencil gStencil5(UWghost.Grid(),shifts5);
        GeneralLocalStencil gStencil7(UWghost.Grid(),shifts7);

        for(int mu=0;mu<Nd;mu++) {
            if(hp.asqtad_c3!=0) { 
                threeLinkDeriv(    Fghost, UWghost, XYCghost, gStencil3, hp.asqtad_c3, mu);
            }
            if(hp.asqtad_c5!=0) {
                fiveLinkDeriv<0>(  Fghost, UWghost, XYCghost, gStencil5, hp.asqtad_c5, mu); 
                fiveLinkDeriv<1>(  Fghost, UWghost, XYCghost, gStencil5, hp.asqtad_c5, mu);
            }
            if(hp.asqtad_c7!=0) {
                sevenLinkDeriv<0>( Fghost, UWghost, XYCghost, gStencil7, hp.asqtad_c7, mu); 
                sevenLinkDeriv<1>( Fghost, UWghost, XYCghost, gStencil7, hp.asqtad_c7, mu); 
                sevenLinkDeriv<2>( Fghost, UWghost, XYCghost, gStencil7, hp.asqtad_c7, mu); 
                sevenLinkDeriv<3>( Fghost, UWghost, XYCghost, gStencil7, hp.asqtad_c7, mu); 
                sevenLinkDeriv<4>( Fghost, UWghost, XYCghost, gStencil7, hp.asqtad_c7, mu); 
                sevenLinkDeriv<5>( Fghost, UWghost, XYCghost, gStencil7, hp.asqtad_c7, mu); 
                sevenLinkDeriv<6>( Fghost, UWghost, XYCghost, gStencil7, hp.asqtad_c7, mu); 
                sevenLinkDeriv<7>( Fghost, UWghost, XYCghost, gStencil7, hp.asqtad_c7, mu); 
                sevenLinkDeriv<8>( Fghost, UWghost, XYCghost, gStencil7, hp.asqtad_c7, mu); 
                sevenLinkDeriv<9>( Fghost, UWghost, XYCghost, gStencil7, hp.asqtad_c7, mu); 
                sevenLinkDeriv<10>(Fghost, UWghost, XYCghost, gStencil7, hp.asqtad_c7, mu); 
                sevenLinkDeriv<11>(Fghost, UWghost, XYCghost, gStencil7, hp.asqtad_c7, mu); 
                sevenLinkDeriv<12>(Fghost, UWghost, XYCghost, gStencil7, hp.asqtad_c7, mu); 
                sevenLinkDeriv<13>(Fghost, UWghost, XYCghost, gStencil7, hp.asqtad_c7, mu);
            } 
        }

        u_force += Ghost.Extract(Fghost);

        // ------------------------------------------- U3 PROJ DERIVATIVE 

        projU3Deriv(temp, _Vmu, u_force, 5e-5);
        u_force = hp.fat7_c1*temp;

        // ------------------------------------ N-LINK DERIVATIVES (FAT7) 

        UWghost  = Ghost.Exchange(_Umu);
        XYCghost = Ghost.Exchange(temp); 
        Fghost   = Zero(); 

        for(int mu=0;mu<Nd;mu++) {
            if(hp.fat7_c3!=0) { 
                threeLinkDeriv(    Fghost, UWghost, XYCghost, gStencil3, hp.fat7_c3, mu);
            }
            if(hp.fat7_c5!=0) {
                fiveLinkDeriv<0>(  Fghost, UWghost, XYCghost, gStencil5, hp.fat7_c5, mu); 
                fiveLinkDeriv<1>(  Fghost, UWghost, XYCghost, gStencil5, hp.fat7_c5, mu);
            }
            if(hp.fat7_c7!=0) {
                sevenLinkDeriv<0>( Fghost, UWghost, XYCghost, gStencil7, hp.fat7_c7, mu); 
                sevenLinkDeriv<1>( Fghost, UWghost, XYCghost, gStencil7, hp.fat7_c7, mu); 
                sevenLinkDeriv<2>( Fghost, UWghost, XYCghost, gStencil7, hp.fat7_c7, mu); 
                sevenLinkDeriv<3>( Fghost, UWghost, XYCghost, gStencil7, hp.fat7_c7, mu); 
                sevenLinkDeriv<4>( Fghost, UWghost, XYCghost, gStencil7, hp.fat7_c7, mu); 
                sevenLinkDeriv<5>( Fghost, UWghost, XYCghost, gStencil7, hp.fat7_c7, mu); 
                sevenLinkDeriv<6>( Fghost, UWghost, XYCghost, gStencil7, hp.fat7_c7, mu); 
                sevenLinkDeriv<7>( Fghost, UWghost, XYCghost, gStencil7, hp.fat7_c7, mu); 
                sevenLinkDeriv<8>( Fghost, UWghost, XYCghost, gStencil7, hp.fat7_c7, mu); 
                sevenLinkDeriv<9>( Fghost, UWghost, XYCghost, gStencil7, hp.fat7_c7, mu); 
                sevenLinkDeriv<10>(Fghost, UWghost, XYCghost, gStencil7, hp.fat7_c7, mu); 
                sevenLinkDeriv<11>(Fghost, UWghost, XYCghost, gStencil7, hp.fat7_c7, mu); 
                sevenLinkDeriv<12>(Fghost, UWghost, XYCghost, gStencil7, hp.fat7_c7, mu); 
                sevenLinkDeriv<13>(Fghost, UWghost, XYCghost, gStencil7, hp.fat7_c7, mu);
            } 
        } // end mu loop

        u_force += Ghost.Extract(Fghost);

        // Close the loop: Multiply on the left by U_mu(x)
        LF force_mu(grid);
        for (int mu = 0; mu < Nd; mu++) {
            force_mu = PeekIndex<LorentzIndex>(_Umu, mu) * PeekIndex<LorentzIndex>(u_force, mu);
            PokeIndex<LorentzIndex>(u_force, force_mu, mu);
        }
    }

};


NAMESPACE_END(Grid);
