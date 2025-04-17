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
int HISQStencilIndex(int mu, int nu) {
    // Nshifts depends on how you built the stencil
    int Nshifts = 5;
    return Nshifts*nu + Nd*Nshifts*mu;
}


/*! @brief Create the mu-nu plane stencil. We allow mu==nu to make indexing the 
    stencil easier, but these entries will not be used. */
inline
std::vector<Coordinate> createHISQStencil() {
    std::vector<Coordinate> shifts;
    for(int mu=0;mu<Nd;mu++)
    for(int nu=0;nu<Nd;nu++) {
        appendShift<Nd>(shifts,mu);
        appendShift<Nd>(shifts,nu);
        appendShift<Nd>(shifts,shiftSignal::NO_SHIFT);
        appendShift<Nd>(shifts,mu,Back(nu));
        appendShift<Nd>(shifts,Back(nu));
    }
    return shifts;
}


/*! @brief Retrieve the stencil entries. */
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
    //          IN--U_v (thin links)
    //              gStencil (HISQ stencil)
    //              mu
    template<class linkRead, class linkWrite, class stencilRead>
    void threeLinkStaple(linkWrite U_fat_v, linkWrite U_3link_v, linkRead U_v, stencilRead gStencil_v, int mu) const {

        SmearingParameters<RealScalar> lt = this->_linkTreatment;
        typedef decltype(getLink(U_v,gStencil_v.GetEntry(0,0),0)) U3matrix;
        int Nsites = U_v.size();

        accelerator_for(site,Nsites,Simd::Nsimd(),{
            U3matrix U0, U1, U2, U3, U4, U5, res;
            for(int nu=0;nu<Nd;nu++) {
                if(nu==mu) continue;
                int s = HISQStencilIndex(mu,nu);

                auto [x_p_mu, x_p_nu, x, x_p_mu_m_nu, x_m_nu] = get3StaplePoints(gStencil_v,s,site);      

                // When you're deciding whether to take an adjoint, the question is: how is the
                // stored link oriented compared to the one you want? If I imagine myself traveling
                // with the to-be-updated link, I have two possible, alternative 3-link paths I can
                // take, one starting by going to the left, the other starting by going to the right.
                U0 = getLink(U_v,x_p_mu     ,nu);
                U1 = getLink(U_v,x_p_nu     ,mu);
                U2 = getLink(U_v,x          ,nu);
                U3 = getLink(U_v,x_p_mu_m_nu,nu);
                U4 = getLink(U_v,x_m_nu     ,mu);
                U5 = getLink(U_v,x_m_nu     ,nu);

                //  "left"          "right"
                res = U2*U1*adj(U0) + adj(U5)*U4*U3;

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


    // Intent: OUT--U_5link (sum of left and right 5-staples attached to U) 
    //              U_fat (accmulates the fat smearing)
    //          IN--U_v (thin links)
    //              U_3link (sum of left and right 3-staples attached to U)
    //              gStencil (HISQ stencil)
    //              mu
    //              updateFatLinks (in the force, you only want U_5link_v)
    template<class linkRead, class linkWrite, class stencilRead>
    void fiveLinkStaple(linkWrite U_fat_v, linkWrite U_5linkA_v, linkWrite U_5linkB_v, linkWrite U_3link_v, 
                        linkRead U_v, stencilRead gStencil_v, int mu) const {

        SmearingParameters<RealScalar> lt = this->_linkTreatment;
        typedef decltype(getLink(U_v,gStencil_v.GetEntry(0,0),0)) U3matrix;
        int Nsites = U_v.size();

        accelerator_for(site,Nsites,Simd::Nsimd(),{
            U3matrix U0, U1, U2, U3, U4, U5, res;
            int sigmaIndex = 0;
            for(int nu=0;nu<Nd;nu++) {
                if(nu==mu) continue;
                int s = HISQStencilIndex(mu,nu);
                for(int rho=0;rho<Nd;rho++) {
                    if (rho == mu || rho == nu) continue;

                    auto [x_p_mu, x_p_nu, x, x_p_mu_m_nu, x_m_nu] = get3StaplePoints(gStencil_v,s,site);      

                    U0 = getLink(      U_v,x_p_mu     ,nu );
                    U1 = getLink(U_3link_v,x_p_nu     ,rho);
                    U2 = getLink(      U_v,x          ,nu );
                    U3 = getLink(      U_v,x_p_mu_m_nu,nu );
                    U4 = getLink(U_3link_v,x_m_nu     ,rho);
                    U5 = getLink(      U_v,x_m_nu     ,nu );

                    res  = U2*U1*adj(U0) + adj(U5)*U4*U3;

                    // Counting 3-link staples: there are three planes attached to the to-be-updated link,
                    // which corresponds to three (forward+backward) staples. For the 5-link staples, for
                    // each plane, there are two remaining directions, so that there are six 5-link staples
                    // altogether. That will not fit in a single GaugeField object, so we use two. You can
                    // think of sigmaIndex and rho together as being the labels that pick out a particular
                    // 5-link staple. They therefore should not be interpreted as directions.
                    if(sigmaIndex<3) {
                        setLink(U_5linkA_v[x->_offset](rho), res);
                    } else {
                        setLink(U_5linkB_v[x->_offset](rho), res);
                    }    

                    setLink(U_fat_v[x->_offset](mu), U_fat_v(x->_offset)(mu) + lt.c_5*res);
                    sigmaIndex++;
                }
            }
        })
        return;
    }


    // Intent: OUT--U_fat (accmulates the fat smearing)
    //          IN--U_v (thin links)
    //              gStencil (HISQ stencil)
    //              mu
    template<class linkRead, class linkWrite, class stencilRead>
    void sevenLinkStaple(linkWrite U_fat_v, linkWrite U_5linkA_v, linkWrite U_5linkB_v, 
                         linkRead U_v, stencilRead gStencil_v, int mu) const {

        SmearingParameters<RealScalar> lt = this->_linkTreatment;
        typedef decltype(getLink(U_v,gStencil_v.GetEntry(0,0),0)) U3matrix;
        int Nsites = U_v.size();

        accelerator_for(site,Nsites,Simd::Nsimd(),{ 
            U3matrix U0, U1, U2, U3, U4, U5, res;
            int sigmaIndex = 0;
            for(int nu=0;nu<Nd;nu++) {
                if(nu==mu) continue;
                int s = HISQStencilIndex(mu,nu);
                for(int rho=0;rho<Nd;rho++) {
                    if (rho == mu || rho == nu) continue;
    
                    auto [x_p_mu, x_p_nu, x, x_p_mu_m_nu, x_m_nu] = get3StaplePoints(gStencil_v,s,site);      
    
                    U0 = getLink(U_v,x_p_mu,nu);
                    if(sigmaIndex<3) {
                        U1 = getLink(U_5linkB_v,x_p_nu,rho);
                    } else {
                        U1 = getLink(U_5linkA_v,x_p_nu,rho);
                    }  
                    U2 = getLink(U_v,x,nu);
                    U3 = getLink(U_v,x_p_mu_m_nu,nu);
                    if(sigmaIndex<3) {
                        U4 = getLink(U_5linkB_v,x_m_nu,rho);
                    } else {
                        U4 = getLink(U_5linkA_v,x_m_nu,rho);
                    }  
                    U5 = getLink(U_v,x_m_nu,nu);
    
                    res  = U2*U1*adj(U0) + adj(U5)*U4*U3;
    
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
        GF Ughost_fat(Ughost.Grid());
        GF Ughost_3link(Ughost.Grid());
        GF Ughost_5linkA(Ughost.Grid());
        GF Ughost_5linkB(Ughost.Grid());

        // mu-nu plane stencil.
        std::vector<Coordinate> shifts = createHISQStencil();

        // A GeneralLocalStencil has two indices: a site and stencil index 
        GeneralLocalStencil gStencil(Ughost.Grid(),shifts);

        // This is where contributions from the smearing get added together
        Ughost_fat=Zero();

        // This loop handles 3-, 5-, and 7-link constructs, minus Lepage and Naik.
        for(int mu=0;mu<Nd;mu++) {

            Ughost_3link =Zero(); Ughost_5linkA=Zero(); Ughost_5linkB=Zero();

            // Create the accessors
            autoView(U_v       , Ughost       , AcceleratorRead);
            autoView(U_fat_v   , Ughost_fat   , AcceleratorWrite);
            autoView(U_3link_v , Ughost_3link , AcceleratorWrite);
            autoView(U_5linkA_v, Ughost_5linkA, AcceleratorWrite);
            autoView(U_5linkB_v, Ughost_5linkB, AcceleratorWrite);

            auto gStencil_v = gStencil.View(AcceleratorRead); 

            typedef decltype(getLink(U_v,gStencil.GetEntry(0,0),0)) U3matrix;
            typedef decltype(U_v)                                   linkRead;
            typedef decltype(U_fat_v)                               linkWrite;
            typedef decltype(gStencil_v)                            stencilRead;

            // CODE IN A MORE CAREFUL TEST FOR THIS
            if((lt.c_7!=0) || (lt.c_5!=0) || (lt.c_3!=0))
                threeLinkStaple<linkRead,linkWrite,stencilRead>(U_fat_v,                         U_3link_v, U_v, gStencil_v, mu);
            if((lt.c_7!=0) || (lt.c_5!=0)) 
                fiveLinkStaple< linkRead,linkWrite,stencilRead>(U_fat_v, U_5linkA_v, U_5linkB_v, U_3link_v, U_v, gStencil_v, mu);
            if( lt.c_7!=0 ) 
                sevenLinkStaple<linkRead,linkWrite,stencilRead>(U_fat_v, U_5linkA_v, U_5linkB_v,            U_v, gStencil_v, mu);
        }

        // c1, c3, c5, c7 construct contributions
        u_smr = Ghost.Extract(Ughost_fat) + lt.c_1*u_thin;

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

            // Naik
            if(lt.c_naik!=0) {
                Vnaik[mu] = lt.c_naik*Gimpl::CovShiftForward(U[mu],mu,
                                        Gimpl::CovShiftForward(U[mu],mu,
                                          Gimpl::CovShiftIdentityForward(U[mu],mu)));
            }

            // LePage
            if(lt.c_lp!=0) {
                for (int nu=0;nu<Nd;nu++) {
                    if(mu==nu) continue; 

                    // nu, nu, mu, Back(nu), Back(nu)
                    V[mu] = V[mu] + lt.c_lp*Gimpl::CovShiftForward(U[nu],nu,
                                              Gimpl::CovShiftForward(U[nu],nu,
                                                Gimpl::CovShiftForward(U[mu],mu,
                                                  Gimpl::CovShiftBackward(U[nu],nu,
                                                    Gimpl::CovShiftIdentityBackward(U[nu],nu)))))
                                    // Back(nu), Back(nu), mu, nu, nu
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
    int _HaloDepth=1; 

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
    void ddVprojectU3(GF& u_deriv, GF& u_mu, GF& u_force, RealScalar const delta=5e-5) {

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

                    u2 = u *u; u3 = u2*u; u4 = u3*u; u5 = u4*u; u6 = u5*u; u7 = u6*u; u8 = u7*u;
                    v2 = v *v; v3 = v2*v; v4 = v3*v; v5 = v4*v; v6 = v5*v;
                    w2 = w *w; w3 = w2*w; w4 = w3*w; w5 = w4*w;
        
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
    GF outerProductHISQ(std::vector<FF>& vecx, int l, int sep) {
        
        auto grid   = this->_grid;
        auto gridRB = this->_gridRB;

        GF XY(grid);
        FF X(grid), Y(grid), RB(gridRB);
        LF XYnu(grid), YXnu(grid);
        X = Zero(); Y=Zero(); 
        
        RB=Zero(); 
        pickCheckerboard(Even,RB,vecx[l]);
        setCheckerboard(X,RB);
        RB=Zero();
        pickCheckerboard(Odd ,RB,vecx[l]);
        setCheckerboard(Y,RB);

        XY = Zero(); XYnu = Zero(); YXnu=Zero(); 
        for (int nu = 0; nu < Nd; nu++) {
            YXnu = outerProduct( Cshift(Y,nu,sep) ,X);
            XYnu = outerProduct( Cshift(X,nu,sep) ,Y);
            PokeIndex<LorentzIndex>(XY,(YXnu-XYnu),nu);
        }
        return XY;   
    }


    // We are calculating the force using the rational approximation. The goal is that we can approximate
    // (Mdag M)^(-nf/4) = alpha_0 + sum_l alpha_l/(M^dag M + beta_l). Hence the index l runs over the
    // introduce a different "Naik epsilon" for each M. Hence, we can think of the total application
    // of this operator as having an index inaik, running over the different Naik epsilons; for each inaik
    // there is a possibly different order_inaik, then the operator has an index l running up to order_inaik.
    // All terms with inaik=0 correspond to epsilon_Naik = 0.
    //
    // Intent: OUT--momentum
    //          IN--vecdt: Monte Carlo separation vector times alpha_{inaik,0}. 
    //              vecx: A vector of fermion fields coming from the MILC code. It is organized so that 
    //                    |X_l> = (Mdag M + beta_l)^-1 |Phi> is on even sites, |Y_l>=D|X_l> is on odd sites.
    //                    All the |X_l> for i=0 come first in memory, followed by all the |X_l> with
    //                    i=1 in memory, and so on.
    //              n_orders_naik: Indexed by unique naik epsilon.
    void force(GF& momentum, std::vector<Real> vecdt, std::vector<FF>& vecx, std::vector<int> n_orders_naik) {

        HISQParameters<Real> hp = this->_linkParams;
        auto grid = this->_grid;

        GF XY(grid);       // outer product field
        GF u_force(grid);  // accumulates the force

        momentum = Zero();

        // These four lines control the loop over rational approximation contributions. As explained above,
        // l indexes over both Naik epsilon and rational approximation order.
        int l = 0;
        for (int inaik = 0; inaik < hp.n_naiks; inaik++) {
            int rat_order = n_orders_naik[inaik];
            for (int i=0; i<rat_order; i++) {

                GF temp(grid);

                // ----------------------------------------- NAIK-LINK DERIVATIVE 
                XY = outerProductHISQ(vecx, l, 3);

                std::vector<LF> Wv(Nd, grid);
                std::vector<LF> XYdag(Nd, grid);
                std::vector<LF> ddW(Nd, grid);
                for (int mu = 0; mu < Nd; mu++) {
                    Wv[mu]    = PeekIndex<LorentzIndex>(_Wmu, mu);
                    XYdag[mu] = PeekIndex<LorentzIndex>(XY, mu);
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

                momentum += hp.asqtad_cnaik*vecdt[l]*temp;


                // -------------------------- ONE-LINK DERIVATIVE (OUTER PRODUCT)
                XY = outerProductHISQ(vecx, l, 1); 

                momentum += hp.asqtad_c1*vecdt[l]*XY; // It's not clear to me whether this should be fat7 or asqtad.


                // -------------------------------------------- LEPAGE DERIVATIVE 
                for (int mu = 0; mu < Nd; mu++) {
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

                momentum += hp.asqtad_clp*vecdt[l]*temp;


                // ------------------------------------------- N-LINK DERIVATIVES 
                PaddedCell Ghost(_HaloDepth,grid);
                GF Ughost  = Ghost.Exchange(_Wmu); 
                GF XYghost = Ghost.Exchange(XY);
                GF Fghost  = Ghost.Exchange(u_force);
                std::vector<Coordinate> shifts = createHISQStencil();
                GeneralLocalStencil gStencil(Ughost.Grid(),shifts);

                Fghost = Zero(); 

                for(int mu=0;mu<Nd;mu++) {

                    autoView(U_v , Ughost , AcceleratorRead);
                    autoView(XY_v, XYghost, AcceleratorRead);
                    autoView(F_v , Fghost , AcceleratorWrite);
        
                    int Nsites = U_v.size();
                    auto gStencil_v = gStencil.View(AcceleratorRead);
            
                    typedef decltype(getLink(U_v,gStencil.GetEntry(0,0),0)) U3matrix;
                    typedef decltype(U_v)                                   linkRead;
                    typedef decltype(F_v)                                   linkWrite;
                    typedef decltype(gStencil_v)                            stencilRead;
        
                    accelerator_for(site,Nsites,Simd::Nsimd(),{ // 3-LINK DERIVATIVE 
                        U3matrix U0, U1, U2, U3, U4, U5, XY0, XY1, XY2, XY3, XY4, XY5, res;
                        for(int nu=0;nu<Nd;nu++) {
                            if(nu==mu) continue;
                            int s = HISQStencilIndex(mu,nu);
        
                            auto [x_p_mu, x_p_nu, x, x_p_mu_m_nu, x_m_nu] = get3StaplePoints(gStencil_v,s,site);      
        
                            U0  = getLink(U_v ,x_p_mu     ,nu);
                            U1  = getLink(U_v ,x_p_nu     ,mu);
                            U2  = getLink(U_v ,x          ,nu);
                            U3  = getLink(U_v ,x_p_mu_m_nu,nu);
                            U4  = getLink(U_v ,x_m_nu     ,mu);
                            U5  = getLink(U_v ,x_m_nu     ,nu);
        
                            XY0 = getLink(XY_v,x_p_mu     ,nu);
                            XY1 = getLink(XY_v,x_p_nu     ,mu);
                            XY2 = getLink(XY_v,x          ,nu);
                            XY3 = getLink(XY_v,x_p_mu_m_nu,nu);
                            XY4 = getLink(XY_v,x_m_nu     ,mu);
                            XY5 = getLink(XY_v,x_m_nu     ,nu);

                            res =   adj(XY2)*U1*adj(U0) +     U2* adj(XY1)*adj(U0) +     U2 *U1*    XY0   
                                  +     XY5 *U4*    U3  + adj(U5)*adj(XY4)*    U3  + adj(U5)*U4*adj(XY3);
        
                            setLink(F_v[x->_offset](mu), F_v(x->_offset)(mu) + hp.asqtad_c3*vecdt[l]*adj(res));
                        }              
                    })


//                    accelerator_for(site,Nsites,Simd::Nsimd(),{ // 5-LINK DERIVATIVE
//                        U3matrix U0, U1, U2, U3, U4, U5, XY0, XY1, XY2, XY3, XY4, XY5, V0, V1, V2, V3, V4, V5, S0, S1, S2, S3, S4, S5, res;
//                        int sigmaIndex = 0;
//                        for(int nu=0;nu<Nd;nu++) {
//                            if(nu==mu) continue;
//                            int s = HISQStencilIndex(mu,nu);
//                            for(int rho=0;rho<Nd;rho++) {
//                                if (rho == mu || rho == nu) continue;
//            
//                                auto [x_p_mu, x_p_nu, x, x_p_mu_m_nu, x_m_nu] = get3StaplePoints(gStencil_v,s,site);      
//        
//                                U0  = getLink(U_v       ,x_p_mu     ,nu );
//                                U1  = getLink(U_v       ,x_p_nu     ,mu );
//                                U2  = getLink(U_v       ,x          ,nu );
//                                U3  = getLink(U_v       ,x_p_mu_m_nu,nu );
//                                U4  = getLink(U_v       ,x_m_nu     ,mu );
//                                U5  = getLink(U_v       ,x_m_nu     ,nu );
//            
//                                XY0 = getLink(XY_v      ,x_p_mu     ,nu );
//                                XY2 = getLink(XY_v      ,x          ,nu );
//                                XY3 = getLink(XY_v      ,x_p_mu_m_nu,nu );
//                                XY5 = getLink(XY_v      ,x_m_nu     ,nu );
//
//                                V0  = getLink(dU_3link_v,x_p_mu     ,rho);
//                                V1  = getLink(dU_3link_v,x_p_nu     ,rho);
//                                V2  = getLink(dU_3link_v,x          ,rho);
//                                V3  = getLink(dU_3link_v,x_p_mu_m_nu,rho);
//                                V4  = getLink(dU_3link_v,x_m_nu     ,rho);
//                                V5  = getLink(dU_3link_v,x_m_nu     ,rho);
//
//                                S0  = getLink( U_3link_v,x_p_mu     ,rho);
//                                S2  = getLink( U_3link_v,x          ,rho);
//                                S3  = getLink( U_3link_v,x_p_mu_m_nu,rho);
//                                S5  = getLink( U_3link_v,x_m_nu     ,rho);
//
//                                res =       V2 *U1*adj(U0) +     U2 *V1*adj(U0) +     U2 *U1*adj(V0) +     S2 *U1*    XY0  + adj(XY2)*U1*adj(S0)
//                                      + adj(V5)*U4*    U3  + adj(U5)*V4*    U3  + adj(U5)*U4*    V3  + adj(S5)*U4*adj(XY3) +     XY5 *U4*    S3 ; 
////                                U0  = getLink(U_v       ,x_p_mu     ,nu );
////                                U1  = getLink(U_3link_v ,x_p_nu     ,rho);
////                                U2  = getLink(U_v       ,x          ,nu );
////                                U3  = getLink(U_v       ,x_p_mu_m_nu,nu );
////                                U4  = getLink(U_3link_v ,x_m_nu     ,rho);
////                                U5  = getLink(U_v       ,x_m_nu     ,nu );
////                  
////                                XY0 = getLink(XY_v      ,x_p_mu     ,nu );
////                                V1  = getLink(dU_3link_v,x_p_nu     ,rho);
////                                XY2 = getLink(XY_v      ,x          ,nu );
////                                XY3 = getLink(XY_v      ,x_p_mu_m_nu,nu );
////                                V4  = getLink(dU_3link_v,x_m_nu     ,rho);
////                                XY5 = getLink(XY_v      ,x_m_nu     ,nu );
////
////                                res  =  adj(XY2)*U1*adj(U0) +     U2 *     V1 *adj(U0) +     U2 *U1*    XY0 
////                                     +      XY5 *U4*    U3  + adj(U5)*     V4 *    U3  + adj(U5)*U4*adj(XY3);
//            
//                                if(sigmaIndex<3) {
//                                    setLink(dU_5linkA_v[x->_offset](rho), res);
//                                } else {
//                                    setLink(dU_5linkB_v[x->_offset](rho), res);
//                                }    
//            
//                                setLink(F_v[x->_offset](mu), F_v(x->_offset)(mu) + hp.asqtad_c5*vecdt[l]*adj(res));
//                                sigmaIndex++;
//                            }
//                        }
//                    })
//            
//
//                    accelerator_for(site,Nsites,Simd::Nsimd(),{ // 7-LINK DERIVATIVE 
//                        U3matrix U0, U1, U2, U3, U4, U5, XY0, V1, XY2, XY3, V4, XY5, W;
//                        int sigmaIndex = 0;
//                        for(int nu=0;nu<Nd;nu++) {
//                            if(nu==mu) continue;
//                            int s = HISQStencilIndex(mu,nu);
//                            for(int rho=0;rho<Nd;rho++) {
//                                if (rho == mu || rho == nu) continue;
//                    
//                                auto [x_p_mu, x_p_nu, x, x_p_mu_m_nu, x_m_nu] = get3StaplePoints(gStencil_v,s,site);      
//                    
//                                U0  = getLink(U_v ,x_p_mu,nu);
//                                XY0 = getLink(XY_v,x_p_mu,nu );
//                                if(sigmaIndex<3) {
//                                    U1 = getLink( U_5linkB_v,x_p_nu,rho);
//                                    V1 = getLink(dU_5linkB_v,x_p_nu,rho);
//                                } else {
//                                    U1 = getLink( U_5linkA_v,x_p_nu,rho);
//                                    V1 = getLink(dU_5linkA_v,x_p_nu,rho);
//                                }  
//                                U2  = getLink(U_v ,x          ,nu);
//                                XY2 = getLink(XY_v,x          ,nu);
//                                U3  = getLink(U_v ,x_p_mu_m_nu,nu);
//                                XY3 = getLink(XY_v,x_p_mu_m_nu,nu);
//                                if(sigmaIndex<3) {
//                                    U4 = getLink( U_5linkB_v,x_m_nu,rho);
//                                    V4 = getLink(dU_5linkB_v,x_m_nu,rho);
//                                } else {
//                                    U4 = getLink( U_5linkA_v,x_m_nu,rho);
//                                    V4 = getLink(dU_5linkA_v,x_m_nu,rho);
//                                }  
//                                U5  = getLink(U_v ,x_m_nu,nu);
//                                XY5 = getLink(XY_v,x_m_nu,nu);
//                  
//                                W  =   adj(XY2)*U1*adj(U0) +     U2 *     V1*adj(U0) +     U2 *U1*    XY0 
//                                     +     XY5 *U4*    U3  + adj(U5)*     V4*    U3  + adj(U5)*U4*adj(XY3);                    
//                    
//                                setLink(F_v[x->_offset](mu), F_v(x->_offset)(mu) + hp.fat7_c7*W*vecdt[l]);
//                                sigmaIndex++;
//                            }
//                        }
//                    })                                                                             
                                                                                  
                } // end mu loop
        
                u_force = Ghost.Extract(Fghost);
                momentum += u_force;

                l++;

            }
        }



        // Close the loop: Multiply on the left by U_mu(x)
        LF mom(grid);
        for (int mu = 0; mu < Nd; mu++) {
            mom = PeekIndex<LorentzIndex>(_Umu, mu) * PeekIndex<LorentzIndex>(momentum, mu);
            PokeIndex<LorentzIndex>(momentum, mom, mu);
        }
    }


};


NAMESPACE_END(Grid);
