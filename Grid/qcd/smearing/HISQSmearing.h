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
    @brief Declares classes related to HISQ smearing 
*/


#pragma once
#include <Grid/Grid.h>
#include <Grid/lattice/PaddedCell.h>
#include <Grid/stencil/GeneralLocalStencil.h>


NAMESPACE_BEGIN(Grid);


// TODO: find a way to fold this into the stencil header. need to access grid to get
// Nd, since you don't want to inherit from QCD.h
/*!  @brief append arbitrary shift path to shifts */
template<typename... Args>
void appendShift(std::vector<Coordinate>& shifts, int dir, Args... args) {
    Coordinate shift(Nd,0);
    generalShift(shift, dir, args...); 
    // push_back creates an element at the end of shifts and
    // assigns the data in the argument to it.
    shifts.push_back(shift);
}


/*!  @brief figure out the stencil index from mu and nu */
inline int stencilIndex(int mu, int nu) {
    // Nshifts depends on how you built the stencil
    int Nshifts = 5;
    return Nshifts*nu + Nd*Nshifts*mu;
}


/*!  @brief structure holding the link treatment */
struct SmearingParameters{
    SmearingParameters(){}
    Real c_1;               // 1 link
    Real c_naik;            // Naik term
    Real c_3;               // 3 link
    Real c_5;               // 5 link
    Real c_7;               // 7 link
    Real c_lp;              // 5 link Lepage
    SmearingParameters(Real c1, Real cnaik, Real c3, Real c5, Real c7, Real clp) 
        : c_1(c1),
          c_naik(cnaik),
          c_3(c3),
          c_5(c5),
          c_7(c7),
          c_lp(clp){}
};


/*!  @brief create fat links from link variables */
template<class LGF, class Gimpl> 
class Smear_HISQ_fat : public Gimpl {

private:
    GridCartesian* const _grid;
    SmearingParameters _linkTreatment;

public:

    INHERIT_GIMPL_TYPES(Gimpl);

    // Don't allow default values here.
    Smear_HISQ_fat(GridCartesian* grid, Real c1, Real cnaik, Real c3, Real c5, Real c7, Real clp) 
        : _grid(grid), 
          _linkTreatment(c1,cnaik,c3,c5,c7,clp) {
        assert(Nc == 3 && "HISQ smearing currently implemented only for Nc==3");
        assert(Nd == 4 && "HISQ smearing only defined for Nd==4");
    }

    // Allow to pass a pointer to a C-style, double array for MILC convenience
    Smear_HISQ_fat(GridCartesian* grid, double* coeff) 
        : _grid(grid), 
          _linkTreatment(coeff[0],coeff[1],coeff[2],coeff[3],coeff[4],coeff[5]) {
        assert(Nc == 3 && "HISQ smearing currently implemented only for Nc==3");
        assert(Nd == 4 && "HISQ smearing only defined for Nd==4");
    }

    ~Smear_HISQ_fat() {}

    void smear(LGF& u_smr, LGF& u_thin) const {

        SmearingParameters lt = this->_linkTreatment;

        // Create a padded cell of extra padding depth=1 and fill the padding.
        int depth = 1;
        PaddedCell Ghost(depth,this->_grid);
        LGF Ughost = Ghost.Exchange(u_thin);

        // This is where auxiliary N-link fields and the final smear will be stored. 
        LGF Ughost_fat(Ughost.Grid());
        LGF Ughost_3link(Ughost.Grid());
        LGF Ughost_5linkA(Ughost.Grid());
        LGF Ughost_5linkB(Ughost.Grid());

        // Create 3-link stencil. We allow mu==nu just to make the indexing easier.
        // Shifts with mu==nu will not be used. 
        std::vector<Coordinate> shifts;
        for(int mu=0;mu<Nd;mu++)
        for(int nu=0;nu<Nd;nu++) {
            appendShift(shifts,mu);
            appendShift(shifts,nu);
            appendShift(shifts,NO_SHIFT);
            appendShift(shifts,mu,Back(nu));
            appendShift(shifts,Back(nu));
        }

        // A GeneralLocalStencil has two indices: a site and stencil index 
        GeneralLocalStencil gStencil(Ughost.Grid(),shifts);

        // This is where contributions from the smearing get added together
        Ughost_fat=Zero();

        // This loop handles 3-, 5-, and 7-link constructs, minus Lepage and Naik.
        for(int mu=0;mu<Nd;mu++) {

            // Create the accessors
            autoView(U_v       , Ughost       , CpuRead);
            autoView(U_fat_v   , Ughost_fat   , CpuWrite);
            autoView(U_3link_v , Ughost_3link , CpuWrite);
            autoView(U_5linkA_v, Ughost_5linkA, CpuWrite);
            autoView(U_5linkB_v, Ughost_5linkB, CpuWrite);

            // TODO: This approach is slightly memory inefficient. It uses 25% extra memory 
            Ughost_3link =Zero();
            Ughost_5linkA=Zero();
            Ughost_5linkB=Zero();

            // We infer some types that will be needed in the calculation.
            typedef decltype(gStencil.GetEntry(0,0)) stencilElement;
            typedef decltype(coalescedReadGeneralPermute(U_v[0](0),gStencil.GetEntry(0,0)->_permute,Nd)) U3matrix;
            stencilElement SE0, SE1, SE2, SE3, SE4;
            U3matrix U0, U1, U2, U3, U4, U5, W;

            for(int site=0;site<U_v.size();site++){ // ----------- 3-link
                for(int nu=0;nu<Nd;nu++) {
                    if(nu==mu) continue;
                    int s = stencilIndex(mu,nu);

                    // The stencil gives us support points in the mu-nu plane that we will use to
                    // grab the links we need.
                    SE0 = gStencil.GetEntry(s+0,site); int x_p_mu      = SE0->_offset;
                    SE1 = gStencil.GetEntry(s+1,site); int x_p_nu      = SE1->_offset;
                    SE2 = gStencil.GetEntry(s+2,site); int x           = SE2->_offset;
                    SE3 = gStencil.GetEntry(s+3,site); int x_p_mu_m_nu = SE3->_offset;
                    SE4 = gStencil.GetEntry(s+4,site); int x_m_nu      = SE4->_offset;

                    // When you're deciding whether to take an adjoint, the question is: how is the
                    // stored link oriented compared to the one you want? If I imagine myself travelling
                    // with the to-be-updated link, I have two possible, alternative 3-link paths I can
                    // take, one starting by going to the left, the other starting by going to the right.
                    U0 = coalescedReadGeneralPermute(U_v[x_p_mu     ](nu),SE0->_permute,Nd);
                    U1 = coalescedReadGeneralPermute(U_v[x_p_nu     ](mu),SE1->_permute,Nd);
                    U2 = coalescedReadGeneralPermute(U_v[x          ](nu),SE2->_permute,Nd);
                    U3 = coalescedReadGeneralPermute(U_v[x_p_mu_m_nu](nu),SE3->_permute,Nd);
                    U4 = coalescedReadGeneralPermute(U_v[x_m_nu     ](mu),SE4->_permute,Nd);
                    U5 = coalescedReadGeneralPermute(U_v[x_m_nu     ](nu),SE4->_permute,Nd);

                    //  "left"          "right"
                    W = U2*U1*adj(U0) + adj(U5)*U4*U3;

                    U_3link_v[site](nu) = W;

                    U_fat_v[site](mu) = U_fat_v[site](mu) + lt.c_3*W;
                }
            }

            for(int site=0;site<U_v.size();site++){ // ----------- 5-link
                int sigmaIndex = 0;
                for(int nu=0;nu<Nd;nu++) {
                    if(nu==mu) continue;
                    int s = stencilIndex(mu,nu);
                    for(int rho=0;rho<Nd;rho++) {
                        if (rho == mu || rho == nu) continue;

                        SE0 = gStencil.GetEntry(s+0,site); int x_p_mu      = SE0->_offset;
                        SE1 = gStencil.GetEntry(s+1,site); int x_p_nu      = SE1->_offset;
                        SE2 = gStencil.GetEntry(s+2,site); int x           = SE2->_offset;
                        SE3 = gStencil.GetEntry(s+3,site); int x_p_mu_m_nu = SE3->_offset;
                        SE4 = gStencil.GetEntry(s+4,site); int x_m_nu      = SE4->_offset;

                        U0 = coalescedReadGeneralPermute(      U_v[x_p_mu     ](nu ),SE0->_permute,Nd);
                        U1 = coalescedReadGeneralPermute(U_3link_v[x_p_nu     ](rho),SE1->_permute,Nd);
                        U2 = coalescedReadGeneralPermute(      U_v[x          ](nu ),SE2->_permute,Nd);
                        U3 = coalescedReadGeneralPermute(      U_v[x_p_mu_m_nu](nu ),SE3->_permute,Nd);
                        U4 = coalescedReadGeneralPermute(U_3link_v[x_m_nu     ](rho),SE4->_permute,Nd);
                        U5 = coalescedReadGeneralPermute(      U_v[x_m_nu     ](nu ),SE4->_permute,Nd);

                        W  = U2*U1*adj(U0) + adj(U5)*U4*U3;

                        if(sigmaIndex<3) {
                            U_5linkA_v[site](rho) = W;
                        } else {
                            U_5linkB_v[site](rho) = W;
                        }    

                        U_fat_v[site](mu) = U_fat_v[site](mu) + lt.c_5*W;

                        sigmaIndex++;
                    }
                }
            }

            for(int site=0;site<U_v.size();site++){ // ----------- 7-link
                int sigmaIndex = 0;
                for(int nu=0;nu<Nd;nu++) {
                    if(nu==mu) continue;
                    int s = stencilIndex(mu,nu);
                    for(int rho=0;rho<Nd;rho++) {
                        if (rho == mu || rho == nu) continue;

                        SE0 = gStencil.GetEntry(s+0,site); int x_p_mu      = SE0->_offset;
                        SE1 = gStencil.GetEntry(s+1,site); int x_p_nu      = SE1->_offset;
                        SE2 = gStencil.GetEntry(s+2,site); int x           = SE2->_offset;
                        SE3 = gStencil.GetEntry(s+3,site); int x_p_mu_m_nu = SE3->_offset;
                        SE4 = gStencil.GetEntry(s+4,site); int x_m_nu      = SE4->_offset;

                        U0 = coalescedReadGeneralPermute(U_v[x_p_mu](nu),SE0->_permute,Nd);
                        if(sigmaIndex<3) {
                            U1 = coalescedReadGeneralPermute(U_5linkB_v[x_p_nu](rho),SE1->_permute,Nd);
                        } else {
                            U1 = coalescedReadGeneralPermute(U_5linkA_v[x_p_nu](rho),SE1->_permute,Nd);
                        }  
                        U2 = coalescedReadGeneralPermute(U_v[x](nu),SE2->_permute,Nd);
                        U3 = coalescedReadGeneralPermute(U_v[x_p_mu_m_nu](nu),SE3->_permute,Nd);
                        if(sigmaIndex<3) {
                            U4 = coalescedReadGeneralPermute(U_5linkB_v[x_m_nu](rho),SE4->_permute,Nd);
                        } else {
                            U4 = coalescedReadGeneralPermute(U_5linkA_v[x_m_nu](rho),SE4->_permute,Nd);
                        }  
                        U5 = coalescedReadGeneralPermute(U_v[x_m_nu](nu),SE4->_permute,Nd);

                        W  = U2*U1*adj(U0) + adj(U5)*U4*U3;

                        U_fat_v[site](mu) = U_fat_v[site](mu) + lt.c_7*W;

                        sigmaIndex++;
                    }
                }
            }

        } // end mu loop

        // c1, c3, c5, c7 construct contributions
        u_smr = Ghost.Extract(Ughost_fat) + lt.c_1*u_thin;

        // Load up U and V std::vectors to access thin and smeared links.
        std::vector<LatticeColourMatrix> U(Nd, u_thin.Grid());
        std::vector<LatticeColourMatrix> V(Nd, u_smr.Grid());
        for (int mu = 0; mu < Nd; mu++) {
            U[mu] = PeekIndex<LorentzIndex>(u_thin, mu);
            V[mu] = PeekIndex<LorentzIndex>(u_smr, mu);
        }

        // Compute LePage term from U_thin:
        for(int mu=0;mu<Nd;mu++) {
            for (int nu_h=1;nu_h<Nd;nu_h++) {
                int nu=(mu+nu_h)%Nd;
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

        // Put V back into u_smr.
        for (int mu = 0; mu < Nd; mu++) {
            PokeIndex<LorentzIndex>(u_smr, V[mu], mu);
        }
    };

//    void derivative(const GaugeField& Gauge) const {
//    };
};


/*!  @brief create long links from link variables. */
template<class LGF>
class Smear_HISQ_Naik {

private:
    GridCartesian* const _grid;

public:

    // Eventually this will take, e.g., coefficients as argument 
    Smear_HISQ_Naik(GridCartesian* grid) : _grid(grid) {
        assert(Nc == 3 && "HISQ smearing currently implemented only for Nc==3");
        assert(Nd == 4 && "HISQ smearing only defined for Nd==4");
    }

    ~Smear_HISQ_Naik() {}

//    void smear(LGF& u_smr, const LGF& U) const {
//    };

//    void derivative(const GaugeField& Gauge) const {
//    };
};


NAMESPACE_END(Grid);