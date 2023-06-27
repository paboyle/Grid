/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/smearing/StoutSmearing.h

Copyright (C) 2019

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

#define BACKWARD_CONST 16
#define NO_SHIFT -1

NAMESPACE_BEGIN(Grid);


// This is to optimize the SIMD (will also need to be in the class, at least for now)
template<class vobj> void gpermute(vobj & inout,int perm) {
    vobj tmp=inout;
    if (perm & 0x1) {permute(inout,tmp,0); tmp=inout;}
    if (perm & 0x2) {permute(inout,tmp,1); tmp=inout;}
    if (perm & 0x4) {permute(inout,tmp,2); tmp=inout;}
    if (perm & 0x8) {permute(inout,tmp,3); tmp=inout;}
}


/*!  @brief signals that you want to go backwards in direction dir */
inline int Back(const int dir) {
    // generalShift will use BACKWARD_CONST to determine whether we step forward or 
    // backward. Should work as long as BACKWARD_CONST > Nd. 
    return dir + BACKWARD_CONST;
}


/*!  @brief shift one unit in direction dir */
void generalShift(Coordinate& shift, int dir) {
    if (dir >= BACKWARD_CONST) {
        dir -= BACKWARD_CONST;
        shift[dir]+=-1;
    } else if (dir == NO_SHIFT) {
        ; // do nothing
    } else {
        shift[dir]+=1;
    }
}


/*!  @brief follow a path of directions, shifting one unit in each direction */
template<typename... Args>
void generalShift(Coordinate& shift, int dir, Args... args) {
    if (dir >= BACKWARD_CONST) {
        dir -= BACKWARD_CONST;
        shift[dir]+=-1;
    } else if (dir == NO_SHIFT) {
        ; // do nothing
    } else {
        shift[dir]+=1;
    }
    generalShift(shift, args...);
}


/*!  @brief append arbitrary shift path to shifts */
template<typename... Args>
void appendShift(std::vector<Coordinate>& shifts, int dir, Args... args) {
    Coordinate shift(Nd,0);
    generalShift(shift, dir, args...); 
    // push_back creates an element at the end of shifts and
    // assigns the data in the argument to it.
    shifts.push_back(shift);
}


/*!  @brief structure holding the link treatment */
struct SmearingParameters {
    SmearingParameters(){}
    Real c_1;               // 1 link
    Real c_naik;            // Naik term
    Real c_3;               // 3 link
    Real c_5;               // 5 link
    Real c_7;               // 7 link
    Real c_lp;              // 5 link Lepage
    SmearingParameters(Real c1, Real cnaik, Real c3, Real c5, Real c7, Real clp) :
        c_1(c1),
        c_naik(cnaik),
        c_3(c3),
        c_5(c5),
        c_7(c7),
        c_lp(clp){}
};


/*!  @brief create fat links from link variables */
template<class LGF> // TODO: change to Gimpl?
class Smear_HISQ_fat {

private:
    GridCartesian* const _grid;
    SmearingParameters _linkTreatment;

public:

    // Don't allow default values here.
    Smear_HISQ_fat(GridCartesian* grid, Real c1, Real cnaik, Real c3, Real c5, Real c7, Real clp) 
        : _grid(grid), 
          _linkTreatment(c1,cnaik,c3,c5,c7,clp) {
        assert(Nc == 3 && "HISQ smearing currently implemented only for Nc==3");
    }

    // Allow to pass a pointer to a C-style, double array for MILC convenience
    Smear_HISQ_fat(GridCartesian* grid, double* coeff) 
        : _grid(grid), 
          _linkTreatment(coeff[0],coeff[1],coeff[2],coeff[3],coeff[4],coeff[5]) {
        assert(Nc == 3 && "HISQ smearing currently implemented only for Nc==3");
    }

    ~Smear_HISQ_fat() {}

    void smear(LGF& u_smr, LGF& u_thin) const {

        SmearingParameters lt = this->_linkTreatment;

        // We create a cell with extra padding 2. This allows us to capture the LePage
        // term without needing to save intermediate gauge fields or extra halo exchanges.
        // The tradeoff is that we compute extra constructs in the padding. 
        int depth = 2;
        PaddedCell Ghost(depth,this->_grid);
        LGF Ughost = Ghost.Exchange(u_thin);

        // Array for <tr U_mu_nu>(x)
        GridBase *GhostGrid = Ughost.Grid();
        LatticeComplex gplaq(GhostGrid); 

        // This is where the 3-link constructs will be stored
        LGF Ughost_fat(Ughost.Grid());

        // Create a stencil, which is a collection of sites neighboring some initial site.
        std::vector<Coordinate> shifts;
        for(int mu=0;mu<Nd;mu++)
        for(int nu=0;nu<Nd;nu++) {
            if(mu==nu) continue;
            appendShift(shifts,mu);
            appendShift(shifts,nu);
            appendShift(shifts,NO_SHIFT);
            appendShift(shifts,mu,Back(nu));
            appendShift(shifts,Back(nu));
        }
        GeneralLocalStencil gStencil(GhostGrid,shifts);

        Ughost_fat=Zero();

        // Create the accessors, here U_v and U_fat_v 
        autoView(U_v    , Ughost    , CpuRead);
        autoView(U_fat_v, Ughost_fat, CpuWrite);

        // This is a loop over local sites.
        for(int ss=0;ss<U_v.size();ss++){

            // This is the stencil index. It increases as we make our way through the spacetime sites,
            // plaquette orientations, and as we travel around a plaquette.
            int s=0;

            for(int mu=0;mu<Nd;mu++)
            for(int nu=0;nu<Nd;nu++) {
                if(mu==nu) continue;

                // x+mu 
                // m+nu 
                // x
                // shift_munu; shift_munu[mu]= 1; shift_munu[nu]=-1;
                // x-nu 
                auto SE0 = gStencil.GetEntry(s+0,ss);
                auto SE1 = gStencil.GetEntry(s+1,ss);
                auto SE2 = gStencil.GetEntry(s+2,ss);
                auto SE3 = gStencil.GetEntry(s+3,ss);
                auto SE4 = gStencil.GetEntry(s+4,ss);

                // Each offset corresponds to a site around the plaquette.
                int o0 = SE0->_offset;
                int o1 = SE1->_offset;
                int o2 = SE2->_offset;
                int o3 = SE3->_offset;
                int o4 = SE4->_offset;

                // When you're deciding whether to take an adjoint, the question is: how is the
                // stored link oriented compared to the one you want? If I imagine myself travelling
                // with the to-be-updated link, I have two possible, alternative 3-link paths I can
                // take, one starting by going to the left, the other starting by going to the right.
                auto U0 = adj(U_v[o0](nu));
                auto U1 = U_v[o1](mu);
                auto U2 = U_v[o2](nu);

                gpermute(U0,SE0->_permute);
                gpermute(U1,SE1->_permute);
                gpermute(U2,SE2->_permute);

                auto U3 = U_v[o3](nu);
                auto U4 = U_v[o4](mu);
                auto U5 = adj(U_v[o4](nu)); 

                gpermute(U3,SE3->_permute);
                gpermute(U4,SE4->_permute);
                gpermute(U4,SE4->_permute);

                //       "left"     "right"
                auto W = U2*U1*U0 + U5*U4*U3;
                U_fat_v[ss](mu) = U_fat_v[ss](mu) + W;

                s=s+5;
            }
        }

        u_smr = lt.c_3*Ghost.Extract(Ughost_fat) + lt.c_1*u_thin;
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
    }

    ~Smear_HISQ_Naik() {}

    void smear(LGF& u_smr, const LGF& U) const {
        
        int depth = 1;
        PaddedCell Ghost(depth,this->_grid);
        LGF Ughost = Ghost.Exchange(u_smr);
    
        GridBase *GhostGrid = Ughost.Grid();
        LatticeComplex gplaq(GhostGrid); 
    
        LGF Ughost_naik(Ughost.Grid());

        std::vector<Coordinate> shifts;
        for(int mu=0;mu<Nd;mu++){
            for(int nu=mu+1;nu<Nd;nu++){
                // forward shifts
                Coordinate x(Nd,0);
                Coordinate shift_mu(Nd,0); shift_mu[mu]=1;
                Coordinate shift_nu(Nd,0); shift_nu[nu]=1;
                // push_back creates an element at the end of shifts and
                // assigns the data in the argument to it.
                shifts.push_back(shift_mu);
                shifts.push_back(shift_nu);
                shifts.push_back(x);
                // reverse shifts
                shift_nu[nu]=-1;
                Coordinate shift_munu(Nd,0); shift_munu[mu]=1; shift_munu[nu]=-1;
                shifts.push_back(shift_munu);
                shifts.push_back(shift_nu); // in principle you don't need both of these grid points,
                shifts.push_back(shift_nu); // but it helps the reader keep track of offsets
            }
        }
        GeneralLocalStencil gStencil(GhostGrid,shifts);

        Ughost_naik=Zero();

        // Create the accessors, here U_v and U_fat_v 
        autoView(U_v     , Ughost     , CpuRead);
        autoView(U_naik_v, Ughost_naik, CpuWrite);

        // This is a loop over local sites.
        for(int ss=0;ss<U_v.size();ss++){

            // This is the stencil index. It increases as we make our way through the spacetime sites,
            // plaquette orientations, and as we travel around a plaquette.
            int s=0;

            for(int mu=0;mu<Nd;mu++){
                for(int nu=mu+1;nu<Nd;nu++){

                    // shift_mu; shift_mu[mu]=1
                    // shift_nu; shift_nu[nu]=1
                    // x
                    // shift_munu; shift_munu[mu]= 1; shift_munu[nu]=-1;
                    // shift_nu  ;   shift_nu[nu]=-1;
                    // shift_nu  ;   shift_nu[nu]=-1;
                    auto SE0 = gStencil.GetEntry(s+0,ss);
                    auto SE1 = gStencil.GetEntry(s+1,ss);
                    auto SE2 = gStencil.GetEntry(s+2,ss);
                    auto SE3 = gStencil.GetEntry(s+3,ss);
                    auto SE4 = gStencil.GetEntry(s+4,ss);
                    auto SE5 = gStencil.GetEntry(s+5,ss);

                    // Each offset corresponds to a site around the plaquette.
                    int o0 = SE0->_offset;
                    int o1 = SE1->_offset;
                    int o2 = SE2->_offset;
                    int o3 = SE3->_offset;
                    int o4 = SE4->_offset;
                    int o5 = SE5->_offset;

                    auto U0 = U_v[o0](nu);
                    auto U1 = adj(U_v[o1](mu));
                    auto U2 = adj(U_v[o2](nu));

                    gpermute(U0,SE0->_permute);
                    gpermute(U1,SE1->_permute);
                    gpermute(U2,SE2->_permute);

                    auto U3 = adj(U_v[o3](nu));
                    auto U4 = adj(U_v[o4](mu));
                    auto U5 = U_v[o5](nu); 

                    gpermute(U3,SE3->_permute);
                    gpermute(U4,SE4->_permute);
                    gpermute(U5,SE5->_permute);

                    // Forward contribution from this orientation
                    auto W = U0*U1*U2;
                    U_naik_v[ss](mu) = U_naik_v[ss](mu) + W;

                    // Backward contribution from this orientation
                    W = U3*U4*U5;
                    U_naik_v[ss](mu) = U_naik_v[ss](mu) + W;

                    s=s+6;
                }
            }
        }

        // Here is my understanding of this part: The padded cell has its own periodic BCs, so
        // if I take a step to the right at the right-most side of the cell, I end up on the
        // left-most side. This means that the plaquettes in the padding are wrong. Luckily
        // all we care about are the plaquettes in the cell, which we obtain from Extract.
        u_smr = Ghost.Extract(Ughost_naik);
    };

//    void derivative(const GaugeField& Gauge) const {
//    };
};


NAMESPACE_END(Grid);
