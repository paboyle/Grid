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

// things like @brief are seen by things like doxygen and javadocs

#pragma once

#include <Grid/Grid.h>
#include <Grid/lattice/PaddedCell.h>
#include <Grid/stencil/GeneralLocalStencil.h>

NAMESPACE_BEGIN(Grid);


// This is to optimize the SIMD (will also need to be in the class, at least for now)
template<class vobj> void gpermute(vobj & inout,int perm) {
    vobj tmp=inout;
    if (perm & 0x1) {permute(inout,tmp,0); tmp=inout;}
    if (perm & 0x2) {permute(inout,tmp,1); tmp=inout;}
    if (perm & 0x4) {permute(inout,tmp,2); tmp=inout;}
    if (perm & 0x8) {permute(inout,tmp,3); tmp=inout;}
}



/*!  @brief create fat links from link variables. */
class Smear_HISQ_fat {

private:
    GridCartesian* const _grid;

public:

    // Eventually this will take, e.g., coefficients as argument 
    Smear_HISQ_fat(GridCartesian* grid) : _grid(grid) {
        assert(Nc == 3 && "HISQ smearing currently implemented only for Nc==3");
    }

    ~Smear_HISQ_fat() {}

    void smear(LatticeGaugeField& u_smr, const LatticeGaugeField& U) const {
        
        // Create a padded cell of extra padding depth=1
        int depth = 1;
        PaddedCell Ghost(depth,this->_grid);
        LatticeGaugeField Ughost = Ghost.Exchange(u_smr);
    
        // Array for <tr U_mu_nu>(x)
        GridBase *GhostGrid = Ughost.Grid();
        LatticeComplex gplaq(GhostGrid); 
    
        // This is where the 3-link constructs will be stored
        LatticeGaugeField Ughost_fat(Ughost.Grid());

        // Create 3-link stencil (class will build its own stencils)
        // writing your own stencil, you're hard-coding the periodic BCs, so you don't need
        // the policy-based stuff, at least for now
        std::vector<Coordinate> shifts;
        for(int mu=0;mu<Nd;mu++){
            for(int nu=mu+1;nu<Nd;nu++){
                // forward shifts
                Coordinate shift_0(Nd,0);
                Coordinate shift_mu(Nd,0); shift_mu[mu]=1;
                Coordinate shift_nu(Nd,0); shift_nu[nu]=1;
                // push_back creates an element at the end of shifts and
                // assigns the data in the argument to it.
                shifts.push_back(shift_mu);
                shifts.push_back(shift_nu);
                shifts.push_back(shift_0);
                // reverse shifts
                shift_nu[nu]=-1;
                Coordinate shift_munu(Nd,0); shift_munu[mu]=1; shift_munu[nu]=-1;
                shifts.push_back(shift_munu);
                shifts.push_back(shift_nu); // in principle you don't need both of these grid points,
                shifts.push_back(shift_nu); // but it helps the reader keep track of offsets
            }
        }
        GeneralLocalStencil gStencil(GhostGrid,shifts);
    
        Ughost_fat=Zero();
    
        // Create the accessors, here U_v and U_fat_v 
        autoView(U_v      , Ughost      , CpuRead);
        autoView(U_fat_v, Ughost_fat, CpuWrite);
    
        // This is a loop over local sites.
        for(int ss=0;ss<U_v.size();ss++){
    
            // This is the stencil index. It increases as we make our way through the spacetime sites,
            // plaquette orientations, and as we travel around a plaquette.
            int s=0;
        
            for(int mu=0;mu<Nd;mu++){
                for(int nu=mu+1;nu<Nd;nu++){
    
                    // shift_mu; shift_mu[mu]=1
                    // shift_nu; shift_nu[nu]=1
                    // shift_0
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
                    U_fat_v[ss](mu) = U_fat_v[ss](mu) + W;
    
                    // Backward contribution from this orientation
                    W = U3*U4*U5;
                    U_fat_v[ss](mu) = U_fat_v[ss](mu) + W;
    
                    s=s+6;
                }
            }
        }
    
        // Here is my understanding of this part: The padded cell has its own periodic BCs, so
        // if I take a step to the right at the right-most side of the cell, I end up on the
        // left-most side. This means that the plaquettes in the padding are wrong. Luckily
        // all we care about are the plaquettes in the cell, which we obtain from Extract.
        u_smr = Ghost.Extract(Ughost_fat);
    };


    // I guess the way this will go is:
    // 1. 3-link smear
    // 2. exchange
    // 3. 5-link calculated from 3-link

//    void derivative(const GaugeField& Gauge) const {
//    };
};



/*!  @brief create long links from link variables. */
class Smear_HISQ_Naik {

private:
    GridCartesian* const _grid;

public:

    // Eventually this will take, e.g., coefficients as argument 
    Smear_HISQ_Naik(GridCartesian* grid) : _grid(grid) {
        assert(Nc == 3 && "HISQ smearing currently implemented only for Nc==3");
    }

    ~Smear_HISQ_Naik() {}

    void smear(LatticeGaugeField& u_smr, const LatticeGaugeField& U) const {
        
        int depth = 1;
        PaddedCell Ghost(depth,this->_grid);
        LatticeGaugeField Ughost = Ghost.Exchange(u_smr);
    
        GridBase *GhostGrid = Ughost.Grid();
        LatticeComplex gplaq(GhostGrid); 
    
        LatticeGaugeField Ughost_naik(Ughost.Grid());

        std::vector<Coordinate> shifts;
        for(int mu=0;mu<Nd;mu++){
            for(int nu=mu+1;nu<Nd;nu++){
                // forward shifts
                Coordinate shift_0(Nd,0);
                Coordinate shift_mu(Nd,0); shift_mu[mu]=1;
                Coordinate shift_nu(Nd,0); shift_nu[nu]=1;
                // push_back creates an element at the end of shifts and
                // assigns the data in the argument to it.
                shifts.push_back(shift_mu);
                shifts.push_back(shift_nu);
                shifts.push_back(shift_0);
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
                    // shift_0
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
