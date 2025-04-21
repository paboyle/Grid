/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./tests/smearing/Test_HISQ_force.cc

Copyright (C) 2024

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
    @file Test_HISQ_force.cc
    @brief test of the HISQ fermion force 
*/


#include <Grid/Grid.h>
#include <Grid/lattice/PaddedCell.h>
#include <Grid/stencil/GeneralLocalStencil.h>
#include <Grid/qcd/smearing/HISQSmearing.h>
using namespace Grid;


//#define USE_DOUBLE true 
#define USE_DOUBLE true 

#if USE_DOUBLE 
    #define PREC double
    typedef LatticeGaugeFieldD LGF;
    typedef StaggeredImplD GIMPL;
    typedef vComplexD COMP;
#else
    #define PREC float
    typedef LatticeGaugeFieldF LGF;
    typedef StaggeredImplF GIMPL;
    typedef vComplexF COMP;
#endif 

typedef typename GIMPL::FermionField FF;




// This is a sort of contrived test situation. The goal is to make sure the fermion force
// code is stable against future changes and get an idea how the HISQ force interface works.
bool testForce(GridCartesian& GRID, LGF Umu, LGF Ucontrol,
               Real fat7_c1  , Real fat7_c3  , Real fat7_c5  , Real fat7_c7  , Real cnaik,
               Real asqtad_c1, Real asqtad_c3, Real asqtad_c5, Real asqtad_c7, Real asqtad_clp) {

    LGF Vmu(&GRID), Wmu(&GRID), Nmu(&GRID), Umom(&GRID);

    // The n_orders_naik array is indexed according to the unique Naik epsilon values. By convention, 
    // index 0 always corresponds to zero Naik epsilon. n_orders_naik[0] is the sum of the RHMC 
    // molecular dynamics orders of all the pseudofermion fermions with zero Naik epsilon at the 
    // head of the rational function parameter file. n_orders_naik[1] is the sum of orders for the 
    // first nonzero Naik epsilon. The sum of all the n_orders_naik elements equals the size of the 
    // vecx and vecdt arrays. Elements of those arrays are grouped according to n_orders_naik, so 
    // the first n_orders_naik[0] elements of vecx and epsv correspond to the pseudofermions with 
    // nonzero Naik epsilon. That group lumps together terms for each of the zero-epsilon 
    // pseudofermions.
    int n_naiks = 1; // Just a charm 
    std::array<Real,GRID_MAX_NAIK> eps_naik = {0,0,0}; 
    std::vector<int> n_orders_naik = {1,1};  
    std::vector<Real> vecdt = {0.1,0.1};  

    HISQParameters<Real> hisq_param(n_naiks  , eps_naik ,
                                    fat7_c1  , fat7_c3  , fat7_c5  , fat7_c7  , 0.,
                                    asqtad_c1, asqtad_c3, asqtad_c5, asqtad_c7, asqtad_clp,
                                    cnaik    , 0.       , 0.);
    HISQReunitSVDParameters<Real> hisq_reunit_svd(false, false, 1, 1, 1);

    Smear_HISQ<GIMPL> fat7(&GRID,fat7_c1,0.,fat7_c3,fat7_c5,fat7_c7,0.);
    fat7.smear(Vmu,Nmu,Umu); // Populate fat7 and Naik links Vmu and Nmu
    fat7.projectU3(Wmu,Vmu); // Populate U(3) projection Wmu
    Force_HISQ<GIMPL> hisq_force(&GRID, hisq_param, Wmu, Vmu, Umu, hisq_reunit_svd);

    std::vector<int> seeds({1,2,3,4});
    GridParallelRNG pRNG(&GRID);  pRNG.SeedFixedIntegers(seeds);

    // Construct a vecx. Eventually it would be nice if this generated the correct distribution,
    // but for now use Gaussian random variables as a placeholder.
    std::vector<FF> vecx;
    int l = 0;
    for (int inaik = 0; inaik < hisq_param.n_naiks; inaik++) {
        int rat_order = n_orders_naik[inaik];
        FF PHI(&GRID);
        PHI = Zero();
        for (int i=0; i<rat_order; i++) {
            gaussian(pRNG,PHI);
            vecx.push_back(PHI); 
        }
    }

    hisq_force.force(Umom,vecdt,vecx,n_orders_naik);
    LGF diff(&GRID);
    diff = Ucontrol-Umom;
    auto absDiff = norm2(diff)/norm2(Ucontrol);
    if (absDiff < 1e-30) {
        Grid_pass(" |Umu-Umom|/|Umu| = ",absDiff);
        return true;
    } else {
        Grid_error(" |Umu-Umom|/|Umu| = ",absDiff);
        return false;
    }
//    NerscIO::writeConfiguration(Umom,"nersc.l8t4b3360.Umom.XY.control");
//    NerscIO::writeConfiguration(Umom,"nersc.l8t4b3360.Umom.3link.control");
//    NerscIO::writeConfiguration(Umom,"nersc.l8t4b3360.Umom.lp.control");
//    NerscIO::writeConfiguration(Umom,"nersc.l8t4b3360.Umom.naik.control");
//    NerscIO::writeConfiguration(Umom,"nersc.l8t4b3360.Umom.5link.control");
    return true;
}


// Intent: IN--Umu: thin links
bool testddUProj(GridCartesian& GRID, LGF Umu, LGF Ucontrol) {

    LGF Uforce(&GRID);

    int n_naiks = 1;
    std::array<Real,GRID_MAX_NAIK> eps_naik = {0,0,0}; 

    Real fat7_c1  = 1/8.;                  Real asqtad_c1  = 1.;
    Real fat7_c3  = 1/16.;                 Real asqtad_c3  = 1/16.;
    Real fat7_c5  = 1/64.;                 Real asqtad_c5  = 1/64.;
    Real fat7_c7  = 1/384.;                Real asqtad_c7  = 1/384.;
    Real cnaik    = -1/24.+eps_naik[0]/8;  Real asqtad_clp = -1/8.;


    LGF Vmu(&GRID), Wmu(&GRID), Nmu(&GRID);
    Smear_HISQ<GIMPL> fat7(&GRID,fat7_c1,0.,fat7_c3,fat7_c5,fat7_c7,0.);

    fat7.smear(Vmu,Nmu,Umu); // Populate fat7 and Naik links Vmu and Nmu
    fat7.projectU3(Wmu,Vmu); // Populate U(3) projection Wmu

    HISQParameters<Real> hisq_param(n_naiks  , eps_naik ,
                                    fat7_c1  , fat7_c3  , fat7_c5  , fat7_c7  , 0.,
                                    asqtad_c1, asqtad_c3, asqtad_c5, asqtad_c7, asqtad_clp,
                                    cnaik    , 0.       , 0.);

    HISQReunitSVDParameters<Real> hisq_reunit_svd(false, false, 1, 1, 1);

    Force_HISQ<GIMPL> hisq_force(&GRID, hisq_param, Wmu, Vmu, Umu, hisq_reunit_svd);

    LGF diff(&GRID);
    hisq_force.ddVprojectU3(Uforce, Umu, Umu, 5e-5);
    diff = Ucontrol-Uforce;
    auto absDiff = norm2(diff)/norm2(Ucontrol);
    if (absDiff < 1e-30) {
        Grid_pass(" |Umu-Usmr|/|Umu| = ",absDiff);
        return true;
    } else {
        Grid_error(" |Umu-Usmr|/|Umu| = ",absDiff);
        return false;
    }
//    NerscIO::writeConfiguration(Uforce,"nersc.l8t4b3360.ddVU3");
}


int main (int argc, char** argv) {

    // Params for the test.
    int Ns = 8;
    int Nt = 4;
    Coordinate latt_size(Nd,0); latt_size[0]=Ns; latt_size[1]=Ns; latt_size[2]=Ns; latt_size[3]=Nt;
    std::string conf_in  = "nersc.l8t4b3360";
    int threads          = GridThread::GetThreads();

    // Initialize the Grid
    Grid_init(&argc,&argv);
    Coordinate simd_layout = GridDefaultSimd(Nd,COMP::Nsimd());
    Coordinate mpi_layout  = GridDefaultMpi();
    Grid_log("mpi     = ",mpi_layout);
    Grid_log("simd    = ",simd_layout);
    Grid_log("latt    = ",latt_size);
    Grid_log("threads = ",threads);
    GridCartesian GRID(latt_size,simd_layout,mpi_layout);

    LGF Umu(&GRID), Ucontrol(&GRID);

#if USE_DOUBLE // NerscIO is hard-coded to double.

    // Read the configuration into Umu
    FieldMetaData header;
    NerscIO::readConfiguration(Umu, header, conf_in);

    bool pass=true;

    // Check derivative of projection
    NerscIO::readConfiguration(Ucontrol, header, "nersc.l8t4b3360.ddVU3.control");
    pass *= testddUProj(GRID,Umu,Ucontrol);

    // Check the 1-link (outer product) 
    NerscIO::readConfiguration(Ucontrol, header, "nersc.l8t4b3360.Umom.XY.control");
    pass *= testForce(GRID, Umu, Ucontrol,
                      1, 0, 0, 0, 0,
                      1, 0, 0, 0, 0 );

    // Check the 3-link
    NerscIO::readConfiguration(Ucontrol, header, "nersc.l8t4b3360.Umom.3link.control");
    pass *= testForce(GRID, Umu, Ucontrol,
                      1, 0, 0, 0, 0,
                      0, 1, 0, 0, 0 );

    // Check the LePage-link 
    NerscIO::readConfiguration(Ucontrol, header, "nersc.l8t4b3360.Umom.lp.control");
    pass *= testForce(GRID, Umu, Ucontrol,
                      1, 0, 0, 0, 0,
                      0, 0, 0, 0, 1 );

    // Check the Naik-link 
    NerscIO::readConfiguration(Ucontrol, header, "nersc.l8t4b3360.Umom.naik.control");
    pass *= testForce(GRID, Umu, Ucontrol,
                      1, 0, 0, 0, 1,
                      0, 0, 0, 0, 0 );

//    // Check the 5-link
//    NerscIO::readConfiguration(Ucontrol, header, "nersc.l8t4b3360.Umom.5link.control");
//    pass *= testForce(GRID, Umu, Ucontrol,
//                      0, 0, 1, 0, 0,
//                      0, 0, 0, 0, 0 );

    if(pass){
        Grid_pass("All tests passed.");
    } else {
        Grid_error("At least one test failed.");
    }

#endif

    Grid_finalize();
}