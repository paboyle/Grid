/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./tests/smearing/Test_fatLinks.cc

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
    @file Test_fatLinks.cc
    @brief test of the HISQ smearing 
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
    typedef PeriodicGimplD GIMPL;
    typedef vComplexD COMP;
#else
    #define PREC float
    typedef LatticeGaugeFieldF LGF;
    typedef PeriodicGimplF GIMPL;
    typedef vComplexF COMP;
#endif 


bool testSmear(GridCartesian& GRID, LGF Umu, LGF Usmr, LGF Unaik, 
               LGF Ucontrol, PREC c1, PREC cnaik, PREC c3, PREC c5, PREC c7, PREC clp) {
    Smear_HISQ<GIMPL> hisq_fat(&GRID,c1,cnaik,c3,c5,c7,clp);
    LGF diff(&GRID), Uproj(&GRID);
    hisq_fat.smear(Usmr, Unaik, Umu);
    bool result;
    if (cnaik < 1e-30) { // Testing anything but Naik term
        diff = Ucontrol-Usmr;
        auto absDiff = norm2(diff)/norm2(Ucontrol);
        if (absDiff < 1e-30) {
            Grid_pass(" |Umu-Usmr|/|Umu| = ",absDiff);
            result = true;
        } else {
            Grid_error(" |Umu-Usmr|/|Umu| = ",absDiff);
            result = false;
        }
    } else { // Testing Naik specifically
        diff = Ucontrol-Unaik;
        auto absDiff = norm2(diff)/norm2(Ucontrol);
        if (absDiff < 1e-30) {
            Grid_pass(" |Umu-Unaik|/|Umu| = ",absDiff);
            result = true;
        } else {
            Grid_error(" |Umu-Unaik|/|Umu| = ",absDiff);
            result = false;
        }
        hisq_fat.projectU3(Uproj,Ucontrol);
//        NerscIO::writeConfiguration(Unaik,"nersc.l8t4b3360.naik");
    }
    return result;
}


void hotStartSmear(GridCartesian& GRID) {
    LGF Uproj(&GRID), Uhot(&GRID);
    GridParallelRNG pRNG(&GRID); pRNG.SeedFixedIntegers(std::vector<int>({111,222,333,444}));
    SU<Nc>::HotConfiguration(pRNG,Uhot);
    Smear_HISQ<GIMPL> hisq_fat(&GRID,1/8.,0.,1/16.,1/64.,1/384.,-1/8.);
    hisq_fat.projectU3(Uproj,Uhot);
    Grid_log("norm2(Uproj) = ",norm2(Uproj)/(Nc*Nd*GRID.gSites()));
}


int main (int argc, char** argv) {

    // Params for the test.
    int Ns = 8;
    int Nt = 4;
    Coordinate latt_size(Nd,0); latt_size[0]=Ns; latt_size[1]=Ns; latt_size[2]=Ns; latt_size[3]=Nt;
    std::string conf_in  = "nersc.l8t4b3360";
    int threads          = GridThread::GetThreads();

    if (sizeof(PREC)==4) {
        Grid_log("Run in single precision.");
    } else { 
        Grid_log("Run in double precision.");
    }

    // Initialize the Grid
    Grid_init(&argc,&argv);
    Coordinate simd_layout = GridDefaultSimd(Nd,COMP::Nsimd());
    Coordinate mpi_layout  = GridDefaultMpi();
    Grid_log("mpi     = ",mpi_layout);
    Grid_log("simd    = ",simd_layout);
    Grid_log("latt    = ",latt_size);
    Grid_log("threads = ",threads);
    GridCartesian GRID(latt_size,simd_layout,mpi_layout);

    LGF Umu(&GRID), Usmr(&GRID), Unaik(&GRID), Ucontrol(&GRID);

#if USE_DOUBLE // NerscIO is hard-coded to double.

    // Read the configuration into Umu
    FieldMetaData header;
    NerscIO::readConfiguration(Umu, header, conf_in);

    bool pass=true;

    // Carry out various tests    
    NerscIO::readConfiguration(Ucontrol, header, "nersc.l8t4b3360.357lplink.control");
    pass *= testSmear(GRID,Umu,Usmr,Unaik,Ucontrol,1/8.,0.,1/16.,1/64.,1/384.,-1/8.);
    NerscIO::readConfiguration(Ucontrol, header, "nersc.l8t4b3360.357link.control");
    pass *= testSmear(GRID,Umu,Usmr,Unaik,Ucontrol,1/8.,0.,1/16.,1/64.,1/384.,0.);
    NerscIO::readConfiguration(Ucontrol, header, "nersc.l8t4b3360.35link.control");
    pass *= testSmear(GRID,Umu,Usmr,Unaik,Ucontrol,1/8.,0.,1/16.,1/64.,0.,0.);
    NerscIO::readConfiguration(Ucontrol, header, "nersc.l8t4b3360.3link.control");
    pass *= testSmear(GRID,Umu,Usmr,Unaik,Ucontrol,1/8.,0.,1/16.,0.,0.,0.);
    NerscIO::readConfiguration(Ucontrol, header, "nersc.l8t4b3360.naik.control");
    pass *= testSmear(GRID,Umu,Usmr,Unaik,Ucontrol,0.,0.8675309,0.,0.,0.,0.);

    if(pass){
        Grid_pass("All tests passed.");
    } else {
        Grid_error("At least one test failed.");
    }

#endif 

    // Does a small hot start cause an issue?
    hotStartSmear(GRID);
    latt_size[0]=16; latt_size[1]=16; latt_size[2]=16; latt_size[3]=16;
    GridCartesian SYMM(latt_size,simd_layout,mpi_layout);
    hotStartSmear(SYMM);

    // Test a C-style instantiation 
    double path_coeff[6] = {1, 2, 3, 4, 5, 6};
    Smear_HISQ<GIMPL> hisq_fat_Cstyle(&GRID,path_coeff);

    Grid_finalize();
}