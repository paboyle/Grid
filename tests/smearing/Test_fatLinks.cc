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


/*!  @brief parameter file to easily adjust Nloop */
struct ConfParameters: Serializable {
    GRID_SERIALIZABLE_CLASS_MEMBERS(
        ConfParameters,
        int, benchmark, 
        int, Nloop);

    template <class ReaderClass>
    ConfParameters(Reader<ReaderClass>& Reader){
        read(Reader, "parameters", *this);
    }
};

//
// one method: input --> fat
// another   : input --> long (naik)
// another   : input --> unitarize
//

int main (int argc, char** argv) {

    // Params for the test.
    int Ns = 8;
    int Nt = 4;
    Coordinate latt_size(Nd,0); latt_size[0]=Ns; latt_size[1]=Ns; latt_size[2]=Ns; latt_size[3]=Nt;
    std::string conf_in  = "nersc.l8t4b3360";
    std::string conf_out = "nersc.l8t4b3360.357link";
    int threads          = GridThread::GetThreads();

    // Initialize the Grid
    Grid_init(&argc,&argv);
    Coordinate simd_layout = GridDefaultSimd(Nd,vComplexD::Nsimd());
    Coordinate mpi_layout  = GridDefaultMpi();
    Grid_log("mpi     = ",mpi_layout);
    Grid_log("simd    = ",simd_layout);
    Grid_log("latt    = ",latt_size);
    Grid_log("threads = ",threads);
    GridCartesian GRID(latt_size,simd_layout,mpi_layout);

    XmlReader Reader("fatParams.xml",false,"grid");
    ConfParameters param(Reader);
    if(param.benchmark) Grid_log("  Nloop = ",param.Nloop);

    // Instantiate the LatticeGaugeField objects holding thin (Umu) and fat (U_smr) links
    LatticeGaugeField Umu(&GRID);
    LatticeGaugeField U_smr(&GRID);

    // Read the configuration into Umu
    FieldMetaData header;
    NerscIO::readConfiguration(Umu, header, conf_in);

    // Smear Umu and store result in U_smr
    Smear_HISQ_fat<LatticeGaugeField> hisq_fat(&GRID,1/8.,0.,1/16.,1/64.,1/384.,0.);
    hisq_fat.smear(U_smr,Umu);

    NerscIO::writeConfiguration(U_smr,conf_out,"HISQ");

    // Test a C-style instantiation 
    double path_coeff[6] = {1, 2, 3, 4, 5, 6};
    Smear_HISQ_fat<LatticeGaugeField> hisq_fat_Cstyle(&GRID,path_coeff);

    // Make sure result doesn't change w.r.t. a trusted lattice
    NerscIO::readConfiguration(Umu, header, "nersc.l8t4b3360.357link.control");
    LatticeGaugeField diff(&GRID);
    diff = Umu-U_smr;
    auto absDiff = norm2(diff)/norm2(Umu);
    Grid_log(" |Umu-U|/|Umu| = ",absDiff);


    if (param.benchmark) {

        autoView(U_v, Umu, CpuRead); // Gauge accessor

        // Read in lattice sequentially, Nloop times 
        double lookupTime = 0.; 
        for(int i=0;i<param.Nloop;i++) {
            double start = usecond();
            for(int ss=0;ss<U_v.size();ss++)
                for(int mu=0;mu<Nd;mu++) {
                    auto U1 = U_v[ss](mu);
            }
            double stop  = usecond();
        	lookupTime += stop-start; // microseconds
        }
        Grid_log("Time to lookup: ",lookupTime,"[ms]");

        // Raise a matrix to the power nmat, for each link. 
        auto U1 = U_v[0](0);
        for(int nmat=1;nmat<8;nmat++) {
            double multTime = 0.; 
            for(int i=0;i<param.Nloop;i++) {
                double start=usecond();
                for(int ss=0;ss<U_v.size();ss++)
                    for(int mu=0;mu<Nd;mu++) {
                        auto U2 = U1;
                        for(int j=1;j<nmat;j++) {
                            U2 *= U1;
                        }
                }
                double stop=usecond();
                multTime += stop-start;
            }
            Grid_log("Time to multiply ",nmat," matrices: ",multTime," [ms]");
        }
    }

    Grid_finalize();
}