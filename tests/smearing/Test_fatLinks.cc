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


/*!  @brief make the logger work like python print */
template<typename... Args>
inline std::string sjoin(Args&&... args) noexcept {
    std::ostringstream msg;
    (msg << ... << args);
    return msg.str();
}
template <typename... Args>
inline void Grid_log(Args&&... args) {
    std::string msg = sjoin(std::forward<Args>(args)...);
    std::cout << GridLogMessage << msg << std::endl;
}

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
    std::string conf_out = "nersc.l8t4b3360.3link";

    // Initialize the Grid
    Grid_init(&argc,&argv);
    Coordinate simd_layout = GridDefaultSimd(Nd,vComplexD::Nsimd());
    Coordinate mpi_layout  = GridDefaultMpi();
    Grid_log(" mpi = ",mpi_layout);
    Grid_log("simd = ",simd_layout);
    Grid_log("latt = ",latt_size);
    GridCartesian GRID(latt_size,simd_layout,mpi_layout);

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
    NerscIO::readConfiguration(Umu, header, "nersc.l8t4b3360.3link.control");
    LatticeGaugeField diff(&GRID);
    diff = Umu-U_smr;
    auto absDiff = norm2(diff)/norm2(Umu);
    Grid_log(" |Umu-U|/|Umu| = ",absDiff);

    Grid_finalize();
}