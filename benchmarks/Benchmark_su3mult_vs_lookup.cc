    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./benchmarks/Benchmark_su3mult_vs_lookup.cc

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

    See the full license in the file "LICENSE" in the top level distribution directory
    *************************************************************************************/
/*
    @file Benchmark_su3mult_vs_lookup.cc
    @brief check to see whether su3 multiplication or lookup tables is faster 
*/

#include <Grid/Grid.h>
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

/*!  @brief parameter file to easily adjust Nloop */
struct ConfParameters: Serializable {
    GRID_SERIALIZABLE_CLASS_MEMBERS(
        ConfParameters,
        int, Nloop);

    template <class ReaderClass>
    ConfParameters(Reader<ReaderClass>& Reader){
        read(Reader, "parameters", *this);
    }
};

int main (int argc, char** argv) {

    // Params for the test.
    int Ns               = 8;
    int Nt               = 4;
    int threads          = GridThread::GetThreads();
    std::string conf_in  = "nersc.l8t4b3360";
    Coordinate latt_size(Nd,0); latt_size[0]=Ns; latt_size[1]=Ns; latt_size[2]=Ns; latt_size[3]=Nt;

    Grid_init(&argc,&argv);

    Coordinate simd_layout = GridDefaultSimd(Nd,vComplexD::Nsimd());
    Coordinate mpi_layout  = GridDefaultMpi();
    GridCartesian GRID(latt_size,simd_layout,mpi_layout);

    Grid_log("    mpi = ",mpi_layout);
    Grid_log("   simd = ",simd_layout);
    Grid_log("   latt = ",latt_size);
    Grid_log("threads = ",threads);

    XmlReader Reader("mult_vs_lookup.xml",false, "grid");
    ConfParameters param(Reader); 
    Grid_log("  Nloop = ",param.Nloop);

    // Gauge field and accessor
    LatticeGaugeField Umu(&GRID);
    autoView(U_v, Umu, CpuRead);

    // Read the configuration into Umu
    FieldMetaData header;
    NerscIO::readConfiguration(Umu, header, conf_in);

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
        Grid_log("Time to multiply ",nmat," matrices: ",lookupTime," [ms]");
    }

    Grid_finalize();
}
