/*************************************************************************************
 
 Grid physics library, www.github.com/paboyle/Grid
 
 Source file: Hadrons/Modules/MDistil/DistilCommon.hpp
 
 Copyright (C) 2015-2019
 
 Author: Felix Erben <ferben@ed.ac.uk>
 Author: Michael Marshall <Michael.Marshall@ed.ac.uk>
 
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
/*  END LEGAL */

#ifndef Hadrons_MDistil_DistilCommon_hpp_
#define Hadrons_MDistil_DistilCommon_hpp_

#include <Hadrons/Distil.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Solver.hpp>
#include <Hadrons/A2AVectors.hpp>
#include <Hadrons/DilutedNoise.hpp>

BEGIN_HADRONS_NAMESPACE
BEGIN_MODULE_NAMESPACE(MDistil)

/******************************************************************************
 Distillation code that is common across modules
 ******************************************************************************/

struct DistilParameters: Serializable {
    GRID_SERIALIZABLE_CLASS_MEMBERS(DistilParameters,
                                    int, nnoise,
                                    int, tsrc,
                                    std::string, TI,
                                    std::string, LI,
                                    std::string, SI )
    DistilParameters() = default;
    template <class ReaderClass> DistilParameters(Reader<ReaderClass>& Reader){read(Reader,"Distil",*this);}
    
    // Numeric parameter is allowed to be empty (in which case it = Default),
    // but assert during setup() if specified but not numeric
    
    static int ParameterDefault( const std::string & s, int Default, bool bCalledFromSetup )
    {
        int i = Default;
        if( s.length() > 0 ) {
            std::istringstream ss( s );
            ss >> i;
            if( bCalledFromSetup )
                assert( !ss.fail() && "Parameter should either be empty or integer" );
        }
        return i;
    }
};

#define DISTIL_PARAMETERS_DEFINE( inSetup ) \
const int Nt{env().getDim(Tdir)}; \
const int nvec{par().nvec}; \
const int nnoise{par().Distil.nnoise}; \
const int tsrc{par().Distil.tsrc}; \
const int TI{Hadrons::MDistil::DistilParameters::ParameterDefault(par().Distil.TI, Nt,   inSetup)}; \
const int LI{Hadrons::MDistil::DistilParameters::ParameterDefault(par().Distil.LI, nvec, inSetup)}; \
const int SI{Hadrons::MDistil::DistilParameters::ParameterDefault(par().Distil.SI, Ns,   inSetup)}; \
const bool full_tdil{ TI == Nt }; \
const bool exact_distillation{ full_tdil && LI == nvec }; \
const int Nt_inv{ full_tdil ? 1 : TI }

/******************************************************************************
 Make a lower dimensional grid in preparation for local slice operations
 ******************************************************************************/

inline GridCartesian * MakeLowerDimGrid( GridCartesian * gridHD )
{
    int nd{static_cast<int>(gridHD->_ndimension)};
    Coordinate latt_size   = gridHD->_gdimensions;
    latt_size[nd-1] = 1;
    Coordinate simd_layout = GridDefaultSimd(nd-1, vComplex::Nsimd());
    simd_layout.push_back( 1 );
    Coordinate mpi_layout  = gridHD->_processors;
    mpi_layout[nd-1] = 1;
    GridCartesian * gridLD = new GridCartesian(latt_size,simd_layout,mpi_layout,*gridHD);
    return gridLD;
}

/*************************************************************************************
 Rotate eigenvectors into our phase convention
 First component of first eigenvector is real and positive
 TODO: Should this be in Distil.hpp?
 *************************************************************************************/

inline void RotateEigen(std::vector<LatticeColourVector> & evec)
{
    ColourVector cv0;
    auto grid = evec[0].Grid();
    Coordinate siteFirst(grid->Nd(),0);
    peekSite(cv0, evec[0], siteFirst);
    Grid::Complex cplx0 = cv0()()(0);
    if( cplx0.imag() == 0 )
        std::cout << GridLogMessage << "RotateEigen() : Site 0 : " << cplx0 << " => already meets phase convention" << std::endl;
    else {
        const Real cplx0_mag = Grid::abs(cplx0);
#ifdef GRID_NVCC
        const Grid::Complex phase = thrust::conj(cplx0 / cplx0_mag);
        const Real argphase = thrust::arg(phase);
#else
        const Grid::Complex phase = std::conj(cplx0 / cplx0_mag);
        const Real argphase = std::arg(phase);
#endif
        std::cout << GridLogMessage << "RotateEigen() : Site 0 : |" << cplx0 << "|=" << cplx0_mag << " => phase=" << (argphase / 3.14159265) << " pi" << std::endl;
        {
            // TODO: Only really needed on the master slice
            for( int k = 0 ; k < evec.size() ; k++ )
                evec[k] *= phase;
            if(grid->IsBoss()){
                for( int c = 0 ; c < Nc ; c++ )
                    cv0()()(c) *= phase;
                cplx0.imag(0); // This assumes phase convention is real, positive (so I get rid of rounding error)
                //pokeSite(cv0, evec[0], siteFirst);
                pokeLocalSite(cv0, evec[0], siteFirst);
            }
        }
    }
}

END_MODULE_NAMESPACE
END_HADRONS_NAMESPACE
#endif // Hadrons_MDistil_DistilCommon_hpp_
