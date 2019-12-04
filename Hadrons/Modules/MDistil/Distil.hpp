/*************************************************************************************
 
 Grid physics library, www.github.com/paboyle/Grid
 
 Source file: Hadrons/Modules/MDistil/Distil.hpp
 
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

#ifndef Hadrons_MDistil_Distil_hpp_
#define Hadrons_MDistil_Distil_hpp_

#include <Hadrons/NamedTensor.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Solver.hpp>
#include <Hadrons/A2AVectors.hpp>
#include <Hadrons/DilutedNoise.hpp>

BEGIN_HADRONS_NAMESPACE
BEGIN_MODULE_NAMESPACE(MDistil)

/******************************************************************************
 Distillation code that is common across modules

 Documentation on how to use this code available at

 *  https://aportelli.github.io/Hadrons-doc/#/mdistil  *
 
 Notation for (stochastic) DistilParameters taken from 1104.3870:

 TI is interlaced dilution in time (corresponding to Nt = time-dimension of the lattice)
 LI is interlaced dilution in laplacian-eigenvector space (corresponding to nvec)
 SI is interlaced dilution in spin (corresponding to Ns, taken from Grid, usually Ns=4)

 This code automatically computes perambulators using exact distillation if
 *   (TI,LI,SI) = (Nt,nvec,Ns)   *
 In this case, nnoise=1 and Noises is set to an array of values =1 as well.
 tsrc then specifies the only timeslice on which the sources are supported.
 (( for stochastic distillation, the vaue of tsrc has no meaning in this code ))

 ******************************************************************************/

struct DistilParameters: Serializable {
    GRID_SERIALIZABLE_CLASS_MEMBERS(DistilParameters,
                                    int, nvec,
                                    int, nnoise,
                                    int, tsrc,
                                    int, TI,
                                    int, LI,
                                    int, SI )
};

/******************************************************************************
 Make a lower dimensional grid in preparation for local slice operations
 ******************************************************************************/

inline void MakeLowerDimGrid( std::unique_ptr<GridCartesian> &up, GridCartesian * gridHD )
{
    int nd{static_cast<int>(gridHD->_ndimension)};
    Coordinate latt_size   = gridHD->_gdimensions;
    latt_size[nd-1] = 1;
    Coordinate simd_layout = GridDefaultSimd(nd-1, vComplex::Nsimd());
    simd_layout.push_back( 1 );
    Coordinate mpi_layout  = gridHD->_processors;
    mpi_layout[nd-1] = 1;
    up.reset( new GridCartesian(latt_size,simd_layout,mpi_layout,*gridHD) );
}

/*************************************************************************************
 Rotate eigenvectors into our phase convention
 First component of first eigenvector is real and positive
 *************************************************************************************/

inline void RotateEigen(std::vector<LatticeColourVector> & evec)
{
    ColourVector cv0;
    auto grid = evec[0].Grid();
    Coordinate siteFirst(grid->Nd(),0);
    peekSite(cv0, evec[0], siteFirst);
    const std::complex<Real> cplx0{cv0()()(0).real(), cv0()()(0).imag()};
    if( cplx0.imag() == 0 )
        LOG(Message) << "RotateEigen() : Site 0 : " << cplx0 << " => already meets phase convention" << std::endl;
    else
    {
        const Real cplx0_mag{ std::abs(cplx0) };
        const std::complex<Real> std_phase{std::conj(cplx0/cplx0_mag)};
        LOG(Message) << "RotateEigen() : Site 0 : |" << cplx0 << "|=" << cplx0_mag
                     << " => phase=" << (std::arg(std_phase) / M_PI) << " pi" << std::endl;
        {
            const Grid::Complex phase{std_phase.real(),std_phase.imag()};
            for( int k = 0 ; k < evec.size() ; k++ )
                evec[k] *= phase;
            // Get rid of the rounding error in imaginary phase on the very first site
            peekSite(cv0, evec[0], siteFirst);
            cv0()()(0).imag(0); // this should be zero after the phase multiply - force it to be so
            pokeSite(cv0, evec[0], siteFirst);
        }
    }
}

END_MODULE_NAMESPACE
END_HADRONS_NAMESPACE
#endif
