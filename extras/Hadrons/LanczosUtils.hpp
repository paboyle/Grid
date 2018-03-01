/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: extras/Hadrons/LanczosUtils.hpp

Copyright (C) 2015-2018

Author: Antonin Portelli <antonin.portelli@me.com>

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
#ifndef Hadrons_LanczosUtils_hpp_
#define Hadrons_LanczosUtils_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/algorithms/iterative/LocalCoherenceLanczos.h>

BEGIN_HADRONS_NAMESPACE

// Lanczos type
#ifndef HADRONS_DEFAULT_LANCZOS_NBASIS
#define HADRONS_DEFAULT_LANCZOS_NBASIS 60
#endif

template <typename T>
struct EigenPack
{
    typedef T VectorType;
    std::vector<RealD> eval;
    std::vector<T>     evec;
    
    EigenPack(void) = default;

    EigenPack(const size_t size, GridBase *grid)
    {
        resize(size, grid);
    }

    void resize(const size_t size, GridBase *grid)
    {
        eval.resize(size);
        evec.resize(size, grid);
    }

    void read(const std::string fileStem)
    {
        std::string     evecFilename = fileStem + "_evec.bin";
        std::string     evalFilename = fileStem + "_eval.xml";
        emptyUserRecord record;
        ScidacReader    binReader;
        XmlReader       xmlReader(evalFilename);

        LOG(Message) << "Reading " << evec.size() << " eigenvectors from '" 
                     << evecFilename << "'" << std::endl;
        binReader.open(evecFilename);
        for(int k = 0; k < evec.size(); ++k) 
        {
            binReader.readScidacFieldRecord(evec[k], record);
        }
        binReader.close();
        LOG(Message) << "Reading " << eval.size() << " eigenvalues from '" 
                     << evalFilename << "'" << std::endl;
        Grid::read(xmlReader, "evals", eval);
    }

    void write(const std::string fileStem)
    {
        std::string     evecFilename = fileStem + "_evec.bin";
        std::string     evalFilename = fileStem + "_eval.xml";
        emptyUserRecord record;
        ScidacWriter    binWriter;
        XmlWriter       xmlWriter(evalFilename);

        LOG(Message) << "Writing " << evec.size() << " eigenvectors to '" 
                     << evecFilename << "'" << std::endl;
        binWriter.open(fileStem + "_evec.bin");
        for(int k = 0; k < evec.size(); ++k) 
        {
            binWriter.writeScidacFieldRecord(evec[k], record);
        }
        binWriter.close();
        LOG(Message) << "Writing " << eval.size() << " eigenvalues to '" 
                     << evalFilename << "'" << std::endl;
        Grid::write(xmlWriter, "evals", eval);
    }
};

template <typename FImpl>
using FineEigenPack = EigenPack<typename FImpl::FermionField>;

template <typename FImpl, int nBasis>
using CoarseEigenPack = EigenPack<
    typename LocalCoherenceLanczos<typename FImpl::SiteSpinor, 
                                   typename FImpl::SiteComplex, 
                                   nBasis>::CoarseField>;

END_HADRONS_NAMESPACE

#endif // Hadrons_LanczosUtils_hpp_