/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MIO/LoadCoarseEigenPack.hpp

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
#ifndef Hadrons_MIO_LoadCoarseEigenPack_hpp_
#define Hadrons_MIO_LoadCoarseEigenPack_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/EigenPack.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *              Load local coherence eigen vectors/values package             *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MIO)

class LoadCoarseEigenPackPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(LoadCoarseEigenPackPar,
                                    std::string, filestem,
                                    bool,         multiFile,
                                    unsigned int, sizeFine,
                                    unsigned int, sizeCoarse,
                                    unsigned int, Ls,
                                    std::vector<int>, blockSize);
};

template <typename Pack>
class TLoadCoarseEigenPack: public Module<LoadCoarseEigenPackPar>
{
public:
    typedef CoarseEigenPack<typename Pack::Field, typename Pack::CoarseField> BasePack;
    template <typename vtype> 
    using iImplScalar = iScalar<iScalar<iScalar<vtype>>>;
    typedef iImplScalar<typename Pack::Field::vector_type> SiteComplex;
public:
    // constructor
    TLoadCoarseEigenPack(const std::string name);
    // destructor
    virtual ~TLoadCoarseEigenPack(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(LoadCoarseFermionEigenPack, ARG(TLoadCoarseEigenPack<CoarseFermionEigenPack<FIMPL, HADRONS_DEFAULT_LANCZOS_NBASIS>>), MIO);

/******************************************************************************
 *                 TLoadCoarseEigenPack implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename Pack>
TLoadCoarseEigenPack<Pack>::TLoadCoarseEigenPack(const std::string name)
: Module<LoadCoarseEigenPackPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename Pack>
std::vector<std::string> TLoadCoarseEigenPack<Pack>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename Pack>
std::vector<std::string> TLoadCoarseEigenPack<Pack>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename Pack>
void TLoadCoarseEigenPack<Pack>::setup(void)
{
    env().createGrid(par().Ls);
    env().createCoarseGrid(par().blockSize, par().Ls);
    envCreateDerived(BasePack, Pack, getName(), par().Ls, par().sizeFine,
                     par().sizeCoarse, env().getRbGrid(par().Ls), 
                     env().getCoarseGrid(par().blockSize, par().Ls));
}

// execution ///////////////////////////////////////////////////////////////////
template <typename Pack>
void TLoadCoarseEigenPack<Pack>::execute(void)
{
    auto                 cg     = env().getCoarseGrid(par().blockSize, par().Ls);
    auto                 &epack = envGetDerived(BasePack, Pack, getName());
    Lattice<SiteComplex> dummy(cg);

    epack.read(par().filestem, par().multiFile, vm().getTrajectory());
    LOG(Message) << "Block Gramm-Schmidt pass 1"<< std::endl;
    blockOrthogonalise(dummy, epack.evec);
    LOG(Message) << "Block Gramm-Schmidt pass 2"<< std::endl;
    blockOrthogonalise(dummy, epack.evec);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MIO_LoadCoarseEigenPack_hpp_
