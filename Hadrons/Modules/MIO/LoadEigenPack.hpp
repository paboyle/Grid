/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MIO/LoadEigenPack.hpp

Copyright (C) 2015-2019

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
#ifndef Hadrons_MIO_LoadEigenPack_hpp_
#define Hadrons_MIO_LoadEigenPack_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/EigenPack.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                   Load eigen vectors/values package                        *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MIO)

class LoadEigenPackPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(LoadEigenPackPar,
                                    std::string, filestem,
                                    bool, multiFile,
                                    unsigned int, size,
                                    unsigned int, Ls);
};

template <typename Pack>
class TLoadEigenPack: public Module<LoadEigenPackPar>
{
public:
    typedef typename Pack::Field   Field;
    typedef typename Pack::FieldIo FieldIo;
    typedef BaseEigenPack<Field>   BasePack;
public:
    // constructor
    TLoadEigenPack(const std::string name);
    // destructor
    virtual ~TLoadEigenPack(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(LoadFermionEigenPack, TLoadEigenPack<FermionEigenPack<FIMPL>>, MIO);
#ifdef GRID_DEFAULT_PRECISION_DOUBLE
MODULE_REGISTER_TMP(LoadFermionEigenPackIo32, ARG(TLoadEigenPack<FermionEigenPack<FIMPL, FIMPLF>>), MIO);
#endif

/******************************************************************************
 *                    TLoadEigenPack implementation                           *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename Pack>
TLoadEigenPack<Pack>::TLoadEigenPack(const std::string name)
: Module<LoadEigenPackPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename Pack>
std::vector<std::string> TLoadEigenPack<Pack>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename Pack>
std::vector<std::string> TLoadEigenPack<Pack>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename Pack>
void TLoadEigenPack<Pack>::setup(void)
{
    GridBase *gridIo = nullptr;

    if (typeHash<Field>() != typeHash<FieldIo>())
    {
        gridIo = envGetRbGrid(FieldIo, par().Ls);
    }
    envCreateDerived(BasePack, Pack, getName(), par().Ls, par().size, 
                     envGetRbGrid(Field, par().Ls), gridIo);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename Pack>
void TLoadEigenPack<Pack>::execute(void)
{
    auto &epack = envGetDerived(BasePack, Pack, getName());

    epack.read(par().filestem, par().multiFile, vm().getTrajectory());
    epack.eval.resize(par().size);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MIO_LoadEigenPack_hpp_
