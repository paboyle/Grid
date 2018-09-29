/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MUtilities/RandomVectors.hpp

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
#ifndef Hadrons_MUtilities_RandomVectors_hpp_
#define Hadrons_MUtilities_RandomVectors_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *           Module generating random lattices for testing purposes           *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)

class RandomVectorsPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(RandomVectorsPar,
                                    unsigned int, size,
                                    unsigned int, Ls);
};

template <typename Field>
class TRandomVectors: public Module<RandomVectorsPar>
{
public:
    // constructor
    TRandomVectors(const std::string name);
    // destructor
    virtual ~TRandomVectors(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(RandomFermions, TRandomVectors<FIMPL::FermionField>, MUtilities);

/******************************************************************************
 *                      TRandomVectors implementation                         *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename Field>
TRandomVectors<Field>::TRandomVectors(const std::string name)
: Module<RandomVectorsPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename Field>
std::vector<std::string> TRandomVectors<Field>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename Field>
std::vector<std::string> TRandomVectors<Field>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename Field>
void TRandomVectors<Field>::setup(void)
{
    if (par().Ls > 1)
    {
        envCreate(std::vector<Field>, getName(), par().Ls, par().size, 
                  envGetGrid(Field, par().Ls));
    }
    else
    {
        envCreate(std::vector<Field>, getName(), 1, par().size, envGetGrid(Field));
    }
}

// execution ///////////////////////////////////////////////////////////////////
template <typename Field>
void TRandomVectors<Field>::execute(void)
{
    LOG(Message) << "Generating " << par().size << " random vectors" << std::endl;

    auto &vec = envGet(std::vector<Field>, getName());
    
    for (unsigned int i = 0; i < vec.size(); ++i)
    {
        random(rng4d(), vec[i]);
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MUtilities_RandomVectors_hpp_
