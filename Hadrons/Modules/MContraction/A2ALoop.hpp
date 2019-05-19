/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MContraction/A2ALoop.hpp

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
#ifndef Hadrons_MContraction_A2ALoop_hpp_
#define Hadrons_MContraction_A2ALoop_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                  From closed loop from all-to-all vectors                  *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class A2ALoopPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2ALoopPar,
                                    std::string, left,
                                    std::string, right);
};

template <typename FImpl>
class TA2ALoop: public Module<A2ALoopPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TA2ALoop(const std::string name);
    // destructor
    virtual ~TA2ALoop(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(A2ALoop, TA2ALoop<FIMPL>, MContraction);

/******************************************************************************
 *                        TA2ALoop implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TA2ALoop<FImpl>::TA2ALoop(const std::string name)
: Module<A2ALoopPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TA2ALoop<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().left, par().right};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TA2ALoop<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2ALoop<FImpl>::setup(void)
{
    envCreateLat(PropagatorField, getName());
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2ALoop<FImpl>::execute(void)
{
    auto &loop  = envGet(PropagatorField, getName());
    auto &left  = envGet(std::vector<FermionField>, par().left);
    auto &right = envGet(std::vector<FermionField>, par().right);

    loop = zero;
    for (unsigned int i = 0; i < left.size(); ++i)
    {
        loop += outerProduct(adj(left[i]), right[i]);
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_A2ALoop_hpp_
