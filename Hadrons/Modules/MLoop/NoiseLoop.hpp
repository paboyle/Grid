/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MLoop/NoiseLoop.hpp

Copyright (C) 2015-2018

Author: Antonin Portelli <antonin.portelli@me.com>
Author: Lanny91 <andrew.lawson@gmail.com>

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

#ifndef Hadrons_MLoop_NoiseLoop_hpp_
#define Hadrons_MLoop_NoiseLoop_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/*
 
 Noise loop propagator
 -----------------------------
 * loop_x = q_x * adj(eta_x)
 
 * options:
 - q = Result of inversion on noise source.
 - eta = noise source.

 */


/******************************************************************************
 *                         NoiseLoop                                          *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MLoop)

class NoiseLoopPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(NoiseLoopPar,
                                    std::string, q,
                                    std::string, eta);
};

template <typename FImpl>
class TNoiseLoop: public Module<NoiseLoopPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TNoiseLoop(const std::string name);
    // destructor
    virtual ~TNoiseLoop(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(NoiseLoop, TNoiseLoop<FIMPL>, MLoop);

/******************************************************************************
 *                 TNoiseLoop implementation                                  *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TNoiseLoop<FImpl>::TNoiseLoop(const std::string name)
: Module<NoiseLoopPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TNoiseLoop<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().q, par().eta};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TNoiseLoop<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TNoiseLoop<FImpl>::setup(void)
{
    envCreateLat(PropagatorField, getName());
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TNoiseLoop<FImpl>::execute(void)
{
    auto &loop = envGet(PropagatorField, getName());
    auto &q    = envGet(PropagatorField, par().q);
    auto &eta  = envGet(PropagatorField, par().eta);
    loop = q*adj(eta);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MLoop_NoiseLoop_hpp_
