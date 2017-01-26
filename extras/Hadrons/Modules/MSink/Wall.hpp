/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: extras/Hadrons/Modules/MSink/Wall.hpp

Copyright (C) 2016

Author: Andrew Lawson <andrew.lawson1991@gmail.com>

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

#ifndef Hadrons_Wall_hpp_
#define Hadrons_Wall_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/*
 
 Wall sink smearing
 -----------------------------
 * prop_x_3 = sum_x_i (prop_x * exp(-i x.mom))
 
 * options:
 - q: input propagator
 - mom: momentum insertion, space-separated float sequence (e.g ".1 .2 1. 0.")
 
 */

/******************************************************************************
 *                         Wall                                               *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSink)

class WallPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(WallPar,
                                    std::string, q,
                                    std::string, mom);
};

template <typename FImpl>
class TWall: public Module<WallPar>
{
public:
    TYPE_ALIASES(FImpl,);
public:
    // constructor
    TWall(const std::string name);
    // destructor
    virtual ~TWall(void) = default;
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_NS(Wall, TWall<FIMPL>, MSink);

/******************************************************************************
 *                 TWall implementation                                       *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TWall<FImpl>::TWall(const std::string name)
: Module<WallPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TWall<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().q};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TWall<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TWall<FImpl>::setup(void)
{
    
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TWall<FImpl>::execute(void)
{
    LOG(Message) << "Wall smearing " << par().q << std::endl;
    
    PropagatorField &q = *env().template getObject<PropagatorField>(par().q);
    std::vector<typename SitePropagator::scalar_object> prop;
    LatticeComplex              ph(env().getGrid()), coor(env().getGrid());
    std::vector<Real>           p;
    Complex                     i(0.0,1.0);
    
    p  = strToVec<Real>(par().mom);
    ph = zero;
    for(unsigned int mu = 0; mu < Nd; mu++)
    {
        LatticeCoordinate(coor, mu);
        ph = ph + p[mu]*coor;
    }
    ph = exp(-i*ph);
    sliceSum<SitePropagator>(ph*q, prop, Tp);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_Wall_hpp_
