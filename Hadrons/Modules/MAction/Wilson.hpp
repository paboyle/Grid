/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MAction/Wilson.hpp

Copyright (C) 2015-2019

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

#ifndef Hadrons_MAction_Wilson_hpp_
#define Hadrons_MAction_Wilson_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                            TWilson quark action                            *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MAction)

class WilsonPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(WilsonPar,
                                    std::string, gauge,
                                    double     , mass,
                                    std::string, boundary,
                                    std::string, string,
                                    std::string, twist);
};

template <typename FImpl>
class TWilson: public Module<WilsonPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TWilson(const std::string name);
    // destructor
    virtual ~TWilson(void) {};
    // dependencies/products
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(Wilson, TWilson<FIMPL>, MAction);
#ifdef GRID_DEFAULT_PRECISION_DOUBLE
MODULE_REGISTER_TMP(WilsonF, TWilson<FIMPLF>, MAction);
#endif

/******************************************************************************
 *                     TWilson template implementation                        *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TWilson<FImpl>::TWilson(const std::string name)
: Module<WilsonPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TWilson<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().gauge};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TWilson<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TWilson<FImpl>::setup(void)
{
    LOG(Message) << "Setting up Wilson fermion matrix with m= " << par().mass
                 << " using gauge field '" << par().gauge << "'" << std::endl;
                 
    auto &U      = envGet(GaugeField, par().gauge);
    auto &grid   = *envGetGrid(FermionField);
    auto &gridRb = *envGetRbGrid(FermionField);
    typename WilsonFermion<FImpl>::ImplParams implParams;
    if (!par().boundary.empty())
    {
        implParams.boundary_phases = strToVec<Complex>(par().boundary);
    }
    if (!par().twist.empty())
    {
        implParams.twist_n_2pi_L   = strToVec<Real>(par().twist);
    }
    LOG(Message) << "Fermion boundary conditions: " << implParams.boundary_phases << std::endl;
    LOG(Message) << "Twists: " << implParams.twist_n_2pi_L << std::endl;
    if (implParams.boundary_phases.size() != env().getNd())
    {
        HADRONS_ERROR(Size, "Wrong number of boundary phase");
    }
    if (implParams.twist_n_2pi_L.size() != env().getNd())
    {
        HADRONS_ERROR(Size, "Wrong number of twist");
    }
    envCreateDerived(FMat, WilsonFermion<FImpl>, getName(), 1, U, grid, gridRb,
                     par().mass, implParams);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TWilson<FImpl>::execute()
{}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_Wilson_hpp_
