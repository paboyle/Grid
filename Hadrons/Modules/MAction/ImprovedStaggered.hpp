/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MAction/ImprovedStaggered.hpp

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

#ifndef Hadrons_MAction_ImprovedStaggered_hpp_
#define Hadrons_MAction_ImprovedStaggered_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                            TImprovedStaggered quark action                            *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MAction)

class ImprovedStaggeredPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(ImprovedStaggeredPar,
                                    std::string, gauge,
                                    std::string, gaugefat,
                                    std::string, gaugelong,
                                    double     , mass,
                                    double     , c1,
                                    double     , c2,
                                    double     , tad,
                                    std::string, boundary,
                                    std::string, string,
                                    std::string, twist);
};

template <typename FImpl>
class TImprovedStaggered: public Module<ImprovedStaggeredPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TImprovedStaggered(const std::string name);
    // destructor
    virtual ~TImprovedStaggered(void) {};
    // dependencies/products
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(ImprovedStaggered, TImprovedStaggered<STAGIMPL>, MAction);
#ifdef GRID_DEFAULT_PRECISION_DOUBLE
MODULE_REGISTER_TMP(ImprovedStaggeredF, TImprovedStaggered<STAGIMPLF>, MAction);
#endif

/******************************************************************************
 *                     TImprovedStaggered template implementation                        *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TImprovedStaggered<FImpl>::TImprovedStaggered(const std::string name)
: Module<ImprovedStaggeredPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TImprovedStaggered<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().gauge, par().gaugefat, par().gaugelong };
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TImprovedStaggered<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TImprovedStaggered<FImpl>::setup(void)
{
    LOG(Message) << "Setting up ImprovedStaggered fermion matrix." << std::endl;
    LOG(Message) << "Using m=" << par().mass << std::endl;
    LOG(Message) << "Using c1=" << par().c1 << std::endl;
    LOG(Message) << "Using c2=" << par().c2 << std::endl;
    LOG(Message) << "Using tadpole u0=" << par().tad << std::endl;
    LOG(Message) << "Using thin links: " << par().gauge << std::endl;
    LOG(Message) << "Using fat links: " << par().gaugefat << std::endl;
    LOG(Message) << "Using long links: " << par().gaugelong << std::endl;
                 
    auto &U      = envGet(GaugeField, par().gauge);
    auto &Ufat   = envGet(GaugeField, par().gaugefat);
    auto &Ulong  = envGet(GaugeField, par().gaugelong);
    auto &grid   = *envGetGrid(FermionField);
    auto &gridRb = *envGetRbGrid(FermionField);
    typename ImprovedStaggeredFermion<FImpl>::ImplParams implParams;
    
    envCreateDerived(FMat, ImprovedStaggeredFermion<FImpl>, getName(), 1,
                     U, Ufat, Ulong,
                     grid, gridRb,
                     par().mass, par().c1, par().c2, par().tad, implParams);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TImprovedStaggered<FImpl>::execute()
{}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_ImprovedStaggered_hpp_
