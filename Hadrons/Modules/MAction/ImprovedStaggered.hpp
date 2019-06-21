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
                                    //std::string, gauge,
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

MODULE_REGISTER_TMP(ImprovedStaggered, TImprovedStaggered<StagIMPL>, MAction);
#ifdef GRID_DEFAULT_PRECISION_DOUBLE
MODULE_REGISTER_TMP(ImprovedStaggeredF, TImprovedStaggered<StagIMPLF>, MAction);
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
    std::vector<std::string> in = {par().gauge};
    
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
    LOG(Message) << "Setting up ImprovedStaggered fermion matrix with m= " << par().mass
                 << " using gauge field '" << par().gauge << "'" << std::endl;
                 
    auto &U      = envGet(GaugeField, par().gauge);
    auto &grid   = *envGetGrid(FermionField);
    auto &gridRb = *envGetRbGrid(FermionField);
    typename ImprovedStaggeredFermion<FImpl>::ImplParams implParams;
    //if (!par().boundary.empty())
    //{
    //    implParams.boundary_phases = strToVec<Complex>(par().boundary);
    //}
    //if (!par().twist.empty())
    //{
    //    implParams.twist_n_2pi_L   = strToVec<Real>(par().twist);
    //}
    //LOG(Message) << "Fermion boundary conditions: " << implParams.boundary_phases
    //             << std::endl;
    //LOG(Message) << "Twists: " << implParams.twist_n_2pi_L
    //             << std::endl;
    //if (implParams.boundary_phases.size() != env().getNd())
    //{
    //    HADRONS_ERROR(Size, "Wrong number of boundary phase");
    //}
    //if (implParams.twist_n_2pi_L.size() != env().getNd())
    //{
    //    HADRONS_ERROR(Size, "Wrong number of twist");
    //}
    envCreateDerived(FMat, ImprovedStaggeredFermion<FImpl>, getName(), 1, U, U, grid, gridRb,
                     par().mass, par().c1, par().c2, par().tad, implParams);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TImprovedStaggered<FImpl>::execute()
{}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_ImprovedStaggered_hpp_