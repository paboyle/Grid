/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: extras/Hadrons/Modules/MAction/DWF.hpp

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

#ifndef Hadrons_MAction_DWF_hpp_
#define Hadrons_MAction_DWF_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                     Domain wall quark action                               *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MAction)

class DWFPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(DWFPar,
                                    std::string, gauge,
                                    unsigned int, Ls,
                                    double      , mass,
                                    double      , M5,
                                    std::string , boundary);
};

template <typename FImpl>
class TDWF: public Module<DWFPar>
{
public:
    FG_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TDWF(const std::string name);
    // destructor
    virtual ~TDWF(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

extern template class TDWF<FIMPL>;
MODULE_REGISTER_TMP(DWF, TDWF<FIMPL>, MAction);

/******************************************************************************
 *                        DWF template implementation                         *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TDWF<FImpl>::TDWF(const std::string name)
: Module<DWFPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TDWF<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().gauge};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TDWF<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDWF<FImpl>::setup(void)
{
    LOG(Message) << "Setting up domain wall fermion matrix with m= "
                 << par().mass << ", M5= " << par().M5 << " and Ls= "
                 << par().Ls << " using gauge field '" << par().gauge << "'"
                 << std::endl;
    LOG(Message) << "Fermion boundary conditions: " << par().boundary
                 << std::endl;
                 
    env().createGrid(par().Ls);
    auto &U    = envGet(LatticeGaugeField, par().gauge);
    auto &g4   = *env().getGrid();
    auto &grb4 = *env().getRbGrid();
    auto &g5   = *env().getGrid(par().Ls);
    auto &grb5 = *env().getRbGrid(par().Ls);
    std::vector<Complex> boundary = strToVec<Complex>(par().boundary);
    typename DomainWallFermion<FImpl>::ImplParams implParams(boundary);
    envCreateDerived(FMat, DomainWallFermion<FImpl>, getName(), par().Ls, U, g5,
                     grb5, g4, grb4, par().mass, par().M5, implParams);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDWF<FImpl>::execute(void)
{}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MAction_DWF_hpp_
