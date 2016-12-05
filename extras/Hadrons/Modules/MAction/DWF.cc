/*******************************************************************************
Grid physics library, www.github.com/paboyle/Grid 

Source file: programs/Hadrons/DWF.cc

Copyright (C) 2016

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

See the full license in the file "LICENSE" in the top level distribution 
directory.
*******************************************************************************/

#include <Grid/Hadrons/Modules/MAction/DWF.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MAction;

/******************************************************************************
*                          DWF implementation                                *
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
DWF::DWF(const std::string name)
: Module<DWFPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
std::vector<std::string> DWF::getInput(void)
{
    std::vector<std::string> in = {par().gauge};
    
    return in;
}

std::vector<std::string> DWF::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
void DWF::setup(void)
{
    unsigned int size;
    
    size = 3*env().lattice4dSize<WilsonFermionR::DoubledGaugeField>();
    env().registerObject(getName(), size, par().Ls);
}

// execution ///////////////////////////////////////////////////////////////////
void DWF::execute(void)
{
    LOG(Message) << "Setting up domain wall fermion matrix with m= "
                 << par().mass << ", M5= " << par().M5 << " and Ls= "
                 << par().Ls << " using gauge field '" << par().gauge << "'"
                 << std::endl;
    env().createGrid(par().Ls);
    auto &U      = *env().getObject<LatticeGaugeField>(par().gauge);
    auto &g4     = *env().getGrid();
    auto &grb4   = *env().getRbGrid();
    auto &g5     = *env().getGrid(par().Ls);
    auto &grb5   = *env().getRbGrid(par().Ls);
    FMat *fMatPt = new DomainWallFermion<FIMPL>(U, g5, grb5, g4, grb4,
                                                par().mass, par().M5);
    env().setObject(getName(), fMatPt);
}
