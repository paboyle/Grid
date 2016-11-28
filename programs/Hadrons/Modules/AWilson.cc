/*******************************************************************************
Grid physics library, www.github.com/paboyle/Grid 

Source file: programs/Hadrons/AWilson.cc

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

#include <Grid/Hadrons/Modules/AWilson.hpp>

using namespace Grid;
using namespace Hadrons;

/******************************************************************************
*                         AWilson implementation                              *
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
AWilson::AWilson(const std::string name)
: Module<AWilsonPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
std::vector<std::string> AWilson::getInput(void)
{
    std::vector<std::string> in = {par().gauge};
    
    return in;
}

std::vector<std::string> AWilson::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
void AWilson::setup(void)
{
    unsigned int size;
    
    size = 3*env().lattice4dSize<WilsonFermionR::DoubledGaugeField>();
    env().registerObject(getName(), size);
}

// execution ///////////////////////////////////////////////////////////////////
void AWilson::execute()
{
    LOG(Message) << "Setting up Wilson fermion matrix with m= " << par().mass
                 << " using gauge field '" << par().gauge << "'" << std::endl;
    auto &U      = *env().getObject<LatticeGaugeField>(par().gauge);
    auto &grid   = *env().getGrid();
    auto &gridRb = *env().getRbGrid();
    FMat *fMatPt = new WilsonFermionR(U, grid, gridRb, par().mass);
    env().setObject(getName(), fMatPt);
}
