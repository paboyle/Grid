/*******************************************************************************
Grid physics library, www.github.com/paboyle/Grid 

Source file: programs/Hadrons/GUnit.cc

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

#include <Hadrons/GUnit.hpp>

using namespace Grid;
using namespace Hadrons;

/******************************************************************************
*                            GUnit implementation                             *
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
GUnit::GUnit(const std::string name)
: Module(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
std::vector<std::string> GUnit::getInput(void)
{
    return std::vector<std::string>();
}

std::vector<std::string> GUnit::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
void GUnit::setup(void)
{
    env().registerLattice<LatticeGaugeField>(getName());
}

// execution ///////////////////////////////////////////////////////////////////
void GUnit::execute(void)
{
    LOG(Message) << "Creating unit gauge configuration" << std::endl;
    LatticeGaugeField &U = *env().create<LatticeGaugeField>(getName());
    SU3::ColdConfiguration(*env().get4dRng(), U);
}
