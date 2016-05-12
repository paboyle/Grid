/*******************************************************************************
Grid physics library, www.github.com/paboyle/Grid 

Source file: programs/Hadrons/SrcPoint.cc

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

#include <Hadrons/Modules/SrcPoint.hpp>

using namespace Grid;
using namespace Hadrons;

/******************************************************************************
*                         SrcPoint implementation                             *
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
SrcPoint::SrcPoint(const std::string name)
: Module<SrcPointPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
std::vector<std::string> SrcPoint::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

std::vector<std::string> SrcPoint::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
void SrcPoint::setup(void)
{
    env().registerLattice<LatticePropagator>(getName());
}

// execution ///////////////////////////////////////////////////////////////////
void SrcPoint::execute(void)
{
    std::vector<int> position = strToVec<int>(par().position);
    SpinColourMatrix id;
    
    LOG(Message) << "Creating point source at position [" << par().position
                 << "]" << std::endl;
    LatticePropagator &src = *env().create<LatticePropagator>(getName());
    id  = 1.;
    src = zero;
    pokeSite(id, src, position);
}
