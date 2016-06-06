/*******************************************************************************
Grid physics library, www.github.com/paboyle/Grid 

Source file: programs/Hadrons/GLoad.cc

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

#include <Hadrons/Modules/GLoad.hpp>

using namespace Grid;
using namespace Hadrons;

/******************************************************************************
*                          GLoad implementation                               *
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
GLoad::GLoad(const std::string name)
: Module<GLoadPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
std::vector<std::string> GLoad::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

std::vector<std::string> GLoad::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
void GLoad::setup(void)
{
    env().registerLattice<LatticeGaugeField>(getName());
}

// execution ///////////////////////////////////////////////////////////////////
void GLoad::execute(void)
{
    NerscField  header;
    std::string fileName = par().file + "."
                           + std::to_string(env().getTrajectory());
    
    LOG(Message) << "Loading NERSC configuration from file '" << fileName
                 << "'" << std::endl;
    LatticeGaugeField &U = *env().createLattice<LatticeGaugeField>(getName());
    NerscIO::readConfiguration(U, header, fileName);
    LOG(Message) << "NERSC header:" << std::endl;
    dump_nersc_header(header, LOG(Message));
}
