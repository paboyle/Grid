/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: extras/Hadrons/Modules/MGauge/Load.cc

Copyright (C) 2015
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

See the full license in the file "LICENSE" in the top level distribution directory
*************************************************************************************/
/*  END LEGAL */

#include <Grid/Hadrons/Modules/MGauge/Load.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MGauge;

/******************************************************************************
*                           TLoad implementation                               *
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
TLoad::TLoad(const std::string name)
: Module<LoadPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
std::vector<std::string> TLoad::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

std::vector<std::string> TLoad::getReference(void)
{
    std::vector<std::string> ref;
    
    return ref;
}

std::vector<std::string> TLoad::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
void TLoad::setup(void)
{
    envCreateLat(LatticeGaugeField, getName());
}

// execution ///////////////////////////////////////////////////////////////////
void TLoad::execute(void)
{
    FieldMetaData header;
    std::string   fileName = par().file + "."
                             + std::to_string(vm().getTrajectory());
    LOG(Message) << "Loading NERSC configuration from file '" << fileName
                 << "'" << std::endl;

    auto &U = envGet(LatticeGaugeField, getName());
    NerscIO::readConfiguration(U, header, fileName);
    LOG(Message) << "NERSC header:" << std::endl;
    dump_meta_data(header, LOG(Message));
}
