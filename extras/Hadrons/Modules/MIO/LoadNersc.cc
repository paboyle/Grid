/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: extras/Hadrons/Modules/MIO/LoadNersc.cc

Copyright (C) 2015-2018

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
#include <Grid/Hadrons/Modules/MIO/LoadNersc.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MIO;

/******************************************************************************
*                       TLoadNersc implementation                             *
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
TLoadNersc::TLoadNersc(const std::string name)
: Module<LoadNerscPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
std::vector<std::string> TLoadNersc::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

std::vector<std::string> TLoadNersc::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
void TLoadNersc::setup(void)
{
    envCreateLat(LatticeGaugeField, getName());
}

// execution ///////////////////////////////////////////////////////////////////
void TLoadNersc::execute(void)
{
    FieldMetaData header;
    std::string   fileName = par().file + "."
                             + std::to_string(vm().getTrajectory());
    LOG(Message) << "Loading NERSC configuration from file '" << fileName
                 << "'" << std::endl;

    auto &U = envGet(LatticeGaugeField, getName());
    NerscIO::readConfiguration(U, header, fileName);
}
