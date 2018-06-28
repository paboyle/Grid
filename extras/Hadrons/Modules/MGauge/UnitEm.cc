/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: extras/Hadrons/Modules/MGauge/StochEm.cc

Copyright (C) 2015
Copyright (C) 2016

Author: James Harrison <j.harrison@soton.ac.uk>

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
#include <Grid/Hadrons/Modules/MGauge/UnitEm.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MGauge;

/******************************************************************************
*                  TStochEm implementation                             *
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
TUnitEm::TUnitEm(const std::string name)
: Module<NoPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
std::vector<std::string> TUnitEm::getInput(void)
{
    return std::vector<std::string>();
}

std::vector<std::string> TUnitEm::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
void TUnitEm::setup(void)
{
    envCreateLat(EmField, getName());
}

// execution ///////////////////////////////////////////////////////////////////
void TUnitEm::execute(void)
{
    PhotonR photon(0, 0); // Just chose arbitrary input values here
    auto    &a = envGet(EmField, getName());
    LOG(Message) << "Generating unit EM potential..." << std::endl;
    photon.UnitField(a);
}
