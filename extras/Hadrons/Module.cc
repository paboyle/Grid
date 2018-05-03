/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: extras/Hadrons/Module.cc

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

#include <Grid/Hadrons/Module.hpp>

using namespace Grid;
using namespace QCD;
using namespace Hadrons;

/******************************************************************************
 *                       ModuleBase implementation                            *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
ModuleBase::ModuleBase(const std::string name)
: name_(name)
{}

// access //////////////////////////////////////////////////////////////////////
std::string ModuleBase::getName(void) const
{
    return name_;
}

// get factory registration name if available
std::string ModuleBase::getRegisteredName(void)
{
    HADRONS_ERROR(Definition, "module '" + getName() + "' has no registered type"
                 + " in the factory");
}

// execution ///////////////////////////////////////////////////////////////////
void ModuleBase::operator()(void)
{
    setup();
    execute();
}
