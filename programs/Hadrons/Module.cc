/*******************************************************************************
Grid physics library, www.github.com/paboyle/Grid 

Source file: programs/Hadrons/Module.cc

Copyright (C) 2015

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

#include <Hadrons/Module.hpp>

using namespace Grid;
using namespace QCD;
using namespace Hadrons;

/******************************************************************************
 *                           Module implementation                            *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
Module::Module(const std::string name)
: name_(name)
, env_(Environment::getInstance())
{}

// access //////////////////////////////////////////////////////////////////////
std::string Module::getName(void) const
{
    return name_;
}

Environment & Module::env(void) const
{
    return env_;
}

// execution ///////////////////////////////////////////////////////////////////
void Module::operator()(void)
{
    setup();
    if (!env().isDryRun())
    {
        execute();
    }
}
