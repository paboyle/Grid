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

#include <Hadrons/AWilson.hpp>
#include <Hadrons/Environment.hpp>

using namespace Grid;
using namespace QCD;
using namespace Hadrons;

/******************************************************************************
*                          AWilson implementation                             *
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
AWilson::AWilson(const std::string name)
: FermionAction(name)
{}

// parse parameters ////////////////////////////////////////////////////////////
void AWilson::parseParameters(XmlReader &reader, const std::string name)
{
   read(reader, name, par_);
}

// create operator /////////////////////////////////////////////////////////////
void AWilson::create(Environment &env)
{
    auto &U      = *env.getGauge();
    auto &grid   = *env.getGrid();
    auto &gridRb = *env.getRbGrid();
    
    setFMat(new WilsonFermionR(U, grid, gridRb, par_.mass));
}
