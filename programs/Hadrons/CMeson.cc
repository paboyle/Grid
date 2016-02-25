/*******************************************************************************
Grid physics library, www.github.com/paboyle/Grid 

Source file: programs/Hadrons/CMeson.cc

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

#include <Hadrons/CMeson.hpp>

using namespace Grid;
using namespace Hadrons;

/******************************************************************************
 *                          CMeson implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
CMeson::CMeson(const std::string &name)
: Module(name)
{}

// parse parameters ////////////////////////////////////////////////////////////
void CMeson::parseParameters(XmlReader &reader, const std::string &name)
{
    read(reader, name, par_);
}

// dependency relation /////////////////////////////////////////////////////////
std::vector<std::string> CMeson::getInput(void)
{
    std::vector<std::string> input = {par_.q1, par_.q2};
    
    return input;
}

std::vector<std::string> CMeson::getOutput(void)
{
    std::vector<std::string> output = {getName()};
    
    return output;
}

// memory footprint ////////////////////////////////////////////////////////////
void CMeson::allocate(Environment &env)
{}

// execution ///////////////////////////////////////////////////////////////////
void CMeson::execute(Environment &env)
{
    LOG(Message) << "computing meson contraction '" << getName() << "'"
                 << std::endl;
}
