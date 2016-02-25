/*******************************************************************************
Grid physics library, www.github.com/paboyle/Grid 

Source file: programs/Hadrons/ModuleFactory.cc

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

#include <Hadrons/ModuleFactory.hpp>

using namespace Grid;
using namespace Hadrons;

/******************************************************************************
 *                      ModuleFactory implementation                          *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
ModuleFactory::ModuleFactory(void)
{
}

// registration ////////////////////////////////////////////////////////////////
void ModuleFactory::registerModule(const std::string &type,
                                   const FactoryFunc &f)
{
    factory_[type] = f;
}

// get module list /////////////////////////////////////////////////////////////
std::vector<std::string> ModuleFactory::getModuleList(void) const
{
    std::vector<std::string> list;
    
    for (auto &f: factory_)
    {
        list.push_back(f.first);
    }
    
    return list;
}

// factory /////////////////////////////////////////////////////////////////////
std::unique_ptr<Module> ModuleFactory::create(const std::string &type,
                                              const std::string &name) const
{
    FactoryFunc func;
    
    try
    {
        func = factory_.at(type);
    }
    catch (std::out_of_range)
    {
        HADRON_ERROR("module type '" + type + "' unknown");
    }
    
    return func(name);
}
