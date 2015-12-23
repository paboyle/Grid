/*
 * ModuleFactory.cc, part of Grid
 *
 * Copyright (C) 2015 Antonin Portelli
 *
 * Grid is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Grid is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Grid.  If not, see <http://www.gnu.org/licenses/>.
 */

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
