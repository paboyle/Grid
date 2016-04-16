/*******************************************************************************
Grid physics library, www.github.com/paboyle/Grid 

Source file: programs/Hadrons/ModuleFactory.hpp

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

#ifndef Hadrons_ModuleFactory_hpp_
#define Hadrons_ModuleFactory_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                            ModuleFactory                                   *
 ******************************************************************************/
class ModuleFactory
{
    SINGLETON(ModuleFactory)
public:
    typedef std::function<std::unique_ptr<Module>(const std::string )>
        FactoryFunc;
public:
    // registration
    void registerModule(const std::string type, const FactoryFunc &f);
    // get module list
    std::vector<std::string> getModuleList(void) const;
    // factory
    std::unique_ptr<Module> create(const std::string type,
                                   const std::string name) const;
private:
    std::map<std::string, FactoryFunc> factory_;
};

END_HADRONS_NAMESPACE

#endif // Hadrons_ModuleFactory_hpp_
