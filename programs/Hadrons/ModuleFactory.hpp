/*
 * ModuleFactory.hpp, part of Grid
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
    typedef std::function<std::unique_ptr<Module>(const std::string &)>
        FactoryFunc;
public:
    // registration
    void registerModule(const std::string &type, const FactoryFunc &f);
    // get module list
    std::vector<std::string> getModuleList(void) const;
    // factory
    std::unique_ptr<Module> create(const std::string &type,
                                   const std::string &name) const;
private:
    std::map<std::string, FactoryFunc> factory_;
};

END_HADRONS_NAMESPACE

#endif // Hadrons_ModuleFactory_hpp_
