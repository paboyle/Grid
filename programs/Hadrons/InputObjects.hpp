/*
 * InputObjects.hpp, part of Grid
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

#ifndef Hadrons_InputObjects_hpp_
#define Hadrons_InputObjects_hpp_

#include <Hadrons/Global.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                       Serializable input classes                           *
 ******************************************************************************/
class Module: Serializable
{
public:
    // constructor
    Module(void) = default;
    // destructor
    virtual ~Module(void) = default;
    // serializable members
    GRID_DECL_CLASS_MEMBERS(Module,
                            std::string             , name,
                            std::vector<std::string>, in
                            );
};

class Parameters: Serializable
{
public:
    // constructor
    Parameters(void) = default;
    // destructor
    virtual ~Parameters(void) = default;
    // serializable members
    GRID_DECL_CLASS_MEMBERS(Parameters,
                            std::vector<unsigned int>, latticeSize,
                            std::vector<Module>      , modules
                            );
};

END_HADRONS_NAMESPACE

#endif // Hadrons_InputObjects_hpp_
