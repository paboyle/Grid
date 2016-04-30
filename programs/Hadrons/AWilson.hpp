/*******************************************************************************
Grid physics library, www.github.com/paboyle/Grid 

Source file: programs/Hadrons/AWilson.hpp

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

#ifndef Hadrons_AWilson_hpp_
#define Hadrons_AWilson_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/FermionAction.hpp>
#include <Hadrons/FermionActionFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                           Wilson fermions                                  *
 ******************************************************************************/
class AWilson: public FermionAction
{
public:
    class Par: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Par, double, mass);
    };
public:
    // constructor
    AWilson(const std::string name);
    // destructor
    virtual ~AWilson(void) = default;
    // parse parameters
    virtual void parseParameters(XmlReader &reader, const std::string name);
    // create operator
    virtual void create(Environment &env);
private:
    Par par_;
};

ACTION_REGISTER(AWilson);

END_HADRONS_NAMESPACE

#endif // Hadrons_AWilson_hpp_
