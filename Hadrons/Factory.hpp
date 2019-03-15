/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Factory.hpp

Copyright (C) 2015-2019

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

#ifndef Hadrons_Factory_hpp_
#define Hadrons_Factory_hpp_

#include <Hadrons/Global.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                        abstract factory class                              *
 ******************************************************************************/
template <typename T>
class Factory
{
public:
    typedef std::function<std::unique_ptr<T>(const std::string)> Func;
public:
    // constructor
    Factory(void) = default;
    // destructor
    virtual ~Factory(void) = default;
    // registration
    void registerBuilder(const std::string type, const Func &f);
    // get builder list
    std::vector<std::string> getBuilderList(void) const;
    // factory
    std::unique_ptr<T> create(const std::string type,
                              const std::string name) const;
private:
    std::map<std::string, Func> builder_;
};

/******************************************************************************
 *                         template implementation                            *
 ******************************************************************************/
// registration ////////////////////////////////////////////////////////////////
template <typename T>
void Factory<T>::registerBuilder(const std::string type, const Func &f)
{
    builder_[type] = f;
}

// get module list /////////////////////////////////////////////////////////////
template <typename T>
std::vector<std::string> Factory<T>::getBuilderList(void) const
{
    std::vector<std::string> list;
    
    for (auto &b: builder_)
    {
        list.push_back(b.first);
    }
    
    return list;
}

// factory /////////////////////////////////////////////////////////////////////
template <typename T>
std::unique_ptr<T> Factory<T>::create(const std::string type,
                                      const std::string name) const
{
    Func func;
    
    try
    {
        func = builder_.at(type);
    }
    catch (std::out_of_range &)
    {
        HADRONS_ERROR(Argument, "object of type '" + type + "' unknown");
    }
    
    return func(name);
}

END_HADRONS_NAMESPACE

#endif // Hadrons_Factory_hpp_
