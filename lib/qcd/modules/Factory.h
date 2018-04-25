/*************************************************************************************
Grid physics library, www.github.com/paboyle/Grid 
Source file: Factory.h

Copyright (C) 2015
Copyright (C) 2016

Author: Antonin Portelli <antonin.portelli@me.com>
Author: Guido Cossu <guido.cossu@ed.ac.uk>

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

#ifndef Factory_hpp_
#define Factory_hpp_


namespace Grid{




/******************************************************************************
 *                        abstract factory class                              *
 ******************************************************************************/
template <typename T, typename CreatorInput>
class Factory
{
public:
    typedef std::function< std::unique_ptr<T> (const CreatorInput&) > Func;

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
              								const CreatorInput& input) const;
private:
    std::map<std::string, Func> builder_;
    virtual std::string obj_type() const = 0;
};

/******************************************************************************
 *                         template implementation                            *
 ******************************************************************************/
// registration ////////////////////////////////////////////////////////////////
template <typename T, typename CreatorInput>
void Factory<T, CreatorInput>::registerBuilder(const std::string type, const Func &f)
{
    builder_[type] = f;
}

// get module list /////////////////////////////////////////////////////////////
template <typename T, typename CreatorInput>
std::vector<std::string> Factory<T, CreatorInput>::getBuilderList(void) const
{
    std::vector<std::string> list;
    
    for (auto &b: builder_)
    {
        list.push_back(b.first);
    }
    
    return list;
}

// factory /////////////////////////////////////////////////////////////////////
template <typename T, typename CreatorInput>
std::unique_ptr<T> Factory<T, CreatorInput>::create(const std::string type,
                                      							const CreatorInput& input) const
{
    Func func;
    
    std::cout << GridLogDebug << "Creating object of type "<< type << std::endl;
    try
    {
        func = builder_.at(type);
    }
    catch (std::out_of_range &)
    {
      //HADRONS_ERROR("object of type '" + type + "' unknown");
    	std::cout << GridLogError << "Error" << std::endl;
    	std::cout << GridLogError << obj_type() << " object of name [" << type << "] unknown" << std::endl;
    	exit(1);
    }
    
    return func(input);
}

}

#endif // Factory_hpp_
