/*
 * CMeson.hpp, part of Grid
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

#ifndef Hadrons_CMeson_hpp_
#define Hadrons_CMeson_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                               CMeson                                       *
 ******************************************************************************/
class CMeson: public Module
{
public:
    class Par: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Par,
                                        std::string, q1,
                                        std::string, q2,
                                        std::string, output);
    };
public:
    // constructor
    CMeson(const std::string &name);
    // destructor
    virtual ~CMeson(void) = default;
    // parse parameters
    virtual void parseParameters(XmlReader &reader, const std::string &name);
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // allocation
    virtual void allocate(Environment &env);
    // execution
    virtual void execute(Environment &env);
private:
    Par par_;
};

MODULE_REGISTER(CMeson);

END_HADRONS_NAMESPACE

#endif // Hadrons_CMeson_hpp_
