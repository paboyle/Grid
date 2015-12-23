/*
 * MQuark.hpp, part of Grid
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

#ifndef Hadrons_MQuark_hpp_
#define Hadrons_MQuark_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                               MQuark                                       *
 ******************************************************************************/
class MQuark: public Module
{
public:
    class Par: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Par, unsigned int, Ls);
    };
public:
    // constructor
    MQuark(const std::string &name);
    // destructor
    virtual ~MQuark(void) = default;
    // parse parameters
    virtual void parseParameters(XmlReader &reader, const std::string &name);
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // memory footprint
    virtual double nCreatedProp(void);
    // execution
    virtual void operator()(Environment &env);
private:
    Par par_;
};

MODULE_REGISTER(MQuark);

END_HADRONS_NAMESPACE

#endif // Hadrons_MQuark_hpp_
