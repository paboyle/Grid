/*******************************************************************************
Grid physics library, www.github.com/paboyle/Grid 

Source file: programs/Hadrons/MSource.hpp

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

/************
 * Sources  *
 ************
 
 Description of all source types.
 Convention: the discrete Heavyside function verifies theta(0) = 1.
 
 point: Point source
 -------------------
 * src(x) = delta_x,o
 
 * arguments: o
    - o: origin, space-separated integer sequence (e.g. "0 1 1 0")
 
 z2Band: Z_2 stochastic source
 -----------------------------
 * src(x) = eta_x * theta(x_0 - ta) * theta(tb - x_0)
 
 * arguments: ta tb
    - ta: begin timeslice (integer)
    - tb: end timesilce (integer)
 
 */


#ifndef Hadrons_MSource_hpp_
#define Hadrons_MSource_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

namespace Grid{
    GRID_SERIALIZABLE_ENUM(SourceType, undef,
                                       point,  1,
                                       z2Band, 2);
}

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                            Source module                                   *
 ******************************************************************************/
class MSource: public Module
{
public:
    class Par: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Par, SourceType, sourceType,
                                        std::vector<std::string>, arguments);
    };
public:
    // constructor
    MSource(const std::string name);
    // destructor
    virtual ~MSource(void) = default;
    // parse parameters
    virtual void parseParameters(XmlReader &reader, const std::string name);
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // allocation
    virtual void allocate(Environment &env);
    // execution
    virtual void execute(Environment &env);
private:
    Par               par_;
    LatticePropagator *src_{nullptr};
};

MODULE_REGISTER(MSource);

END_HADRONS_NAMESPACE

#endif // Hadrons_MSource_hpp_
