/*******************************************************************************
Grid physics library, www.github.com/paboyle/Grid 

Source file: programs/Hadrons/SrcZ2.hpp

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

#ifndef Hadrons_SrcZ2_hpp_
#define Hadrons_SrcZ2_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/*
 
 Z_2 stochastic source
 -----------------------------
 * src_x = eta_x * theta(x_3 - ta) * theta(tb - x_3)
 
 * options:
 - tA: begin timeslice (integer)
 - tB: end timesilce (integer)
 
 */
 
/******************************************************************************
 *                               SrcZ2                                        *
 ******************************************************************************/
class SrcZ2Par: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(SrcZ2Par,
                                    unsigned int, tA,
                                    unsigned int, tB);
};

class SrcZ2: public Module<SrcZ2Par>
{
public:
    // constructor
    SrcZ2(const std::string name);
    // destructor
    virtual ~SrcZ2(void) = default;
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER(SrcZ2);

END_HADRONS_NAMESPACE

#endif // Hadrons_SrcZ2_hpp_
