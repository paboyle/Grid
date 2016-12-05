/*******************************************************************************
Grid physics library, www.github.com/paboyle/Grid 

Source file: programs/Hadrons/Z2.hpp

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

#ifndef Hadrons_Z2_hpp_
#define Hadrons_Z2_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>

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
 *                                 Z2                                         *
 ******************************************************************************/
namespace MSource
{
    class Z2Par: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Z2Par,
                                        unsigned int, tA,
                                        unsigned int, tB);
    };

    class Z2: public Module<Z2Par>
    {
    public:
        // constructor
        Z2(const std::string name);
        // destructor
        virtual ~Z2(void) = default;
        // dependency relation
        virtual std::vector<std::string> getInput(void);
        virtual std::vector<std::string> getOutput(void);
        // setup
        virtual void setup(void);
        // execution
        virtual void execute(void);
    };
}

MODULE_REGISTER_NS(Z2, MSource);

END_HADRONS_NAMESPACE

#endif // Hadrons_Z2_hpp_
