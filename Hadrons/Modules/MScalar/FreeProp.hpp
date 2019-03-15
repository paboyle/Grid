/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MScalar/FreeProp.hpp

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
#ifndef Hadrons_MScalar_FreeProp_hpp_
#define Hadrons_MScalar_FreeProp_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                               FreeProp                                     *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MScalar)

class FreePropPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(FreePropPar,
                                    std::string, source,
                                    double,      mass,
                                    std::string, output);
};

class TFreeProp: public Module<FreePropPar>
{
public:
    BASIC_TYPE_ALIASES(SIMPL,);
public:
    // constructor
    TFreeProp(const std::string name);
    // destructor
    virtual ~TFreeProp(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    std::string freeMomPropName_;
    bool        freePropDone_;
};

MODULE_REGISTER(FreeProp, TFreeProp, MScalar);

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MScalar_FreeProp_hpp_
