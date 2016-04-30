/*******************************************************************************
Grid physics library, www.github.com/paboyle/Grid 

Source file: programs/Hadrons/FermionAction.hpp

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

#ifndef Hadrons_FermionAction_hpp_
#define Hadrons_FermionAction_hpp_

#include <Hadrons/Global.hpp>

BEGIN_HADRONS_NAMESPACE

// action registration macro
#define ACTION_REGISTER(action)\
class action##ActionRegistrar\
{\
public:\
    action##ActionRegistrar(void)\
    {\
        FermionActionFactory &actionFac = FermionActionFactory::getInstance();\
        actionFac.registerBuilder(#action, [&](const std::string name)\
                                  {\
                                      return std::unique_ptr<action>(\
                                          new action(name));\
                                  });\
    }\
};\
static action##ActionRegistrar action##ActionRegistrarInstance;

/******************************************************************************
 *                             FermionAction                                  *
 ******************************************************************************/
class Environment;

class FermionAction
{
public:
    typedef FermionOperator<WilsonImplR> FMat;
    typedef std::unique_ptr<FMat>        FMatPt;
public:
    // constructor
    FermionAction(const std::string name);
    // destructor
    virtual ~FermionAction(void) = default;
    // access
            std::string  getName(void) const;
    virtual unsigned int getLs(void) const;
            void         setFMat(FMat *fMat);
            FMat *       getFMat(void);
    // parse parameters
    virtual void parseParameters(XmlReader &reader, const std::string name) = 0;
    // create operator
    virtual void create(Environment &env) = 0;
private:
    std::string name_;
    FMatPt      fMat_;
};

END_HADRONS_NAMESPACE

#endif // Hadrons_FermionAction_hpp_
