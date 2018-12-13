/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MIO/LoadNersc.hpp

Copyright (C) 2015-2018

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
#ifndef Hadrons_MIO_LoadNersc_hpp_
#define Hadrons_MIO_LoadNersc_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                       Load a NERSC configuration                           *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MIO)

class LoadNerscPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(LoadNerscPar,
                                    std::string, file);
};

template <typename GImpl>
class TLoadNersc: public Module<LoadNerscPar>
{
public:
    GAUGE_TYPE_ALIASES(GImpl,);
public:
    // constructor
    TLoadNersc(const std::string name);
    // destructor
    virtual ~TLoadNersc(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(LoadNersc,  TLoadNersc<GIMPL>,  MIO);

/******************************************************************************
*                       TLoadNersc implementation                             *
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename GImpl>
TLoadNersc<GImpl>::TLoadNersc(const std::string name)
: Module<LoadNerscPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename GImpl>
std::vector<std::string> TLoadNersc<GImpl>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename GImpl>
std::vector<std::string> TLoadNersc<GImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename GImpl>
void TLoadNersc<GImpl>::setup(void)
{
    envCreateLat(GaugeField, getName());
}

// execution ///////////////////////////////////////////////////////////////////
template <typename GImpl>
void TLoadNersc<GImpl>::execute(void)
{
    FieldMetaData header;
    std::string   fileName = par().file + "."
                             + std::to_string(vm().getTrajectory());
    LOG(Message) << "Loading NERSC configuration from file '" << fileName
                 << "'" << std::endl;

    auto &U = envGet(GaugeField, getName());
    NerscIO::readConfiguration(U, header, fileName);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MIO_LoadNersc_hpp_
