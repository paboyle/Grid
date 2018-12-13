/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MIO/LoadCosmHol.hpp

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
#ifndef Hadrons_MIO_LoadCosmHol_hpp_
#define Hadrons_MIO_LoadCosmHol_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                    Load scalar SU(N) configurations                        *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MIO)

class LoadCosmHolPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(LoadCosmHolPar,
                                    std::string, file);
};

class ScalarActionParameters: Serializable 
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(ScalarActionParameters,
                                    double, mass_squared,
                                    double, lambda,
                                    double, g);
};

template <typename SImpl>
class TLoadCosmHol: public Module<LoadCosmHolPar>
{
public:
    typedef typename SImpl::Field Field;
public:
    // constructor
    TLoadCosmHol(const std::string name);
    // destructor
    virtual ~TLoadCosmHol(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(LoadCosmHolSU2, TLoadCosmHol<ScalarNxNAdjImplR<2>>, MIO);
MODULE_REGISTER_TMP(LoadCosmHolSU3, TLoadCosmHol<ScalarNxNAdjImplR<3>>, MIO);
MODULE_REGISTER_TMP(LoadCosmHolSU4, TLoadCosmHol<ScalarNxNAdjImplR<4>>, MIO);
MODULE_REGISTER_TMP(LoadCosmHolSU5, TLoadCosmHol<ScalarNxNAdjImplR<5>>, MIO);
MODULE_REGISTER_TMP(LoadCosmHolSU6, TLoadCosmHol<ScalarNxNAdjImplR<6>>, MIO);

/******************************************************************************
 *                       TLoadCosmHol implementation                          *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename SImpl>
TLoadCosmHol<SImpl>::TLoadCosmHol(const std::string name)
: Module<LoadCosmHolPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename SImpl>
std::vector<std::string> TLoadCosmHol<SImpl>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename SImpl>
std::vector<std::string> TLoadCosmHol<SImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename SImpl>
void TLoadCosmHol<SImpl>::setup(void)
{
    envCreateLat(Field, getName());
}

// execution ///////////////////////////////////////////////////////////////////
template <typename SImpl>
void TLoadCosmHol<SImpl>::execute(void)
{
    ScalarActionParameters    md;
    std::string        filename = par().file + "."
                                  + std::to_string(vm().getTrajectory());
    ScidacReader       reader;
    const unsigned int N    = SImpl::Group::Dimension;
    auto               &phi = envGet(Field, getName());

    LOG(Message) << "Loading CosmHol configuration from file '" << filename
                 << "'" << std::endl;
    reader.open(filename);
    reader.readScidacFieldRecord(phi, md);
    reader.close();
    LOG(Message) << "tr(phi^2) = " 
                 << -TensorRemove(sum(trace(phi*phi))).real()/env().getVolume() 
                 << std::endl;
    LOG(Message) << "Configuration parameters:" << std::endl;
    LOG(Message) << "     N = " << N << std::endl;
    LOG(Message) << "   m^2 = " << md.mass_squared << std::endl;
    LOG(Message) << "lambda = " << md.lambda << std::endl;
    LOG(Message) << "     g = " << md.g << std::endl;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MIO_LoadCosmHol_hpp_
