/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MNoise/TimeDilutedSpinColorDiagonal.hpp

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
#ifndef Hadrons_MNoise_TimeDilutedSpinColorDiagonal_hpp_
#define Hadrons_MNoise_TimeDilutedSpinColorDiagonal_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/DilutedNoise.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *             Generate time diluted spin-color diagonal noise                *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MNoise)

template <typename FImpl>
class TTimeDilutedSpinColorDiagonal: public Module<NoPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TTimeDilutedSpinColorDiagonal(const std::string name);
    // destructor
    virtual ~TTimeDilutedSpinColorDiagonal(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(TimeDilutedSpinColorDiagonal, TTimeDilutedSpinColorDiagonal<FIMPL>, MNoise);
MODULE_REGISTER_TMP(ZTimeDilutedSpinColorDiagonal, TTimeDilutedSpinColorDiagonal<ZFIMPL>, MNoise);

/******************************************************************************
 *              TTimeDilutedSpinColorDiagonal implementation                  *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TTimeDilutedSpinColorDiagonal<FImpl>::TTimeDilutedSpinColorDiagonal(const std::string name)
: Module<NoPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TTimeDilutedSpinColorDiagonal<FImpl>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TTimeDilutedSpinColorDiagonal<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TTimeDilutedSpinColorDiagonal<FImpl>::setup(void)
{
    envCreateDerived(DilutedNoise<FImpl>, 
                     TimeDilutedSpinColorDiagonalNoise<FImpl>,
                     getName(), 1, envGetGrid(FermionField));
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TTimeDilutedSpinColorDiagonal<FImpl>::execute(void)
{
    auto &noise = envGet(DilutedNoise<FImpl>, getName());

    LOG(Message) << "Generating time-diluted, spin-color diagonal noise" << std::endl;
    noise.generateNoise(rng4d());
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MNoise_TimeDilutedSpinColorDiagonal_hpp_
