/*************************************************************************************
 
 Grid physics library, www.github.com/paboyle/Grid
 
 Source file: Hadrons/Modules/MDistil/LoadDistilNoise.hpp
 
 Copyright (C) 2019
 
 Author: Felix Erben <ferben@ed.ac.uk>
 Author: Michael Marshall <Michael.Marshall@ed.ac.uk>
 
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

#ifndef Hadrons_MIO_LoadDistilNoise_hpp_
#define Hadrons_MIO_LoadDistilNoise_hpp_

#include <Hadrons/Modules/MDistil/Distil.hpp>

BEGIN_HADRONS_NAMESPACE
BEGIN_MODULE_NAMESPACE(MIO)

/******************************************************************************
 *                         LoadDistilNoise                                 *
 ******************************************************************************/

class LoadDistilNoisePar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(LoadDistilNoisePar,
                                        std::string, NoiseFileName,
                                        std::string, DistilParams);
};

template <typename FImpl>
class TLoadDistilNoise: public Module<LoadDistilNoisePar>
{
public:
    // constructor
    TLoadDistilNoise(const std::string name);
    // destructor
    virtual ~TLoadDistilNoise(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(LoadDistilNoise, TLoadDistilNoise<FIMPL>, MIO);

/******************************************************************************
 *                 TLoadDistilNoise implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TLoadDistilNoise<FImpl>::TLoadDistilNoise(const std::string name) : Module<LoadDistilNoisePar>(name) {}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TLoadDistilNoise<FImpl>::getInput(void)
{
    return {par().DistilParams};
}

template <typename FImpl>
std::vector<std::string> TLoadDistilNoise<FImpl>::getOutput(void)
{
    return {getName()};
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TLoadDistilNoise<FImpl>::setup(void)
{
    const MDistil::DistilParameters &dp{envGet(MDistil::DistilParameters,  par().DistilParams)};
    const int Nt{env().getDim(Tdir)};
    envCreate(MDistil::NoiseTensor, getName(), 1, dp.nnoise, Nt, dp.nvec, Ns);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TLoadDistilNoise<FImpl>::execute(void)
{
  auto &noises = envGet(MDistil::NoiseTensor, getName());
  std::string sNoiseName{ par().NoiseFileName };
  sNoiseName.append( 1, '.' );
  sNoiseName.append( std::to_string( vm().getTrajectory() ) );
  noises.read(sNoiseName.c_str());
}

END_MODULE_NAMESPACE
END_HADRONS_NAMESPACE
#endif
