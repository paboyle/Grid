/*************************************************************************************
 
 Grid physics library, www.github.com/paboyle/Grid
 
 Source file: Hadrons/Modules/MDistil/LoadPerambulator.hpp
 
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

#ifndef Hadrons_MIO_LoadPerambulator_hpp_
#define Hadrons_MIO_LoadPerambulator_hpp_

#include <Hadrons/Modules/MDistil/Distil.hpp>

BEGIN_HADRONS_NAMESPACE
BEGIN_MODULE_NAMESPACE(MIO)

/******************************************************************************
 *                         LoadPerambulator                                 *
 ******************************************************************************/

class LoadPerambulatorPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(LoadPerambulatorPar,
                                        std::string, PerambFileName,
                                        std::string, DistilParams);
};

template <typename FImpl>
class TLoadPerambulator: public Module<LoadPerambulatorPar>
{
public:
    // constructor
    TLoadPerambulator(const std::string name);
    // destructor
    virtual ~TLoadPerambulator(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(LoadPerambulator, TLoadPerambulator<FIMPL>, MIO);

/******************************************************************************
 *                 TLoadPerambulator implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TLoadPerambulator<FImpl>::TLoadPerambulator(const std::string name) : Module<LoadPerambulatorPar>(name) {}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TLoadPerambulator<FImpl>::getInput(void)
{
    return {par().DistilParams};
}

template <typename FImpl>
std::vector<std::string> TLoadPerambulator<FImpl>::getOutput(void)
{
    return {getName()};
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TLoadPerambulator<FImpl>::setup(void)
{
    const MDistil::DistilParameters &dp{envGet(MDistil::DistilParameters,  par().DistilParams)};
    const int Nt{env().getDim(Tdir)}; 
    const bool full_tdil{ dp.TI == Nt };
    const int Nt_inv{ full_tdil ? 1 : dp.TI };
    envCreate(MDistil::PerambTensor, getName(), 1, Nt,dp.nvec,dp.LI,dp.nnoise,Nt_inv,dp.SI);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TLoadPerambulator<FImpl>::execute(void)
{
  auto &perambulator = envGet(MDistil::PerambTensor, getName());
  std::string sPerambName{ par().PerambFileName };
  sPerambName.append( 1, '.' );
  sPerambName.append( std::to_string( vm().getTrajectory() ) );
  perambulator.read(sPerambName.c_str());
}

END_MODULE_NAMESPACE
END_HADRONS_NAMESPACE
#endif
