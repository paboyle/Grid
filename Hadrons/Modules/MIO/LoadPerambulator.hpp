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

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/EigenPack.hpp>
#include <Hadrons/Distil.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         LoadPerambulator                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MIO)

class LoadPerambulatorPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(LoadPerambulatorPar,
                                        std::string, PerambFileName, //stem!!!
	                                int, nvec,
	                                MDistil::DistilParameters, Distil);
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
TLoadPerambulator<FImpl>::TLoadPerambulator(const std::string name)
: Module<LoadPerambulatorPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TLoadPerambulator<FImpl>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TLoadPerambulator<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TLoadPerambulator<FImpl>::setup(void)
{
    const int nvec{par().nvec};
    const MDistil::DistilParameters & Distil{par().Distil};
    const int LI{Distil.LI};
    const int nnoise{Distil.nnoise};
    const int Nt_inv{Distil.Nt_inv}; // TODO: PROBABLY BETTER: if (full_tdil) Nt_inv=1; else Nt_inv = TI;
    const int Ns{Distil.Ns};
    std::array<std::string,6> sIndexNames{"Nt", "nvec", "LI", "nnoise", "Nt_inv", "SI"};

    envCreate(MDistil::Perambulator<SpinVector COMMA 6 COMMA sizeof(Real)>, getName(), 1,
		              sIndexNames,Distil.Nt,nvec,Distil.LI,Distil.nnoise,Distil.Nt_inv,Distil.SI);

}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TLoadPerambulator<FImpl>::execute(void)
{
  auto &perambulator = envGet(MDistil::Perambulator<SpinVector COMMA 6 COMMA sizeof(Real)>, getName());
  const std::string sPerambName{par().PerambFileName + "." + std::to_string(vm().getTrajectory())};
  perambulator.read(sPerambName.c_str());
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MIO_LoadPerambulator_hpp_
