/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MGauge/Electrify.hpp

Copyright (C) 2015-2018

Author: Vera Guelpers <Vera.Guelpers@ed.ac.uk>

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

#ifndef Hadrons_MGauge_Electrify_hpp_
#define Hadrons_MGauge_Electrify_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                              Electrify gauge                               *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MGauge)

/****************************************************************************
*  Electrify a gauge field:
*
*  Ue_mu(x) = U_mu(x)*exp(ieqA_mu(x))
*
*  with
*
*  - gauge: U_mu(x): gauge field
*  - emField: A_mu(x): electromagnetic photon field
*  - e: value for the elementary charge
*  - q: charge in units of e
*
*****************************************************************************/


class ElectrifyPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(ElectrifyPar,
                                    std::string, gauge,
				    std::string, emField,
				    double, e,
				    double, charge);
};

template <typename GImpl>
class TElectrify: public Module<ElectrifyPar>
{
public:
    GAUGE_TYPE_ALIASES(GImpl,);
public:
    typedef PhotonR::GaugeField     EmField;
public:
    // constructor
    TElectrify(const std::string name);
    // destructor
    virtual ~TElectrify(void) {};
    // dependencies/products
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(Electrify, TElectrify<GIMPL>, MGauge);

/******************************************************************************
*                            TElectrify implementation                             *
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename GImpl>
TElectrify<GImpl>::TElectrify(const std::string name)
: Module<ElectrifyPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename GImpl>
std::vector<std::string> TElectrify<GImpl>::getInput(void)
{
    std::vector<std::string> in = {par().gauge, par().emField};

    return in;
}

template <typename GImpl>
std::vector<std::string> TElectrify<GImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename GImpl>
void TElectrify<GImpl>::setup(void)
{
    envCreateLat(GaugeField, getName());
    envTmpLat(LatticeComplex, "eiAmu");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename GImpl>
void TElectrify<GImpl>::execute(void)
{
    LOG(Message) << "Electrify the gauge field " << par().gauge << " using the photon field " 
                  << par().emField << " with charge e*q= " << par().e << "*" << par().charge << std::endl;
    
    auto &Ue = envGet(GaugeField, getName());
    auto &U = envGet(GaugeField, par().gauge);
    auto &A = envGet(EmField,  par().emField);
    envGetTmp(LatticeComplex, eiAmu);

    Complex i(0.0,1.0);

    for(unsigned int mu = 0; mu < env().getNd(); mu++)
    {
	eiAmu = exp(i * (Real)(par().e * par().charge) * PeekIndex<LorentzIndex>(A, mu));
	PokeIndex<LorentzIndex>(Ue, PeekIndex<LorentzIndex>(U, mu) * eiAmu, mu);
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MGauge_Electrify_hpp_
