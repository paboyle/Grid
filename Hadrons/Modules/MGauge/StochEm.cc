/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MGauge/StochEm.cc

Copyright (C) 2015-2019

Author: Antonin Portelli <antonin.portelli@me.com>
Author: James Harrison <jch1g10@soton.ac.uk>
Author: Vera Guelpers <vmg1n14@soton.ac.uk>

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
#include <Hadrons/Modules/MGauge/StochEm.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MGauge;

/******************************************************************************
*                  TStochEm implementation                             *
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
TStochEm::TStochEm(const std::string name)
: Module<StochEmPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
std::vector<std::string> TStochEm::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

std::vector<std::string> TStochEm::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
void TStochEm::setup(void)
{
    weightDone_ = env().hasCreatedObject("_" + getName() + "_weight");
    envCacheLat(EmComp, "_" + getName() + "_weight");
    envCreateLat(EmField, getName());
}

// execution ///////////////////////////////////////////////////////////////////
void TStochEm::execute(void)
{
    LOG(Message) << "Generating stochastic EM potential..." << std::endl;

    std::vector<Real> improvements = strToVec<Real>(par().improvement);
    PhotonR photon(envGetGrid(EmField), par().gauge, par().zmScheme, improvements);
    auto    &a = envGet(EmField, getName());
    auto    &w = envGet(EmComp, "_" + getName() + "_weight");
    
    if (!weightDone_)
    {
        LOG(Message) << "Caching stochastic EM potential weight (gauge: "
                     << par().gauge << ", zero-mode scheme: "
                     << par().zmScheme << ")..." << std::endl;
        photon.StochasticWeight(w);
    }
    photon.StochasticField(a, rng4d(), w);
}
