/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MScalar/FreeProp.cc

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
#include <Hadrons/Modules/MScalar/FreeProp.hpp>
#include <Hadrons/Modules/MScalar/Scalar.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MScalar;

/******************************************************************************
*                        TFreeProp implementation                             *
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
TFreeProp::TFreeProp(const std::string name)
: Module<FreePropPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
std::vector<std::string> TFreeProp::getInput(void)
{
    std::vector<std::string> in = {par().source};
    
    return in;
}

std::vector<std::string> TFreeProp::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
void TFreeProp::setup(void)
{
    freeMomPropName_ = FREEMOMPROP(par().mass);
    
    freePropDone_ = env().hasCreatedObject(freeMomPropName_);
    envCacheLat(ScalarField, freeMomPropName_);
    envCreateLat(ScalarField, getName());
}

// execution ///////////////////////////////////////////////////////////////////
void TFreeProp::execute(void)
{
    auto &freeMomProp = envGet(ScalarField, freeMomPropName_);
    auto &prop        = envGet(ScalarField, getName());
    auto &source      = envGet(ScalarField, par().source);

    if (!freePropDone_)
    {
        LOG(Message) << "Caching momentum space free scalar propagator"
                     << " (mass= " << par().mass << ")..." << std::endl;
        SIMPL::MomentumSpacePropagator(freeMomProp, par().mass);
    }
    LOG(Message) << "Computing free scalar propagator..." << std::endl;
    SIMPL::FreePropagator(source, prop, freeMomProp);
    
    if (!par().output.empty())
    {
        std::vector<TComplex> buf;
        std::vector<Complex>  result;
        
        sliceSum(prop, buf, Tp);
        result.resize(buf.size());
        for (unsigned int t = 0; t < buf.size(); ++t)
        {
            result[t] = TensorRemove(buf[t]);
        }
        saveResult(par().output, "freeprop", result);
    }
}
