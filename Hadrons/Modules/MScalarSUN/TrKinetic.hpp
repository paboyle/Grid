/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MScalarSUN/TrKinetic.hpp

Copyright (C) 2015-2019

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
#ifndef Hadrons_MScalarSUN_TrKinetic_hpp_
#define Hadrons_MScalarSUN_TrKinetic_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Modules/MScalarSUN/Utils.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         Trace of kinetic term                              *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MScalarSUN)

class TrKineticPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(TrKineticPar,
                                    std::string,  field,
                                    DiffType,     type,
                                    std::string,  output);
};

class TrKineticResult: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(TrKineticResult,
                                    std::vector<std::vector<Complex>>, value,
                                    DiffType,                          type);
};

template <typename SImpl>
class TTrKinetic: public Module<TrKineticPar>
{
public:
    typedef typename SImpl::Field        Field;
    typedef typename SImpl::ComplexField ComplexField;
public:
    // constructor
    TTrKinetic(const std::string name);
    // destructor
    virtual ~TTrKinetic(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(TrKineticSU2, TTrKinetic<ScalarNxNAdjImplR<2>>, MScalarSUN);
MODULE_REGISTER_TMP(TrKineticSU3, TTrKinetic<ScalarNxNAdjImplR<3>>, MScalarSUN);
MODULE_REGISTER_TMP(TrKineticSU4, TTrKinetic<ScalarNxNAdjImplR<4>>, MScalarSUN);
MODULE_REGISTER_TMP(TrKineticSU5, TTrKinetic<ScalarNxNAdjImplR<5>>, MScalarSUN);
MODULE_REGISTER_TMP(TrKineticSU6, TTrKinetic<ScalarNxNAdjImplR<6>>, MScalarSUN);

/******************************************************************************
 *                      TTrKinetic implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename SImpl>
TTrKinetic<SImpl>::TTrKinetic(const std::string name)
: Module<TrKineticPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename SImpl>
std::vector<std::string> TTrKinetic<SImpl>::getInput(void)
{
    std::vector<std::string> in = {par().field};
    
    return in;
}

template <typename SImpl>
std::vector<std::string> TTrKinetic<SImpl>::getOutput(void)
{
    std::vector<std::string> out ;

    for (unsigned int mu = 0; mu < env().getNd(); ++mu)
    for (unsigned int nu = mu; nu < env().getNd(); ++nu)
    {
        out.push_back(varName(getName(), mu, nu));
    }
    out.push_back(varName(getName(), "sum"));
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename SImpl>
void TTrKinetic<SImpl>::setup(void)
{
    for (unsigned int mu = 0; mu < env().getNd(); ++mu)
    for (unsigned int nu = mu; nu < env().getNd(); ++nu)
    {
        envCreateLat(ComplexField, varName(getName(), mu, nu));
    }
    envCreateLat(ComplexField, varName(getName(), "sum"));
    envTmp(std::vector<Field>, "der", 1, env().getNd(), envGetGrid(Field));
}

// execution ///////////////////////////////////////////////////////////////////
template <typename SImpl>
void TTrKinetic<SImpl>::execute(void)
{
    LOG(Message) << "Computing tr(d_mu phi*d_nu phi) using " << par().type
                 << " derivative" << std::endl; 

    const unsigned int nd = env().getNd();
    TrKineticResult    result;
    auto               &phi    = envGet(Field, par().field);
    auto               &sumkin = envGet(ComplexField, varName(getName(), "sum"));

    envGetTmp(std::vector<Field>, der);
    sumkin = zero;
    if (!par().output.empty())
    {
        result.type = par().type;
        result.value.resize(nd, std::vector<Complex>(nd));
    }
    for (unsigned int mu = 0; mu < nd; ++mu)
    {
        dmu(der[mu], phi, mu, par().type);
    }
    for (unsigned int mu = 0; mu < nd; ++mu)
    for (unsigned int nu = mu; nu < nd; ++nu)
    {
        auto &out = envGet(ComplexField, varName(getName(), mu, nu));

        out = -trace(der[mu]*der[nu]);
        if (mu == nu)
        {
            sumkin += out;
        }
        if (!par().output.empty())
        {
            result.value[mu][nu] = TensorRemove(sum(out));
            result.value[mu][nu] = result.value[nu][mu];
        }
    }
    saveResult(par().output, "trkinetic", result);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MScalarSUN_TrKinetic_hpp_
