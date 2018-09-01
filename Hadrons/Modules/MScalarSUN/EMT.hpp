/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MScalarSUN/EMT.hpp

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
#ifndef Hadrons_MScalarSUN_EMT_hpp_
#define Hadrons_MScalarSUN_EMT_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Modules/MScalarSUN/Utils.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         Energy-momentum tensor                             *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MScalarSUN)

class EMTPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(EMTPar,
                                    std::string, kinetic,
                                    std::string, phiPow,
                                    std::string, improvement,
                                    double     , m2,
                                    double     , lambda,
                                    double     , g,
                                    double     , xi,
                                    std::string, output);
};

class EMTResult: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(EMTResult,
                                    std::vector<std::vector<Complex>>, value,
                                    double,                            m2,
                                    double,                            lambda,
                                    double,                            g,
                                    double,                            xi);
};

template <typename SImpl>
class TEMT: public Module<EMTPar>
{
public:
    typedef typename SImpl::Field        Field;
    typedef typename SImpl::ComplexField ComplexField;
public:
    // constructor
    TEMT(const std::string name);
    // destructor
    virtual ~TEMT(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(EMTSU2, TEMT<ScalarNxNAdjImplR<2>>, MScalarSUN);
MODULE_REGISTER_TMP(EMTSU3, TEMT<ScalarNxNAdjImplR<3>>, MScalarSUN);
MODULE_REGISTER_TMP(EMTSU4, TEMT<ScalarNxNAdjImplR<4>>, MScalarSUN);
MODULE_REGISTER_TMP(EMTSU5, TEMT<ScalarNxNAdjImplR<5>>, MScalarSUN);
MODULE_REGISTER_TMP(EMTSU6, TEMT<ScalarNxNAdjImplR<6>>, MScalarSUN);

/******************************************************************************
 *                           TEMT implementation                              *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename SImpl>
TEMT<SImpl>::TEMT(const std::string name)
: Module<EMTPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename SImpl>
std::vector<std::string> TEMT<SImpl>::getInput(void)
{
    std::vector<std::string> in;
    
    for (unsigned int mu = 0; mu < env().getNd(); ++mu)
    for (unsigned int nu = mu; nu < env().getNd(); ++nu)
    {
        in.push_back(varName(par().kinetic, mu, nu));
        if (!par().improvement.empty())
        {
            in.push_back(varName(par().improvement, mu, nu));
        }
    }
    in.push_back(varName(par().kinetic, "sum"));
    in.push_back(varName(par().phiPow, 2));
    in.push_back(varName(par().phiPow, 4));

    return in;
}

template <typename SImpl>
std::vector<std::string> TEMT<SImpl>::getOutput(void)
{
    std::vector<std::string> out;
    
    for (unsigned int mu = 0; mu < env().getNd(); ++mu)
    for (unsigned int nu = mu; nu < env().getNd(); ++nu)
    {
        out.push_back(varName(getName(), mu, nu));
    }

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename SImpl>
void TEMT<SImpl>::setup(void)
{
    for (unsigned int mu = 0; mu < env().getNd(); ++mu)
    for (unsigned int nu = mu; nu < env().getNd(); ++nu)
    {
        envCreateLat(ComplexField, varName(getName(), mu, nu));
    }
}

// execution ///////////////////////////////////////////////////////////////////
template <typename SImpl>
void TEMT<SImpl>::execute(void)
{
    LOG(Message) << "Computing energy-momentum tensor" << std::endl;
    LOG(Message) << "  kinetic terms: '" << par().kinetic << "'" << std::endl;
    LOG(Message) << "      tr(phi^n): '" << par().phiPow << "'" << std::endl;
    if (!par().improvement.empty())
    {
        LOG(Message) << "    improvement: '" << par().improvement << "'" << std::endl;
    }
    LOG(Message) << "            m^2= " << par().m2 << std::endl;
    LOG(Message) << "         lambda= " << par().lambda << std::endl;
    LOG(Message) << "              g= " << par().g << std::endl;
    if (!par().improvement.empty())
    {
        LOG(Message) << "             xi= " << par().xi << std::endl;
    }

    const unsigned int N = SImpl::Group::Dimension, nd = env().getNd();
    auto               &trphi2 = envGet(ComplexField, varName(par().phiPow, 2));
    auto               &trphi4 = envGet(ComplexField, varName(par().phiPow, 4));
    auto               &sumkin = envGet(ComplexField, varName(par().kinetic, "sum"));
    EMTResult          result;

    if (!par().output.empty())
    {
        result.m2     = par().m2;
        result.g      = par().g;
        result.lambda = par().lambda;
        result.xi     = par().xi;
        result.value.resize(nd, std::vector<Complex>(nd));
    }
    for (unsigned int mu = 0; mu < nd; ++mu)
    for (unsigned int nu = mu; nu < nd; ++nu)
    {
        auto &out   = envGet(ComplexField, varName(getName(), mu, nu));
        auto &trkin = envGet(ComplexField, varName(par().kinetic, mu, nu));
        
        out = 2.*trkin;
        if (!par().improvement.empty())
        {
            auto &imp = envGet(ComplexField, varName(par().improvement, mu, nu));

            out += par().xi*imp;
        }
        if (mu == nu)
        {
            out -= sumkin + par().m2*trphi2 + par().lambda*trphi4;
        }
        out *= N/par().g;
        if (!par().output.empty())
        {
            result.value[mu][nu] = TensorRemove(sum(out));
            result.value[mu][nu] = result.value[nu][mu];
        }
    }
    if (!par().output.empty())
    {
        saveResult(par().output, "emt", result);
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MScalarSUN_EMT_hpp_
