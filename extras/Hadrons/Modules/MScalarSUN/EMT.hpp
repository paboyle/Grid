/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: extras/Hadrons/Modules/MScalarSUN/EMT.hpp

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

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>
#include <Grid/Hadrons/Modules/MScalarSUN/Utils.hpp>

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
    virtual ~TEMT(void) = default;
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_NS(EMTSU2, TEMT<ScalarNxNAdjImplR<2>>, MScalarSUN);
MODULE_REGISTER_NS(EMTSU3, TEMT<ScalarNxNAdjImplR<3>>, MScalarSUN);
MODULE_REGISTER_NS(EMTSU4, TEMT<ScalarNxNAdjImplR<4>>, MScalarSUN);
MODULE_REGISTER_NS(EMTSU5, TEMT<ScalarNxNAdjImplR<5>>, MScalarSUN);
MODULE_REGISTER_NS(EMTSU6, TEMT<ScalarNxNAdjImplR<6>>, MScalarSUN);

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
        in.push_back(varName(par().improvement, mu, nu));
    }
    in.push_back(varName(par().phiPow, 2));
    in.push_back(varName(par().phiPow, 4));

    return in;
}

template <typename SImpl>
std::vector<std::string> TEMT<SImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
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
    envCreateLat(ComplexField, getName());
    envTmpLat(ComplexField, "sumkin");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename SImpl>
void TEMT<SImpl>::execute(void)
{
    LOG(Message) << "Computing energy-momentum tensor" << std::endl;
    LOG(Message) << "  kinetic terms: '" << par().kinetic << "'" << std::endl;
    LOG(Message) << "      tr(phi^n): '" << par().phiPow << "'" << std::endl;
    LOG(Message) << "    improvement: '" << par().improvement << "'" << std::endl;
    LOG(Message) << "            m^2= " << par().m2 << std::endl;
    LOG(Message) << "         lambda= " << par().lambda << std::endl;
    LOG(Message) << "              g= " << par().g << std::endl;
    LOG(Message) << "             xi= " << par().xi << std::endl;

    const unsigned int N       = SImpl::Group::Dimension;
    auto               &trphi2 = envGet(ComplexField, varName(par().phiPow, 2));
    auto               &trphi4 = envGet(ComplexField, varName(par().phiPow, 4));

    envGetTmp(ComplexField, sumkin);
    sumkin = zero;
    for (unsigned int mu = 0; mu < env().getNd(); ++mu)
    {
        auto &trkin = envGet(ComplexField, varName(par().kinetic, mu, mu));

        sumkin += trkin;
    }
    for (unsigned int mu = 0; mu < env().getNd(); ++mu)
    for (unsigned int nu = mu; nu < env().getNd(); ++nu)
    {
        auto &out   = envGet(ComplexField, varName(getName(), mu, nu));
        auto &trkin = envGet(ComplexField, varName(par().kinetic, mu, nu));
        auto &imp   = envGet(ComplexField, varName(par().improvement, mu, nu));

        out = 2.*trkin + par().xi*imp;
        if (mu == nu)
        {
            out -= sumkin + par().m2*trphi2 + par().lambda*trphi4;
        }
        out *= N/par().g;
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MScalarSUN_EMT_hpp_
