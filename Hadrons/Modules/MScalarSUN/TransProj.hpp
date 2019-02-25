/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MScalarSUN/TransProj.hpp

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
#ifndef Hadrons_MScalarSUN_TransProj_hpp_
#define Hadrons_MScalarSUN_TransProj_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Modules/MScalarSUN/Utils.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         Transverse projection                              *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MScalarSUN)

class TransProjPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(TransProjPar,
                                    std::string,  op,
                                    DiffType,     type,
                                    std::string,  output);
};

class TransProjResult: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(TransProjResult,
                                    std::vector<std::vector<Complex>>, value,
                                    DiffType,                          type);
};

template <typename SImpl>
class TTransProj: public Module<TransProjPar>
{
public:
    typedef typename SImpl::Field        Field;
    typedef typename SImpl::ComplexField ComplexField;
public:
    // constructor
    TTransProj(const std::string name);
    // destructor
    virtual ~TTransProj(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution

    virtual void execute(void);
};

MODULE_REGISTER_TMP(TransProjSU2, TTransProj<ScalarNxNAdjImplR<2>>, MScalarSUN);
MODULE_REGISTER_TMP(TransProjSU3, TTransProj<ScalarNxNAdjImplR<3>>, MScalarSUN);
MODULE_REGISTER_TMP(TransProjSU4, TTransProj<ScalarNxNAdjImplR<4>>, MScalarSUN);
MODULE_REGISTER_TMP(TransProjSU5, TTransProj<ScalarNxNAdjImplR<5>>, MScalarSUN);
MODULE_REGISTER_TMP(TransProjSU6, TTransProj<ScalarNxNAdjImplR<6>>, MScalarSUN);

/******************************************************************************
 *                        TTransProj implementation                           *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename SImpl>
TTransProj<SImpl>::TTransProj(const std::string name)
: Module<TransProjPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename SImpl>
std::vector<std::string> TTransProj<SImpl>::getInput(void)
{
    std::vector<std::string> in = {par().op};
    
    return in;
}

template <typename SImpl>
std::vector<std::string> TTransProj<SImpl>::getOutput(void)
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
void TTransProj<SImpl>::setup(void)
{
    for (unsigned int mu = 0; mu < env().getNd(); ++mu)
    for (unsigned int nu = mu; nu < env().getNd(); ++nu)
    {
        envCreateLat(ComplexField, varName(getName(), mu, nu));
    }
    envTmpLat(ComplexField, "buf1");
    envTmpLat(ComplexField, "buf2");
    envTmpLat(ComplexField, "lap");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename SImpl>
void TTransProj<SImpl>::execute(void)
{
    LOG(Message) << "Computing (delta_mu,nu d^2 - d_mu*d_nu)*op using " 
                 << par().type << " derivatives and op= '" << par().op 
                 << "'" << std::endl; 

    const unsigned int nd = env().getNd();
    TransProjResult    result;
    auto               &op = envGet(ComplexField, par().op);

    envGetTmp(ComplexField, buf1);
    envGetTmp(ComplexField, buf2);
    envGetTmp(ComplexField, lap);
    lap = zero;
    if (!par().output.empty())
    {
        result.type = par().type;
        result.value.resize(nd, std::vector<Complex>(nd));
    }
    for (unsigned int mu = 0; mu < nd; ++mu)
    {
        dmu(buf1, op, mu, par().type);
        dmu(buf2, buf1, mu, par().type);
        lap += buf2;
    }
    for (unsigned int mu = 0; mu < nd; ++mu)
    for (unsigned int nu = mu; nu < nd; ++nu)
    {
        auto &out = envGet(ComplexField, varName(getName(), mu, nu));
        dmu(buf1, op, mu, par().type);
        dmu(buf2, buf1, nu, par().type);
        out = -buf2;
        if (mu == nu)
        {
            out += lap;
        }
        if (!par().output.empty())
        {
            result.value[mu][nu] = TensorRemove(sum(out));
            result.value[mu][nu] = result.value[nu][mu];
        }
    }
    if (!par().output.empty())
    {
        saveResult(par().output, "transproj", result);
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MScalarSUN_TransProj_hpp_
