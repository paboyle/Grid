/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: extras/Hadrons/Modules/MScalarSUN/TrKinetic.hpp

Copyright (C) 2015-2018


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

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         TrKinetic                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MScalarSUN)

class TrKineticPar: Serializable
{
public:
    GRID_SERIALIZABLE_ENUM(DiffType, undef, forward, 1, backward, 2, central, 3);
    GRID_SERIALIZABLE_CLASS_MEMBERS(TrKineticPar,
                                    std::string,  field,
                                    DiffType,     type,
                                    std::string,  output);
};

template <typename SImpl>
class TTrKinetic: public Module<TrKineticPar>
{
public:
    typedef typename SImpl::Field        Field;
    typedef typename SImpl::ComplexField ComplexField;
    class Result: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
                                        std::string, op,
                                        Complex    , value);
    };
public:
    // constructor
    TTrKinetic(const std::string name);
    // destructor
    virtual ~TTrKinetic(void) = default;
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    std::string outName(const unsigned int mu, const unsigned int nu);
    std::string bufName(const unsigned int mu);
};

MODULE_REGISTER_NS(TrKineticSU2, TTrKinetic<ScalarNxNAdjImplR<2>>, MScalarSUN);
MODULE_REGISTER_NS(TrKineticSU3, TTrKinetic<ScalarNxNAdjImplR<3>>, MScalarSUN);
MODULE_REGISTER_NS(TrKineticSU4, TTrKinetic<ScalarNxNAdjImplR<4>>, MScalarSUN);
MODULE_REGISTER_NS(TrKineticSU5, TTrKinetic<ScalarNxNAdjImplR<5>>, MScalarSUN);
MODULE_REGISTER_NS(TrKineticSU6, TTrKinetic<ScalarNxNAdjImplR<6>>, MScalarSUN);

/******************************************************************************
 *                 TTrKinetic implementation                             *
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
        out.push_back(outName(mu, nu));
    }
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename SImpl>
void TTrKinetic<SImpl>::setup(void)
{
    for (unsigned int mu = 0; mu < env().getNd(); ++mu)
    for (unsigned int nu = mu; nu < env().getNd(); ++nu)
    {
        envCreateLat(ComplexField, outName(mu, nu));
    }
    envTmp(std::vector<Field>, "der", 1, env().getNd(), env().getGrid());
}

// execution ///////////////////////////////////////////////////////////////////
template <typename SImpl>
void TTrKinetic<SImpl>::execute(void)
{
    LOG(Message) << "Computing tr(d_mu phi*d_nu phi) using " << par().type
                 << " derivative" << std::endl; 

    std::vector<Result> result;
    auto                &phi = envGet(Field, par().field);

    envGetTmp(std::vector<Field>, der);
    for (unsigned int mu = 0; mu < env().getNd(); ++mu)
    {
        switch(par().type)
        {
            case TrKineticPar::DiffType::backward:
                der[mu] = phi - Cshift(phi, mu, -1);
                break;
            case TrKineticPar::DiffType::forward:
                der[mu] = Cshift(phi, mu, 1) - phi;
                break;
            case TrKineticPar::DiffType::central:
                der[mu] = 0.5*(Cshift(phi, mu, 1) - Cshift(phi, mu, -1));
                break;
        }
    }
    for (unsigned int mu = 0; mu < env().getNd(); ++mu)
    for (unsigned int nu = mu; nu < env().getNd(); ++nu)
    {
        auto &out = envGet(ComplexField, outName(mu, nu));

        out = -trace(der[mu]*der[nu]);
        if (!par().output.empty())
        {
            Result r;

            r.op    = "tr(d_" + std::to_string(mu) + "phi*d_" 
                      + std::to_string(nu) + "phi)";
            r.value = TensorRemove(sum(out));
            result.push_back(r);
        }
    }
    if (result.size() > 0)
    {
        saveResult(par().output, "trkinetic", result);
    }
}

// variable name generators ////////////////////////////////////////////////////
template <typename SImpl>
std::string TTrKinetic<SImpl>::outName(const unsigned int mu, 
                                       const unsigned int nu)
{
    return getName() + "_" + std::to_string(mu) + "_" + std::to_string(nu);
}

template <typename SImpl>
std::string TTrKinetic<SImpl>::bufName(const unsigned int mu)
{
    return "d_" + std::to_string(mu);
}


END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MScalarSUN_TrKinetic_hpp_
