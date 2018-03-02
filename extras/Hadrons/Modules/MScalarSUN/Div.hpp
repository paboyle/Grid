/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: extras/Hadrons/Modules/MScalarSUN/Div.hpp

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
#ifndef Hadrons_MScalarSUN_Div_hpp_
#define Hadrons_MScalarSUN_Div_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                       Divergence of a vector field                         *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MScalarSUN)

class DivPar: Serializable
{
public:
    GRID_SERIALIZABLE_ENUM(DiffType, undef, forward, 1, backward, 2, central, 3);
    GRID_SERIALIZABLE_CLASS_MEMBERS(DivPar,
                                    std::vector<std::string>, op,
                                    DiffType,                 type,
                                    std::string,              output);
};

template <typename SImpl>
class TDiv: public Module<DivPar>
{
public:
    typedef typename SImpl::Field        Field;
    typedef typename SImpl::ComplexField ComplexField;
    class Result: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
                                        DivPar::DiffType, type,
                                        Complex,          value);
    };
public:
    // constructor
    TDiv(const std::string name);
    // destructor
    virtual ~TDiv(void) = default;
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_NS(DivSU2, TDiv<ScalarNxNAdjImplR<2>>, MScalarSUN);
MODULE_REGISTER_NS(DivSU3, TDiv<ScalarNxNAdjImplR<3>>, MScalarSUN);
MODULE_REGISTER_NS(DivSU4, TDiv<ScalarNxNAdjImplR<4>>, MScalarSUN);
MODULE_REGISTER_NS(DivSU5, TDiv<ScalarNxNAdjImplR<5>>, MScalarSUN);
MODULE_REGISTER_NS(DivSU6, TDiv<ScalarNxNAdjImplR<6>>, MScalarSUN);

/******************************************************************************
 *                           TDiv implementation                              *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename SImpl>
TDiv<SImpl>::TDiv(const std::string name)
: Module<DivPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename SImpl>
std::vector<std::string> TDiv<SImpl>::getInput(void)
{
    return par().op;
}

template <typename SImpl>
std::vector<std::string> TDiv<SImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename SImpl>
void TDiv<SImpl>::setup(void)
{
    if (par().op.size() != env().getNd())
    {
        HADRON_ERROR(Size, "the number of components differs from number of dimensions");
    }
    envCreateLat(ComplexField, getName());
}

// execution ///////////////////////////////////////////////////////////////////
template <typename SImpl>
void TDiv<SImpl>::execute(void)
{
    const auto nd = env().getNd();

    LOG(Message) << "Computing the " << par().type << " divergence of [";
    for (unsigned int mu = 0; mu < nd; ++mu)
    {
        std::cout << par().op[mu] << ((mu == nd - 1) ? "]" : ", ");
    }
    std::cout << std::endl;

    auto &div = envGet(ComplexField, getName());
    div = zero;
    for (unsigned int mu = 0; mu < nd; ++mu)
    {
        auto &op = envGet(ComplexField, par().op[mu]);
        switch(par().type)
        {
            case DivPar::DiffType::backward:
                div += op - Cshift(op, mu, -1);
                break;
            case DivPar::DiffType::forward:
                div += Cshift(op, mu, 1) - op;
                break;
            case DivPar::DiffType::central:
                div += 0.5*(Cshift(op, mu, 1) - Cshift(op, mu, -1));
                break;
        }
    }
    if (!par().output.empty())
    {
        Result       r;

        r.type  = par().type;
        r.value = TensorRemove(sum(div));
        saveResult(par().output, "div", r);
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MScalarSUN_Div_hpp_
