/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MScalarSUN/Div.hpp

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

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Modules/MScalarSUN/Utils.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                       Divergence of a vector field                         *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MScalarSUN)

class DivPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(DivPar,
                                    std::vector<std::string>, op,
                                    DiffType,                 type,
                                    std::string,              output);
};

class DivResult: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(DivResult,
                                    DiffType, type,
                                    Complex,  value);
};

template <typename SImpl>
class TDiv: public Module<DivPar>
{
public:
    typedef typename SImpl::Field        Field;
    typedef typename SImpl::ComplexField ComplexField;
public:
    // constructor
    TDiv(const std::string name);
    // destructor
    virtual ~TDiv(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(DivSU2, TDiv<ScalarNxNAdjImplR<2>>, MScalarSUN);
MODULE_REGISTER_TMP(DivSU3, TDiv<ScalarNxNAdjImplR<3>>, MScalarSUN);
MODULE_REGISTER_TMP(DivSU4, TDiv<ScalarNxNAdjImplR<4>>, MScalarSUN);
MODULE_REGISTER_TMP(DivSU5, TDiv<ScalarNxNAdjImplR<5>>, MScalarSUN);
MODULE_REGISTER_TMP(DivSU6, TDiv<ScalarNxNAdjImplR<6>>, MScalarSUN);

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
        HADRONS_ERROR(Size, "the number of components differs from number of dimensions");
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
        std::cout << "'" << par().op[mu] << ((mu == nd - 1) ? "']" : "', ");
    }
    std::cout << std::endl;

    auto &div = envGet(ComplexField, getName());
    div = zero;
    for (unsigned int mu = 0; mu < nd; ++mu)
    {
        auto &op = envGet(ComplexField, par().op[mu]);
        dmuAcc(div, op, mu, par().type);
    }
    if (!par().output.empty())
    {
        DivResult r;

        r.type  = par().type;
        r.value = TensorRemove(sum(div));
        saveResult(par().output, "div", r);
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MScalarSUN_Div_hpp_
