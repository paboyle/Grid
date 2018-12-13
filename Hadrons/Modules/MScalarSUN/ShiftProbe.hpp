/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MScalarSUN/ShiftProbe.hpp

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
#ifndef Hadrons_MScalarSUN_ShiftProbe_hpp_
#define Hadrons_MScalarSUN_ShiftProbe_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Modules/MScalarSUN/Utils.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *         Ward identity phi^n probe with fields at different positions       *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MScalarSUN)

typedef std::pair<int, int> ShiftPair;

class ShiftProbePar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(ShiftProbePar,
                                    std::string, field,
                                    std::string, shifts,
                                    std::string, output);
};

class ShiftProbeResult: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(ShiftProbeResult,
                                    std::string, shifts,
                                    Complex,     value);
};

template <typename SImpl>
class TShiftProbe: public Module<ShiftProbePar>
{
public:
    typedef typename SImpl::Field                          Field;
    typedef typename SImpl::ComplexField                   ComplexField;
public:
    // constructor
    TShiftProbe(const std::string name);
    // destructor
    virtual ~TShiftProbe(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(ShiftProbeSU2, TShiftProbe<ScalarNxNAdjImplR<2>>, MScalarSUN);
MODULE_REGISTER_TMP(ShiftProbeSU3, TShiftProbe<ScalarNxNAdjImplR<3>>, MScalarSUN);
MODULE_REGISTER_TMP(ShiftProbeSU4, TShiftProbe<ScalarNxNAdjImplR<4>>, MScalarSUN);
MODULE_REGISTER_TMP(ShiftProbeSU5, TShiftProbe<ScalarNxNAdjImplR<5>>, MScalarSUN);
MODULE_REGISTER_TMP(ShiftProbeSU6, TShiftProbe<ScalarNxNAdjImplR<6>>, MScalarSUN);

/******************************************************************************
 *                        TShiftProbe implementation                          *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename SImpl>
TShiftProbe<SImpl>::TShiftProbe(const std::string name)
: Module<ShiftProbePar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename SImpl>
std::vector<std::string> TShiftProbe<SImpl>::getInput(void)
{
    std::vector<std::string> in = {par().field};
    
    return in;
}

template <typename SImpl>
std::vector<std::string> TShiftProbe<SImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename SImpl>
void TShiftProbe<SImpl>::setup(void)
{
    envTmpLat(Field, "acc");
    envCreateLat(ComplexField, getName());
}

// execution ///////////////////////////////////////////////////////////////////
template <typename SImpl>
void TShiftProbe<SImpl>::execute(void)
{
    LOG(Message) << "Creating shift probe for shifts " << par().shifts
                 << std::endl;

    std::vector<ShiftPair> shift;
    double                 sign;
    auto                   &phi   = envGet(Field, par().field);
    auto                   &probe = envGet(ComplexField, getName());

    shift = strToVec<ShiftPair>(par().shifts);
    if (shift.size() % 2 != 0)
    {
        HADRONS_ERROR(Size, "the number of shifts is odd");
    }
    sign = (shift.size() % 4 == 0) ? 1. : -1.;
    for (auto &s: shift)
    {
        if (s.first >= env().getNd())
        {
            HADRONS_ERROR(Size, "dimension to large for shift <" 
                               + std::to_string(s.first) + " " 
                               + std::to_string(s.second) + ">" );
        }
    }
    envGetTmp(Field, acc);
    acc = 1.;
    for (unsigned int i = 0; i < shift.size(); ++i)
    {
        if (shift[i].second == 0)
        {
            acc *= phi;
        }
        else
        {
            acc *= Cshift(phi, shift[i].first, shift[i].second);
        }
    }
    probe = sign*trace(acc);
    if (!par().output.empty())
    {
        ShiftProbeResult r;

        r.shifts = par().shifts;
        r.value  = TensorRemove(sum(probe));
        saveResult(par().output, "probe", r);
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MScalarSUN_ShiftProbe_hpp_
