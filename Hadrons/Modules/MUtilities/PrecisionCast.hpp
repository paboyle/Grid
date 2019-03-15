/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MUtilities/PrecisionCast.hpp

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
#ifndef Hadrons_MUtilities_PrecisionCast_hpp_
#define Hadrons_MUtilities_PrecisionCast_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                          Precision cast module                             *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)

class PrecisionCastPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(PrecisionCastPar,
                                    std::string, field);
};

template <typename FieldIn, typename FieldOut>
class TPrecisionCast: public Module<PrecisionCastPar>
{
public:
    // constructor
    TPrecisionCast(const std::string name);
    // destructor
    virtual ~TPrecisionCast(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(GaugeSinglePrecisionCast, 
                    ARG(TPrecisionCast<GIMPLD::GaugeField, GIMPLF::GaugeField>),
                    MUtilities);
MODULE_REGISTER_TMP(FermionSinglePrecisionCast, 
                    ARG(TPrecisionCast<FIMPLD::FermionField, FIMPLF::FermionField>),
                    MUtilities);

/******************************************************************************
 *                     TPrecisionCast implementation                          *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FieldIn, typename FieldOut>
TPrecisionCast<FieldIn, FieldOut>::TPrecisionCast(const std::string name)
: Module<PrecisionCastPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FieldIn, typename FieldOut>
std::vector<std::string> TPrecisionCast<FieldIn, FieldOut>::getInput(void)
{
    std::vector<std::string> in = {par().field};
    
    return in;
}

template <typename FieldIn, typename FieldOut>
std::vector<std::string> TPrecisionCast<FieldIn, FieldOut>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FieldIn, typename FieldOut>
void TPrecisionCast<FieldIn, FieldOut>::setup(void)
{
    envCreateLat(FieldOut, getName());
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FieldIn, typename FieldOut>
void TPrecisionCast<FieldIn, FieldOut>::execute(void)
{
    LOG(Message) << "Casting field '" << par().field << "'" << std::endl;
    LOG(Message) << "In  type: " << typeName<FieldIn>() << std::endl;
    LOG(Message) << "Out type: " << typeName<FieldOut>() << std::endl;

    auto &in  = envGet(FieldIn,  par().field);
    auto &out = envGet(FieldOut, getName());

    precisionChange(out, in);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MUtilities_PrecisionCast_hpp_
