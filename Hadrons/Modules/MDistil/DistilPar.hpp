/*************************************************************************************
 
 Grid physics library, www.github.com/paboyle/Grid
 
 Source file: Hadrons/Modules/MDistil/DistilPar.hpp
 
 Copyright (C) 2019
 
 Author: Felix Erben <ferben@ed.ac.uk>
 Author: Michael Marshall <Michael.Marshall@ed.ac.uk>
 
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

#ifndef Hadrons_MDistil_DistilPar_hpp_
#define Hadrons_MDistil_DistilPar_hpp_

#include <Hadrons/Modules/MDistil/Distil.hpp>

BEGIN_HADRONS_NAMESPACE
BEGIN_MODULE_NAMESPACE(MDistil)

/******************************************************************************
 *                         DistilPar                                 *
 ******************************************************************************/

template <typename FImpl>
class TDistilPar: public Module<DistilParameters>
{
public:
    // constructor
    TDistilPar(const std::string name);
    // destructor
    virtual ~TDistilPar(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(DistilPar, TDistilPar<FIMPL>, MDistil);

/******************************************************************************
 *                 TDistilPar implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TDistilPar<FImpl>::TDistilPar(const std::string name) : Module<DistilParameters>(name) {}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TDistilPar<FImpl>::getInput(void)
{
    return {};
}

template <typename FImpl>
std::vector<std::string> TDistilPar<FImpl>::getOutput(void)
{
    return {getName()};
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDistilPar<FImpl>::setup(void)
{
    envCreate(DistilParameters, getName(), 1, par() );
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDistilPar<FImpl>::execute(void)
{
    // Nothing to do. setup() created and initialised the output object
}

END_MODULE_NAMESPACE
END_HADRONS_NAMESPACE
#endif
