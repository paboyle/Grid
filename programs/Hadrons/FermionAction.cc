/*******************************************************************************
Grid physics library, www.github.com/paboyle/Grid 

Source file: programs/Hadrons/FermionAction.cc

Copyright (C) 2016

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

See the full license in the file "LICENSE" in the top level distribution 
directory.
*******************************************************************************/

#include <Hadrons/FermionAction.hpp>

using namespace Grid;
using namespace QCD;
using namespace Hadrons;

/******************************************************************************
 *                      FermionAction implementation                          *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
FermionAction::FermionAction(const std::string name)
: name_(name)
{}

// access //////////////////////////////////////////////////////////////////////
std::string FermionAction::getName(void) const
{
    return name_;
}

unsigned int FermionAction::getLs(void) const
{
    return 1;
}

void FermionAction::setFMat(FMat *fMat)
{
    fMat_.reset(fMat);
}

FermionAction::FMat * FermionAction::getFMat(void)
{
    return fMat_.get();
}
