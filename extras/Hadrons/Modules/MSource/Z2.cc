/*******************************************************************************
Grid physics library, www.github.com/paboyle/Grid 

Source file: programs/Hadrons/Z2.cc

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

#include <Grid/Hadrons/Modules/MSource/Z2.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MSource;

/******************************************************************************
*                              Z2 implementation                              *
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
Z2::Z2(const std::string name)
: Module<Z2Par>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
std::vector<std::string> Z2::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

std::vector<std::string> Z2::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
void Z2::setup(void)
{
    env().registerLattice<PropagatorField>(getName());
}

// execution ///////////////////////////////////////////////////////////////////
void Z2::execute(void)
{
    Lattice<iScalar<vInteger>> t(env().getGrid());
    LatticeComplex             eta(env().getGrid());
    LatticeFermion             phi(env().getGrid());
    Complex                    shift(1., 1.);
    
    if (par().tA == par().tB)
    {
        LOG(Message) << "Generating Z_2 wall source at t= " << par().tA
                     << std::endl;
    }
    else
    {
        LOG(Message) << "Generating Z_2 band for " << par().tA << " <= t <= "
                     << par().tB << std::endl;
    }
    PropagatorField &src = *env().createLattice<PropagatorField>(getName());
    LatticeCoordinate(t, Tp);
    bernoulli(*env().get4dRng(), eta);
    eta = (2.*eta - shift)*(1./::sqrt(2.));
    eta = where((t >= par().tA) and (t <= par().tB), eta, 0.*eta);
    src = 1.;
    src = src*eta;
}
