/*******************************************************************************
Grid physics library, www.github.com/paboyle/Grid 

Source file: programs/Hadrons/SrcZ2.cc

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

#include <Hadrons/SrcZ2.hpp>

using namespace Grid;
using namespace Hadrons;

/******************************************************************************
*                           SrcZ2 implementation                              *
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
SrcZ2::SrcZ2(const std::string name)
: Module(name)
{}

// parse parameters ////////////////////////////////////////////////////////////
void SrcZ2::parseParameters(XmlReader &reader, const std::string name)
{
   read(reader, name, par_);
}

// dependencies/products ///////////////////////////////////////////////////////
std::vector<std::string> SrcZ2::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

std::vector<std::string> SrcZ2::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
void SrcZ2::setup(void)
{
    env().registerLattice<LatticePropagator>(getName());
}

// execution ///////////////////////////////////////////////////////////////////
void SrcZ2::execute(void)
{
    Lattice<iScalar<vInteger>> t(env().getGrid());
    LatticeComplex             eta(env().getGrid());
    LatticeFermion             phi(env().getGrid());
    Complex                    shift(1., 1.);
    
    if (par_.tA == par_.tB)
    {
        LOG(Message) << "Generating Z_2 wall source at t= " << par_.tA
                     << std::endl;
    }
    else
    {
        LOG(Message) << "Generating Z_2 band for " << par_.tA << " <= t <= "
                     << par_.tB << std::endl;
    }
    LatticePropagator &src = *env().create<LatticePropagator>(getName());
    LatticeCoordinate(t, Tp);
    bernoulli(*env().get4dRng(), eta);
    eta = (2.*eta - shift)*(1./::sqrt(2.));
    eta = where((t >= par_.tA) and (t <= par_.tB), eta, 0.*eta);
    src = 1.;
    src = src*eta;
}
