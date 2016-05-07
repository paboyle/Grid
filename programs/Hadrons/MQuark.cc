/*******************************************************************************
Grid physics library, www.github.com/paboyle/Grid 

Source file: programs/Hadrons/MQuark.cc

Copyright (C) 2015

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

#include <Hadrons/MQuark.hpp>

using namespace Grid;
using namespace QCD;
using namespace Hadrons;

/******************************************************************************
 *                          MQuark implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
MQuark::MQuark(const std::string name)
: Module(name)
{}

// parse parameters ////////////////////////////////////////////////////////////
void MQuark::parseParameters(XmlReader &reader, const std::string name)
{
    read(reader, name, par_);
}

// dependencies/products ///////////////////////////////////////////////////////
std::vector<std::string> MQuark::getInput(void)
{
    std::vector<std::string> in = {par_.source, par_.solver};
    
    return in;
}

std::vector<std::string> MQuark::getOutput(void)
{
    std::vector<std::string> out = {getName(), getName() + "_5d"};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
void MQuark::setup(void)
{
    Ls_ = env().getObjectLs(env().getSolverAction(par_.solver));
    env().registerLattice<LatticePropagator>(getName());
    if (Ls_ > 1)
    {
        env().registerLattice<LatticePropagator>(getName() + "_5d", Ls_);
    }
}

// execution ///////////////////////////////////////////////////////////////////
void MQuark::execute(void)
{
    LatticePropagator *fullSource;
    LatticeFermion    source(env().getGrid(Ls_)), sol(env().getGrid(Ls_));
    std::string       propName = (Ls_ == 1) ? getName() : (getName() + "_5d");
    
    LOG(Message) << "Computing quark propagator '" << getName() << "'"
                 << std::endl;
    LatticePropagator &prop = *env().create<LatticePropagator>(propName);
    // source conversion for 4D sources
    if (!env().isObject5d(par_.source))
    {
        if (Ls_ == 1)
        {
            fullSource = env().get<LatticePropagator>(par_.source);
        }
        else
        {
            HADRON_ERROR("MQuark not implemented with 5D actions");
        }
    }
    // source conversion for 5D sources
    else
    {
        if (Ls_ == 1)
        {
            HADRON_ERROR("MQuark not implemented with 5D actions");
        }
        else if (Ls_ != env().getObjectLs(par_.source))
        {
            HADRON_ERROR("Ls mismatch between quark action and source");
        }
        else
        {
            fullSource = env().get<LatticePropagator>(par_.source);
        }
    }
    LOG(Message) << "Inverting using solver '" << par_.solver
                 << "' on source '" << par_.source << "'" << std::endl;
    for (unsigned int s = 0; s < Ns; ++s)
    for (unsigned int c = 0; c < Nc; ++c)
    {
        PropToFerm(source, *fullSource, s, c);
        sol = zero;
        env().callSolver(par_.solver, sol, source);
        FermToProp(prop, sol, s, c);
    }
    // create 4D propagators from 5D one if necessary
    if (Ls_ > 1)
    {
        HADRON_ERROR("MQuark not implemented with 5D actions");
    }
}
