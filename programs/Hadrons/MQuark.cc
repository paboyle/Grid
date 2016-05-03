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
void MQuark::setup(Environment &env)
{
    auto dim = env.getFermionMatrix(par_.solver)->Grid()->GlobalDimensions();
    
    if (dim.size() == Nd)
    {
        Ls_ = 1;
    }
    else
    {
        Ls_ = dim[0];
    }
}

// allocation //////////////////////////////////////////////////////////////////
void MQuark::allocate(Environment &env)
{
    env.createProp(getName());
    quark_ = env.getProp(getName());
    if (Ls_ > 1)
    {
        env.createProp(getName() + "_5d", Ls_);
        quark5d_ = env.getProp(getName() + "_5d");
    }
}

// execution ///////////////////////////////////////////////////////////////////
void MQuark::execute(Environment &env)
{
    LatticePropagator *fullSource;
    LatticeFermion    source(env.getGrid(Ls_)), sol(env.getGrid(Ls_));
    
    LOG(Message) << "computing quark propagator '" << getName() << "'"
                 << std::endl;
    if (!env.isProp5d(par_.source))
    {
        if (Ls_ == 1)
        {
            fullSource = env.getProp(par_.source);
        }
        else
        {
            HADRON_ERROR("MQuark not implemented with 5D actions");
        }
    }
    else
    {
        if (Ls_ == 1)
        {
            HADRON_ERROR("MQuark not implemented with 5D actions");
        }
        else if (Ls_ != env.getPropLs(par_.source))
        {
            HADRON_ERROR("MQuark not implemented with 5D actions");
        }
        else
        {
            fullSource = env.getProp(par_.source);
        }
    }
    
    LOG(Message) << "inverting using solver '" << par_.solver
                 << "' on source '" << par_.source << "'" << std::endl;
    for (unsigned int s = 0; s < Ns; ++s)
    for (unsigned int c = 0; c < Nc; ++c)
    {
        PropToFerm(source, *fullSource, s, c);
        sol = zero;
        env.callSolver(par_.solver, sol, source);
        if (Ls_ == 1)
        {
            FermToProp(*quark_, sol, s, c);
        }
        else
        {
            HADRON_ERROR("MQuark not implemented with 5D actions");
        }
    }
}
