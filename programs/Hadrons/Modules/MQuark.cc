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

#include <Hadrons/Modules/MQuark.hpp>

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

// dependencies/products ///////////////////////////////////////////////////////
std::vector<std::string> MQuark::getInput(void)
{
    std::vector<std::string> in = {par().source, par().solver};
    
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
    Ls_ = env().getObjectLs(par().solver);
    env().registerLattice<LatticePropagator>(getName());
    if (Ls_ > 1)
    {
        env().registerLattice<LatticePropagator>(getName() + "_5d", Ls_);
    }
}

// execution ///////////////////////////////////////////////////////////////////
void MQuark::execute(void)
{
    LatticeFermion source(env().getGrid(Ls_)), sol(env().getGrid(Ls_)),
                   tmp(env().getGrid());
    std::string    propName = (Ls_ == 1) ? getName() : (getName() + "_5d");
    
    LOG(Message) << "Computing quark propagator '" << getName() << "'"
                 << std::endl;
    LatticePropagator   &prop    = *env().createLattice<LatticePropagator>(propName);
    LatticePropagator   &fullSrc = *env().getObject<LatticePropagator>(par().source);
    Solver              &solver  = *env().getObject<Solver>(par().solver);
    if (Ls_ > 1)
    {
        env().createLattice<LatticePropagator>(getName());
    }

    LOG(Message) << "Inverting using solver '" << par().solver
                 << "' on source '" << par().source << "'" << std::endl;
    for (unsigned int s = 0; s < Ns; ++s)
    for (unsigned int c = 0; c < Nc; ++c)
    {
        LOG(Message) << "Inversion for spin= " << s << ", color= " << c
                     << std::endl;
        // source conversion for 4D sources
        if (!env().isObject5d(par().source))
        {
            if (Ls_ == 1)
            {
                PropToFerm(source, fullSrc, s, c);
            }
            else
            {
                source = zero;
                PropToFerm(tmp, fullSrc, s, c);
                InsertSlice(tmp, source, 0, 0);
                InsertSlice(tmp, source, Ls_-1, 0);
                axpby_ssp_pplus(source, 0., source, 1., source, 0, 0);
                axpby_ssp_pminus(source, 0., source, 1., source, Ls_-1, Ls_-1);
            }
        }
        // source conversion for 5D sources
        else
        {
            if (Ls_ != env().getObjectLs(par().source))
            {
                HADRON_ERROR("Ls mismatch between quark action and source");
            }
            else
            {
                PropToFerm(source, fullSrc, s, c);
            }
        }
        sol = zero;
        solver(sol, source);
        FermToProp(prop, sol, s, c);
        // create 4D propagators from 5D one if necessary
        if (Ls_ > 1)
        {
            LatticePropagator &p4d = *env().getObject<LatticePropagator>(getName());
            
            axpby_ssp_pminus(sol, 0., sol, 1., sol, 0, 0);
            axpby_ssp_pplus(sol, 0., sol, 1., sol, 0, Ls_-1);
            ExtractSlice(tmp, sol, 0, 0);
            FermToProp(p4d, tmp, s, c);
        }
    }
}
