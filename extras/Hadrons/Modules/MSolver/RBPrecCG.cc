/*******************************************************************************
Grid physics library, www.github.com/paboyle/Grid 

Source file: programs/Hadrons/RBPrecCG.cc

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

#include <Grid/Hadrons/Modules/MSolver/RBPrecCG.hpp>

using namespace Grid;
using namespace QCD;
using namespace Hadrons;
using namespace MSolver;

/******************************************************************************
*                          RBPrecCG implementation                            *
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
RBPrecCG::RBPrecCG(const std::string name)
: Module(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
std::vector<std::string> RBPrecCG::getInput(void)
{
    std::vector<std::string> in = {par().action};
    
    return in;
}

std::vector<std::string> RBPrecCG::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
void RBPrecCG::setup(void)
{
    auto Ls = env().getObjectLs(par().action);
    
    env().registerObject(getName(), 0, Ls);
    env().addOwnership(getName(), par().action);
}

// execution ///////////////////////////////////////////////////////////////////
void RBPrecCG::execute(void)
{
    auto &mat   = *(env().getObject<FMat>(par().action));
    auto solver = [&mat, this](LatticeFermion &sol,
                               const LatticeFermion &source)
    {
        ConjugateGradient<LatticeFermion>           cg(par().residual, 10000);
        SchurRedBlackDiagMooeeSolve<LatticeFermion> schurSolver(cg);
        
        schurSolver(mat, source, sol);
    };
    
    LOG(Message) << "setting up Schur red-black preconditioned CG for"
                 << " action '" << par().action << "' with residual "
                 << par().residual << std::endl;
    env().setObject(getName(), new SolverFn(solver));
}
