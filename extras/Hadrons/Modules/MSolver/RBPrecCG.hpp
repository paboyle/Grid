/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: extras/Hadrons/Modules/MSolver/RBPrecCG.hpp

Copyright (C) 2015
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

See the full license in the file "LICENSE" in the top level distribution directory
*************************************************************************************/
/*  END LEGAL */

#ifndef Hadrons_MSolver_RBPrecCG_hpp_
#define Hadrons_MSolver_RBPrecCG_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                     Schur red-black preconditioned CG                      *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSolver)

class RBPrecCGPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(RBPrecCGPar,
                                    std::string, action,
                                    double     , residual);
};

template <typename FImpl>
class TRBPrecCG: public Module<RBPrecCGPar>
{
public:
    FGS_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TRBPrecCG(const std::string name);
    // destructor
    virtual ~TRBPrecCG(void) = default;
    // dependencies/products
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_NS(RBPrecCG, TRBPrecCG<FIMPL>, MSolver);

/******************************************************************************
 *                      TRBPrecCG template implementation                     *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TRBPrecCG<FImpl>::TRBPrecCG(const std::string name)
: Module(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TRBPrecCG<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().action};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TRBPrecCG<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TRBPrecCG<FImpl>::setup(void)
{
    LOG(Message) << "setting up Schur red-black preconditioned CG for"
                 << " action '" << par().action << "' with residual "
                 << par().residual << std::endl;

    auto Ls     = env().getObjectLs(par().action);
    auto &mat   = mGetObj(FMat, par().action);
    auto solver = [&mat, this](FermionField &sol, const FermionField &source)
    {
        ConjugateGradient<FermionField>           cg(par().residual, 10000);
        SchurRedBlackDiagMooeeSolve<FermionField> schurSolver(cg);
        
        schurSolver(mat, source, sol);
    };
    mCreateObj(SolverFn, getName(), Ls, solver);
    env().addOwnership(getName(), par().action);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TRBPrecCG<FImpl>::execute(void)
{}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSolver_RBPrecCG_hpp_
