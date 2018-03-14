/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: extras/Hadrons/Modules/MSolver/RBPrecCG.hpp

Copyright (C) 2015-2018

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
#include <Grid/Hadrons/EigenPack.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                     Schur red-black preconditioned CG                      *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSolver)

class RBPrecCGPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(RBPrecCGPar ,
                                    std::string , action,
                                    unsigned int, maxIteration,
                                    double      , residual,
                                    std::string , eigenPack);
};

template <typename FImpl, int nBasis>
class TRBPrecCG: public Module<RBPrecCGPar>
{
public:
    typedef FermionEigenPack<FImpl>                       EPack;
    typedef CoarseFermionEigenPack<FImpl, nBasis>         CoarseEPack;
    typedef DeflatedGuesser<typename FImpl::FermionField> FineGuesser;
    typedef LocalCoherenceDeflatedGuesser<
        typename FImpl::FermionField,
        typename CoarseEPack::CoarseField> CoarseGuesser;
    FGS_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TRBPrecCG(const std::string name);
    // destructor
    virtual ~TRBPrecCG(void) = default;
    // dependencies/products
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getReference(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_NS(RBPrecCG,  
    ARG(TRBPrecCG<FIMPL, HADRONS_DEFAULT_LANCZOS_NBASIS>), MSolver);
MODULE_REGISTER_NS(ZRBPrecCG, 
    ARG(TRBPrecCG<ZFIMPL, HADRONS_DEFAULT_LANCZOS_NBASIS>), MSolver);

/******************************************************************************
 *                      TRBPrecCG template implementation                     *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, int nBasis>
TRBPrecCG<FImpl, nBasis>::TRBPrecCG(const std::string name)
: Module(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, int nBasis>
std::vector<std::string> TRBPrecCG<FImpl, nBasis>::getInput(void)
{
    std::vector<std::string> in = {};
    
    return in;
}

template <typename FImpl, int nBasis>
std::vector<std::string> TRBPrecCG<FImpl, nBasis>::getReference(void)
{
    std::vector<std::string> ref = {par().action};
    
    return ref;
}

template <typename FImpl, int nBasis>
std::vector<std::string> TRBPrecCG<FImpl, nBasis>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, int nBasis>
void TRBPrecCG<FImpl, nBasis>::setup(void)
{
    if (par().maxIteration == 0)
    {
        HADRON_ERROR(Argument, "zero maximum iteration");
    }

    LOG(Message) << "setting up Schur red-black preconditioned CG for"
                 << " action '" << par().action << "' with residual "
                 << par().residual << ", maximum iteration " 
                 << par().maxIteration << std::endl;

    auto Ls                 = env().getObjectLs(par().action);
    auto &mat               = envGet(FMat, par().action);
    std::string guesserName = getName() + "_guesser";

    if (par().eigenPack.empty())
    {
        env().template createDerivedObject<Guesser<FermionField>, ZeroGuesser<FermionField>>
            (guesserName, Environment::Storage::object, Ls);
    }
    else
    {
        try
        {
            auto &epack = envGetDerived(EPack, CoarseEPack, par().eigenPack);

            envCreateDerived(Guesser<FermionField>, CoarseGuesser,
                             guesserName, Ls, epack.evec, epack.evecCoarse,
                             epack.evalCoarse);
        }
        catch (Exceptions::Definition &e)
        {
            auto &epack = envGet(EPack, par().eigenPack);

            envCreateDerived(Guesser<FermionField>, FineGuesser,
                             guesserName, Ls, epack.evec, epack.eval);
        }
    }
    auto &guesser = envGet(Guesser<FermionField>, guesserName);
    auto solver = [&mat, &guesser, this](FermionField &sol, const FermionField &source)
    {
        ConjugateGradient<FermionField>           cg(par().residual, 
                                                     par().maxIteration);
        HADRONS_DEFAULT_SCHUR_SOLVE<FermionField> schurSolver(cg);
        
        schurSolver(mat, source, sol, guesser);
    };
    envCreate(SolverFn, getName(), Ls, solver);
}


// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, int nBasis>
void TRBPrecCG<FImpl, nBasis>::execute(void)
{}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSolver_RBPrecCG_hpp_
