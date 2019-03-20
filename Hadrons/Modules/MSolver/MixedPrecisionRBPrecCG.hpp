/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MSolver/MixedPrecisionRBPrecCG.hpp

Copyright (C) 2015-2019

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
#ifndef Hadrons_MSolver_MixedPrecisionRBPrecCG_hpp_
#define Hadrons_MSolver_MixedPrecisionRBPrecCG_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Solver.hpp>
#include <Hadrons/EigenPack.hpp>
#include <Hadrons/Modules/MSolver/Guesser.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *              Mixed precision schur red-black preconditioned CG             *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSolver)

class MixedPrecisionRBPrecCGPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(MixedPrecisionRBPrecCGPar,
                                    std::string , innerAction,
                                    std::string , outerAction,
                                    unsigned int, maxInnerIteration,
                                    unsigned int, maxOuterIteration,
                                    double      , residual,
                                    std::string , eigenPack);
};

template <typename FImplInner, typename FImplOuter, int nBasis>
class TMixedPrecisionRBPrecCG: public Module<MixedPrecisionRBPrecCGPar>
{
public:
    FERM_TYPE_ALIASES(FImplInner, Inner);
    FERM_TYPE_ALIASES(FImplOuter, Outer);
    SOLVER_TYPE_ALIASES(FImplOuter,);
    typedef HADRONS_DEFAULT_SCHUR_OP<FMatInner, FermionFieldInner> SchurFMatInner;
    typedef HADRONS_DEFAULT_SCHUR_OP<FMatOuter, FermionFieldOuter> SchurFMatOuter;
private:
    template <typename Field>
    class OperatorFunctionWrapper: public OperatorFunction<Field>
    {
    public:
        OperatorFunctionWrapper(LinearFunction<Field> &fn): fn_(fn) {};
        virtual ~OperatorFunctionWrapper(void) = default;
        virtual void operator()(LinearOperatorBase<Field> &op, 
                                const Field &in, Field &out)
        {
            fn_(in, out);
        }
    private:
        LinearFunction<Field> &fn_;
    };
public:
    // constructor
    TMixedPrecisionRBPrecCG(const std::string name);
    // destructor
    virtual ~TMixedPrecisionRBPrecCG(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getReference(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(MixedPrecisionRBPrecCG, 
    ARG(TMixedPrecisionRBPrecCG<FIMPLF, FIMPLD, HADRONS_DEFAULT_LANCZOS_NBASIS>), MSolver);
MODULE_REGISTER_TMP(ZMixedPrecisionRBPrecCG, 
    ARG(TMixedPrecisionRBPrecCG<ZFIMPLF, ZFIMPLD, HADRONS_DEFAULT_LANCZOS_NBASIS>), MSolver);

/******************************************************************************
 *                 TMixedPrecisionRBPrecCG implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImplInner, typename FImplOuter, int nBasis>
TMixedPrecisionRBPrecCG<FImplInner, FImplOuter, nBasis>
::TMixedPrecisionRBPrecCG(const std::string name)
: Module<MixedPrecisionRBPrecCGPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImplInner, typename FImplOuter, int nBasis>
std::vector<std::string> TMixedPrecisionRBPrecCG<FImplInner, FImplOuter, nBasis>
::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename FImplInner, typename FImplOuter, int nBasis>
std::vector<std::string> TMixedPrecisionRBPrecCG<FImplInner, FImplOuter, nBasis>
::getReference(void)
{
    std::vector<std::string> ref = {par().innerAction, par().outerAction};
    
    if (!par().eigenPack.empty())
    {
        ref.push_back(par().eigenPack);
    }
    
    return ref;
}

template <typename FImplInner, typename FImplOuter, int nBasis>
std::vector<std::string> TMixedPrecisionRBPrecCG<FImplInner, FImplOuter, nBasis>
::getOutput(void)
{
    std::vector<std::string> out = {getName(), getName() + "_subtract"};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImplInner, typename FImplOuter, int nBasis>
void TMixedPrecisionRBPrecCG<FImplInner, FImplOuter, nBasis>
::setup(void)
{
    LOG(Message) << "Setting up Schur red-black preconditioned mixed-precision "
                 << "CG for inner/outer action '" << par().innerAction 
                 << "'/'" << par().outerAction << "', residual "
                 << par().residual << ", and maximum inner/outer iteration " 
                 << par().maxInnerIteration << "/" << par().maxOuterIteration
                 << std::endl;

    auto Ls        = env().getObjectLs(par().innerAction);
    auto &imat     = envGet(FMatInner, par().innerAction);
    auto &omat     = envGet(FMatOuter, par().outerAction);
    auto guesserPt = makeGuesser<FImplOuter, nBasis>(par().eigenPack);

    auto makeSolver = [&imat, &omat, guesserPt, Ls, this](bool subGuess) 
    {
        return [&imat, &omat, guesserPt, subGuess, Ls, this]
        (FermionFieldOuter &sol, const FermionFieldOuter &source) 
        {
            typedef typename FermionFieldInner::vector_type VTypeInner;

            SchurFMatInner simat(imat);
            SchurFMatOuter somat(omat);
            MixedPrecisionConjugateGradient<FermionFieldOuter, FermionFieldInner> 
                mpcg(par().residual, par().maxInnerIteration, 
                     par().maxOuterIteration, 
                     env().template getRbGrid<VTypeInner>(Ls),
                     simat, somat);
            OperatorFunctionWrapper<FermionFieldOuter> wmpcg(mpcg);
            HADRONS_DEFAULT_SCHUR_SOLVE<FermionFieldOuter> schurSolver(wmpcg);
            schurSolver.subtractGuess(subGuess);
            schurSolver(omat, source, sol, *guesserPt);
        };
    };
    auto solver = makeSolver(false);
    envCreate(Solver, getName(), Ls, solver, omat);
    auto solver_subtract = makeSolver(true);
    envCreate(Solver, getName() + "_subtract", Ls, solver_subtract, omat);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImplInner, typename FImplOuter, int nBasis>
void TMixedPrecisionRBPrecCG<FImplInner, FImplOuter, nBasis>
::execute(void)
{}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSolver_MixedPrecisionRBPrecCG_hpp_
