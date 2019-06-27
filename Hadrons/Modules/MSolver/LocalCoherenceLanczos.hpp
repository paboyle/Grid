/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MSolver/LocalCoherenceLanczos.hpp

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
#ifndef Hadrons_MSolver_LocalCoherenceLanczos_hpp_
#define Hadrons_MSolver_LocalCoherenceLanczos_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/EigenPack.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                    Local coherence Lanczos eigensolver                     *
 *****************************************************************************/
BEGIN_MODULE_NAMESPACE(MSolver)

class LocalCoherenceLanczosPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(LocalCoherenceLanczosPar,
                                    std::string,   action,
                                    bool,          doCoarse,
                                    LanczosParams, fineParams,
                                    LanczosParams, coarseParams,
                                    ChebyParams,   smoother,
                                    RealD,         coarseRelaxTol,
                                    std::string,   blockSize,
                                    std::string,   output,
                                    bool,          multiFile);
};

template <typename FImpl, int nBasis, typename FImplIo = FImpl>
class TLocalCoherenceLanczos: public Module<LocalCoherenceLanczosPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    typedef LocalCoherenceLanczos<typename FImpl::SiteSpinor, 
                                  typename FImpl::SiteComplex, 
                                  nBasis>                  LCL;
    typedef BaseFermionEigenPack<FImpl>                    BasePack;
    typedef CoarseFermionEigenPack<FImpl, nBasis, FImplIo> CoarsePack;
    typedef typename CoarsePack::Field                     Field;
    typedef typename CoarsePack::FieldIo                   FieldIo;
    typedef typename CoarsePack::CoarseField               CoarseField;
    typedef typename CoarsePack::CoarseFieldIo             CoarseFieldIo;
    typedef HADRONS_DEFAULT_SCHUR_OP<FMat, FermionField>   SchurFMat;
public:
    // constructor
    TLocalCoherenceLanczos(const std::string name);
    // destructor
    virtual ~TLocalCoherenceLanczos(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(LocalCoherenceLanczos, ARG(TLocalCoherenceLanczos<FIMPL, HADRONS_DEFAULT_LANCZOS_NBASIS>), MSolver);
MODULE_REGISTER_TMP(ZLocalCoherenceLanczos, ARG(TLocalCoherenceLanczos<ZFIMPL, HADRONS_DEFAULT_LANCZOS_NBASIS>), MSolver);
MODULE_REGISTER_TMP(StagLocalCoherenceLanczos, ARG(TLocalCoherenceLanczos<STAGIMPL, HADRONS_DEFAULT_LANCZOS_NBASIS>), MSolver);
#ifdef GRID_DEFAULT_PRECISION_DOUBLE
MODULE_REGISTER_TMP(LocalCoherenceLanczosIo32, ARG(TLocalCoherenceLanczos<FIMPL, HADRONS_DEFAULT_LANCZOS_NBASIS, FIMPLF>), MSolver);
MODULE_REGISTER_TMP(ZLocalCoherenceLanczosIo32, ARG(TLocalCoherenceLanczos<ZFIMPL, HADRONS_DEFAULT_LANCZOS_NBASIS, ZFIMPLF>), MSolver);
MODULE_REGISTER_TMP(StagLocalCoherenceLanczosIo32, ARG(TLocalCoherenceLanczos<STAGIMPL, HADRONS_DEFAULT_LANCZOS_NBASIS, STAGIMPLF>), MSolver);
#endif

/******************************************************************************
 *                 TLocalCoherenceLanczos implementation                      *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, int nBasis, typename FImplIo>
TLocalCoherenceLanczos<FImpl, nBasis, FImplIo>::TLocalCoherenceLanczos(const std::string name)
: Module<LocalCoherenceLanczosPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, int nBasis, typename FImplIo>
std::vector<std::string> TLocalCoherenceLanczos<FImpl, nBasis, FImplIo>::getInput(void)
{
    std::vector<std::string> in = {par().action};
    
    return in;
}

template <typename FImpl, int nBasis, typename FImplIo>
std::vector<std::string> TLocalCoherenceLanczos<FImpl, nBasis, FImplIo>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, int nBasis, typename FImplIo>
void TLocalCoherenceLanczos<FImpl, nBasis, FImplIo>::setup(void)
{
    LOG(Message) << "Setting up local coherence Lanczos eigensolver for"
    << " action '" << par().action << "' (" << nBasis
    << " eigenvectors)..." << std::endl;
    
    unsigned int Ls        = env().getObjectLs(par().action);
    auto         blockSize = strToVec<int>(par().blockSize);
    GridBase     *grid = nullptr, *gridCoarse = nullptr,
                 *gridIo = nullptr, *gridCoarseIo = nullptr;
    
    if (Ls == 1)
    {
        grid = envGetRbGrid(Field);
        gridCoarse = envGetCoarseGrid(CoarseField, blockSize);
        if (typeHash<Field>() != typeHash<FieldIo>())
        {
            gridIo       = envGetRbGrid(FieldIo);
        }
        if (typeHash<CoarseField>() != typeHash<CoarseFieldIo>())
        {
            gridCoarseIo = envGetCoarseGrid(CoarseFieldIo, blockSize);
        }
    }
    else
    {
        grid = envGetRbGrid(Field, Ls);
        gridCoarse = envGetCoarseGrid(CoarseField, blockSize, Ls);
        if (typeHash<Field>() != typeHash<FieldIo>())
        {
            gridIo       = envGetRbGrid(FieldIo, Ls);
        }
        if (typeHash<CoarseField>() != typeHash<CoarseFieldIo>())
        {
            gridCoarseIo = envGetCoarseGrid(CoarseFieldIo, blockSize, Ls);
        }
    }
    
    int  cNm = (par().doCoarse) ? par().coarseParams.Nm : 0;
    
    LOG(Message) << "Coarse grid: " << gridCoarse->GlobalDimensions() << std::endl;
    envCreateDerived(BasePack, CoarsePack, getName(), Ls,
                     par().fineParams.Nm, cNm, grid, gridCoarse,
                     gridIo, gridCoarseIo);
    
    auto &epack = envGetDerived(BasePack, CoarsePack, getName());
    
    envTmp(SchurFMat, "mat", Ls, envGet(FMat, par().action));
    envGetTmp(SchurFMat, mat);
    envTmp(LCL, "solver", Ls, grid, gridCoarse, mat,
           Odd, epack.evec, epack.evecCoarse, epack.eval, epack.evalCoarse);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, int nBasis, typename FImplIo>
void TLocalCoherenceLanczos<FImpl, nBasis, FImplIo>::execute(void)
{
    auto &finePar   = par().fineParams;
    auto &coarsePar = par().coarseParams;
    auto &epack     = envGetDerived(BasePack, CoarsePack, getName());

    epack.record.operatorXml = vm().getModule(par().action)->parString();
    epack.record.solverXml   = parString();
    envGetTmp(LCL, solver);
    LOG(Message) << "Performing fine grid IRL -- Nstop= " 
                 << finePar.Nstop << ", Nk= " << finePar.Nk << ", Nm= " 
                 << finePar.Nm << std::endl;
    solver.calcFine(finePar.Cheby, finePar.Nstop, finePar.Nk, finePar.Nm,
                    finePar.resid,finePar.MaxIt, finePar.betastp, 
                    finePar.MinRes);
    solver.testFine(finePar.resid*100.0);
    if (!par().output.empty())
    {
        epack.writeFine(par().output, par().multiFile, vm().getTrajectory());
    }
    if (par().doCoarse)
    {
        LOG(Message) << "Orthogonalising" << std::endl;
        solver.Orthogonalise();
        LOG(Message) << "Performing coarse grid IRL -- Nstop= " 
                    << coarsePar.Nstop << ", Nk= " << coarsePar.Nk << ", Nm= " 
                    << coarsePar.Nm << std::endl;
        solver.calcCoarse(coarsePar.Cheby, par().smoother, par().coarseRelaxTol,
                          coarsePar.Nstop, coarsePar.Nk, coarsePar.Nm, 
                          coarsePar.resid, coarsePar.MaxIt, coarsePar.betastp, 
                          coarsePar.MinRes);
        solver.testCoarse(coarsePar.resid*100.0, par().smoother, 
                        par().coarseRelaxTol);
        if (!par().output.empty())
        {
            epack.writeCoarse(par().output, par().multiFile, vm().getTrajectory());
        }
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSolver_LocalCoherenceLanczos_hpp_
