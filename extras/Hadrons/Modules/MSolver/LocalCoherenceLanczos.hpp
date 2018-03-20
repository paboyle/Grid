/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: extras/Hadrons/Modules/MSolver/LocalCoherenceLanczos.hpp

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
#ifndef Hadrons_MSolver_LocalCoherenceLanczos_hpp_
#define Hadrons_MSolver_LocalCoherenceLanczos_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>
#include <Grid/Hadrons/LanczosUtils.hpp>

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
                                    int,           doFine,
                                    int,           doCoarse,
                                    LanczosParams, fineParams,
                                    LanczosParams, coarseParams,
                                    ChebyParams,   smoother,
                                    RealD,         coarseRelaxTol,
                                    std::string,   blockSize,
                                    std::string,   output);
};

template <typename FImpl, int nBasis>
class TLocalCoherenceLanczos: public Module<LocalCoherenceLanczosPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    typedef LocalCoherenceLanczos<typename FImpl::SiteSpinor, 
                                  typename FImpl::SiteComplex, 
                                  nBasis>                LCL;
    typedef FineEigenPack<FImpl>                         FinePack;
    typedef CoarseEigenPack<FImpl, nBasis>               CoarsePack; 
    typedef HADRONS_DEFAULT_SCHUR_OP<FMat, FermionField> SchurFMat;
public:
    // constructor
    TLocalCoherenceLanczos(const std::string name);
    // destructor
    virtual ~TLocalCoherenceLanczos(void) = default;
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    void makeCoarseGrid(void);
private:
    std::vector<int>                       coarseDim_;
    int                                    Ls_, cLs_{1};
    std::unique_ptr<GridCartesian>         coarseGrid4_{nullptr};
    std::unique_ptr<GridCartesian>         coarseGrid_{nullptr};
    std::unique_ptr<GridRedBlackCartesian> coarseGrid4Rb_{nullptr};
    std::unique_ptr<GridRedBlackCartesian> coarseGridRb_{nullptr};
    std::string                            fineName_, coarseName_;
};

MODULE_REGISTER_NS(LocalCoherenceLanczos, 
    ARG(TLocalCoherenceLanczos<FIMPL, HADRONS_DEFAULT_LANCZOS_NBASIS>), 
    MSolver);
MODULE_REGISTER_NS(ZLocalCoherenceLanczos, 
    ARG(TLocalCoherenceLanczos<ZFIMPL, HADRONS_DEFAULT_LANCZOS_NBASIS>), 
    MSolver);

/******************************************************************************
 *                 TLocalCoherenceLanczos implementation                      *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, int nBasis>
TLocalCoherenceLanczos<FImpl, nBasis>::TLocalCoherenceLanczos(const std::string name)
: Module<LocalCoherenceLanczosPar>(name)
{
    fineName_   = getName() + "_fine";
    coarseName_ = getName() + "_coarse";
}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, int nBasis>
std::vector<std::string> TLocalCoherenceLanczos<FImpl, nBasis>::getInput(void)
{
    std::vector<std::string> in = {par().action};
    
    return in;
}

template <typename FImpl, int nBasis>
std::vector<std::string> TLocalCoherenceLanczos<FImpl, nBasis>::getOutput(void)
{
    std::vector<std::string> out = {fineName_, coarseName_};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, int nBasis>
void TLocalCoherenceLanczos<FImpl, nBasis>::makeCoarseGrid(void)
{
    int              nd        = env().getNd();
    std::vector<int> blockSize = strToVec<int>(par().blockSize);
    auto             fineDim   = env().getDim();

    Ls_ = env().getObjectLs(par().action);
    env().createGrid(Ls_);
    coarseDim_.resize(nd);
    for (int d = 0; d < coarseDim_.size(); d++)
    {
        coarseDim_[d] = fineDim[d]/blockSize[d];
        if (coarseDim_[d]*blockSize[d] != fineDim[d])
        {
            HADRON_ERROR(Size, "Fine dimension " + std::to_string(d) 
                         + " (" + std::to_string(fineDim[d]) 
                         + ") not divisible by coarse dimension ("
                         + std::to_string(coarseDim_[d]) + ")"); 
        }
    }
    if (blockSize.size() > nd)
    {
        cLs_ = Ls_/blockSize[nd];
        if (cLs_*blockSize[nd] != Ls_)
        {
            HADRON_ERROR(Size, "Fine Ls (" + std::to_string(Ls_) 
                         + ") not divisible by coarse Ls ("
                         + std::to_string(cLs_) + ")");
        }
    }
    if (Ls_ > 1)
    {
        coarseGrid4_.reset(SpaceTimeGrid::makeFourDimGrid(
            coarseDim_, GridDefaultSimd(nd, vComplex::Nsimd()),
            GridDefaultMpi()));
        coarseGrid4Rb_.reset(SpaceTimeGrid::makeFourDimRedBlackGrid(coarseGrid4_.get()));
        coarseGrid_.reset(SpaceTimeGrid::makeFiveDimGrid(cLs_, coarseGrid4_.get()));
        coarseGridRb_.reset(SpaceTimeGrid::makeFiveDimRedBlackGrid(cLs_, coarseGrid4_.get()));
    }
    else
    {
        coarseGrid_.reset(SpaceTimeGrid::makeFourDimGrid(
            coarseDim_, GridDefaultSimd(nd, vComplex::Nsimd()),
            GridDefaultMpi()));
        coarseGridRb_.reset(SpaceTimeGrid::makeFourDimRedBlackGrid(coarseGrid_.get()));
    }
}

template <typename FImpl, int nBasis>
void TLocalCoherenceLanczos<FImpl, nBasis>::setup(void)
{
    LOG(Message) << "Setting up local coherence Lanczos eigensolver for"
                 << " action '" << par().action << "' (" << nBasis
                 << " eigenvectors)..." << std::endl;
    
    if (!coarseGrid_)
    {
        makeCoarseGrid();
    }
    LOG(Message) << "Coarse grid: " << coarseGrid_->GlobalDimensions() << std::endl;
    envCreate(FinePack, fineName_, Ls_, par().fineParams.Nm, env().getRbGrid(Ls_));
    envCreate(CoarsePack, coarseName_, Ls_, par().coarseParams.Nm, coarseGridRb_.get());
    auto &fine   = envGet(FinePack, fineName_);
    auto &coarse = envGet(CoarsePack, coarseName_);
    envTmp(SchurFMat, "mat", Ls_, envGet(FMat, par().action));
    envGetTmp(SchurFMat, mat);
    envTmp(LCL, "solver", Ls_, env().getRbGrid(Ls_), coarseGridRb_.get(), mat, 
           Odd, fine.evec, coarse.evec, fine.eval, coarse.eval);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, int nBasis>
void TLocalCoherenceLanczos<FImpl, nBasis>::execute(void)
{
    auto &finePar   = par().fineParams;
    auto &coarsePar = par().coarseParams;
    auto &fine      = envGet(FinePack, fineName_);
    auto &coarse    = envGet(CoarsePack, coarseName_);

    envGetTmp(LCL, solver);
    if (par().doFine)
    {
        LOG(Message) << "Performing fine grid IRL -- Nstop= " 
                     << finePar.Nstop << ", Nk= " << finePar.Nk << ", Nm= " 
                     << finePar.Nm << std::endl;
        solver.calcFine(finePar.Cheby, finePar.Nstop, finePar.Nk, finePar.Nm,
                        finePar.resid,finePar.MaxIt, finePar.betastp, 
                        finePar.MinRes);
        solver.testFine(finePar.resid*100.0);
        LOG(Message) << "Orthogonalising" << std::endl;
        solver.Orthogonalise();
        if (!par().output.empty())
        {
            fine.write(par().output + "_fine");
        }
    }
    if (par().doCoarse)
    {
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
            coarse.write(par().output + "_coarse");
        }
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSolver_LocalCoherenceLanczos_hpp_
