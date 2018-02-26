#ifndef Hadrons_MSolver_LocalCoherenceLanczos_hpp_
#define Hadrons_MSolver_LocalCoherenceLanczos_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>
#include <Grid/algorithms/iterative/LocalCoherenceLanczos.h>

#ifndef DEFAULT_LANCZOS_NBASIS
#define DEFAULT_LANCZOS_NBASIS 60
#endif

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         LocalCoherenceLanczos                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSolver)

class LocalCoherenceLanczosPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(LocalCoherenceLanczosPar,
                                    std::string,      action,
                                    int,              doFine,
                                    int,              doCoarse,
                                    LanczosParams,    fineParams,
                                    LanczosParams,    coarseParams,
                                    ChebyParams,      smoother,
                                    RealD,            coarseRelaxTol,
                                    std::string,      blockSize);
};

template <typename FImpl, int nBasis>
class TLocalCoherenceLanczos: public Module<LocalCoherenceLanczosPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    typedef LocalCoherenceLanczos<vSpinColourVector, vTComplex, nBasis> LCL;
    typedef typename LCL::FineField                             FineField;
    typedef typename LCL::CoarseField                           CoarseField;
    typedef SchurDiagMooeeOperator<FMat, FermionField> SchurFMat;
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
    std::string fevecName_, cevecName_, fevalName_, cevalName_;
};

MODULE_REGISTER_NS(LocalCoherenceLanczos, 
                   ARG(TLocalCoherenceLanczos<FIMPL, DEFAULT_LANCZOS_NBASIS>), 
                   MSolver);

/******************************************************************************
 *                 TLocalCoherenceLanczos implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, int nBasis>
TLocalCoherenceLanczos<FImpl, nBasis>::TLocalCoherenceLanczos(const std::string name)
: Module<LocalCoherenceLanczosPar>(name)
{
    fevecName_ = getName() + "_fineEvec";
    cevecName_ = getName() + "_coarseEvec";
    fevalName_ = getName() + "_fineEval";
    cevalName_ = getName() + "_coarseEval";
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
    std::vector<std::string> out = {fevecName_, cevecName_, fevalName_, cevalName_};
    
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
    envCreate(std::vector<FineField>, fevecName_, Ls_, par().fineParams.Nm,
              env().getRbGrid(Ls_));
    envCreate(std::vector<CoarseField>, cevecName_, Ls_, par().coarseParams.Nm, 
              coarseGridRb_.get());
    envCreate(std::vector<RealD>, fevalName_, Ls_, par().fineParams.Nm);
    envCreate(std::vector<RealD>, cevalName_, Ls_, par().coarseParams.Nm);
    envTmp(SchurFMat, "mat", Ls_, envGet(FMat, par().action));
    envGetTmp(SchurFMat, mat);
    envTmp(LCL, "solver", Ls_, env().getRbGrid(Ls_), coarseGridRb_.get(), mat, Odd, 
           envGet(std::vector<FineField>, fevecName_),
           envGet(std::vector<CoarseField>, cevecName_),
           envGet(std::vector<RealD>, fevalName_),
           envGet(std::vector<RealD>, cevalName_));
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, int nBasis>
void TLocalCoherenceLanczos<FImpl, nBasis>::execute(void)
{
    auto &fine   = par().fineParams;
    auto &coarse = par().coarseParams;
    
    envGetTmp(LCL, solver);
    if (par().doFine)
    {
        LOG(Message) << "Performing fine grid IRL -- Nstop= " 
                     << fine.Nstop << ", Nk= " << fine.Nk << ", Nm= " 
                     << fine.Nm << std::endl;
        solver.calcFine(fine.Cheby, fine.Nstop, fine.Nk, fine.Nm, fine.resid,
                        fine.MaxIt, fine.betastp, fine.MinRes);
        LOG(Message) << "Orthogonalising" << std::endl;
        solver.Orthogonalise();
    }
    if (par().doCoarse)
    {
        LOG(Message) << "Performing coarse grid IRL -- Nstop= " 
                     << fine.Nstop << ", Nk= " << fine.Nk << ", Nm= " 
                     << fine.Nm << std::endl;
        solver.calcCoarse(coarse.Cheby, par().smoother, par().coarseRelaxTol,
			              coarse.Nstop, coarse.Nk, coarse.Nm, coarse.resid, 
                          coarse.MaxIt, coarse.betastp, coarse.MinRes);
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSolver_LocalCoherenceLanczos_hpp_
