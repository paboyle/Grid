#ifndef Hadrons_MSolver_A2AVectors_hpp_
#define Hadrons_MSolver_A2AVectors_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>
#include <Grid/Hadrons/EigenPack.hpp>
#include <Grid/Hadrons/AllToAllVectors.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         A2AVectors                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSolver)

class A2AVectorsPar: Serializable
{
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(A2AVectorsPar,
                                  bool, return_5d,
                                  int, Nl,
                                  int, N,
                                  std::vector<std::string>, sources,
                                  std::string, action,
                                  std::string, eigenpack,
                                  std::string, solver);
};

template <typename FImpl, int nBasis>
class TA2AVectors : public Module<A2AVectorsPar>
{
  public:
    FERM_TYPE_ALIASES(FImpl,);
    SOLVER_TYPE_ALIASES(FImpl,);

    typedef FermionEigenPack<FImpl> EPack;
    typedef CoarseFermionEigenPack<FImpl, nBasis> CoarseEPack;

    typedef A2AModesSchurDiagTwo<typename FImpl::FermionField, FMat> A2ABase;

  public:
    // constructor
    TA2AVectors(const std::string name);
    // destructor
    virtual ~TA2AVectors(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);

  private:
    unsigned int Ls_;
    std::string retName_;
};

MODULE_REGISTER_TMP(A2AVectors, ARG(TA2AVectors<FIMPL, HADRONS_DEFAULT_LANCZOS_NBASIS>), MSolver);
MODULE_REGISTER_TMP(ZA2AVectors, ARG(TA2AVectors<ZFIMPL, HADRONS_DEFAULT_LANCZOS_NBASIS>), MSolver);

/******************************************************************************
 *                 TA2AVectors implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, int nBasis>
TA2AVectors<FImpl, nBasis>::TA2AVectors(const std::string name)
: Module<A2AVectorsPar>(name)
, retName_ (name + "_ret")
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, int nBasis>
std::vector<std::string> TA2AVectors<FImpl, nBasis>::getInput(void)
{
    std::vector<std::string> in = {par().action, par().solver, par().solver + "_subtract"};

    int n = par().sources.size();

    for (unsigned int t = 0; t < n; t += 1)
    {
        in.push_back(par().sources[t]);
    }

    return in;
}

template <typename FImpl, int nBasis>
std::vector<std::string> TA2AVectors<FImpl, nBasis>::getOutput(void)
{
    std::vector<std::string> out = {getName(), retName_};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, int nBasis>
void TA2AVectors<FImpl, nBasis>::setup(void)
{
    int N = par().N;
    int Nl = par().Nl;
    int Nh = N - Nl;
    bool return_5d = par().return_5d;
    int Ls, Ls_;

    std::string sub_string = "";
    if (Nl > 0) sub_string = "_subtract";
    auto &solver = envGet(SolverFn, par().solver + sub_string);
    Ls_ = env().getObjectLs(par().solver + sub_string);

    auto &action = envGet(FMat, par().action);

    envTmpLat(FermionField, "ferm_src", Ls_);
    envTmpLat(FermionField, "tmp");
    envTmpLat(FermionField, "tmp2");

    std::vector<FermionField> *evec;
    const std::vector<RealD> *eval;

    if (Nl > 0)
    {
        // Low modes
        auto &epack = envGet(EPack, par().eigenpack);

        LOG(Message) << "Creating a2a vectors " << getName() <<
                     " using eigenpack '" << par().eigenpack << "' ("
                     << epack.evec.size() << " modes)" <<
                     " and " << Nh << " high modes." << std::endl;
        evec = &epack.evec;
        eval = &epack.eval;
    }
    else
    {
        LOG(Message) << "Creating a2a vectors " << getName() <<
                     " using " << Nh << " high modes only." << std::endl;
    }
    envCreate(A2ABase, retName_, Ls_,
                     evec, eval,
                     action,
                     solver,
                     Nl, Nh,
                     return_5d);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, int nBasis>
void TA2AVectors<FImpl, nBasis>::execute(void)
{
    auto &action = envGet(FMat, par().action);

    int Nt = env().getDim(Tp);
    int Nc = FImpl::Dimension;
    int Ls_;
    int Nl = par().Nl;

    std::string sub_string = "";
    if (Nl > 0) sub_string = "_subtract";
    Ls_ = env().getObjectLs(par().solver + sub_string);

    auto &a2areturn = envGet(A2ABase, retName_);

    // High modes
    auto sources = par().sources;
    int Nsrc = par().sources.size();

    envGetTmp(FermionField, ferm_src);
    envGetTmp(FermionField, tmp);
    envGetTmp(FermionField, tmp2);

    // TODO: At the moment weighting only applies to the 4d->5d source path
    // similar to how the 5d and 4d srcs are passed in, this needs more work to be less brittle
    double weight = 1.0 / sqrt(Ns*Nc*Nsrc);
    int N_count = 0;
    for (unsigned int s = 0; s < Ns; ++s)
        for (unsigned int c = 0; c < Nc; ++c)
        for (unsigned int T = 0; T < Nsrc; T++)
        {
            auto &prop_src = envGet(PropagatorField, sources[T]);
            LOG(Message) << "A2A src for s = " << s << " , c = " << c << ", T = " << T << std::endl;
            // source conversion for 4D sources
            if (!env().isObject5d(sources[T]))
            {
                if (Ls_ == 1)
                {
                    PropToFerm<FImpl>(ferm_src, prop_src, s, c);
                }
                else
                {
                    PropToFerm<FImpl>(tmp2, prop_src, s, c);
                    tmp = weight*tmp2;
                    action.ImportPhysicalFermionSource(tmp, ferm_src);
                }
            }
            // source conversion for 5D sources
            else
            {
                if (Ls_ != env().getObjectLs(sources[T]))
                {
                    HADRONS_ERROR(Size, "Ls mismatch between quark action and source");
                }
                else
                {
                    PropToFerm<FImpl>(ferm_src, prop_src, s, c);
                }
            }
            LOG(Message) << "a2areturn.high_modes Ncount = " << N_count << std::endl;
            a2areturn.high_modes(ferm_src, tmp, N_count);
            N_count++;
        }
}
END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSolver_A2AVectors_hpp_
