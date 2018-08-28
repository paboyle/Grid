#ifndef Hadrons_MSolver_A2AVectors_hpp_
#define Hadrons_MSolver_A2AVectors_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>
#include <Grid/Hadrons/Solver.hpp>
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
                                  std::string, eigenPack,
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

    typedef A2AModesSchurDiagTwo<typename FImpl::FermionField, FMat, Solver> A2ABase;

  public:
    // constructor
    TA2AVectors(const std::string name);
    // destructor
    virtual ~TA2AVectors(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getReference(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);

  private:
    unsigned int Ls_;
    std::string className_;
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
, className_ (name + "_class")
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, int nBasis>
std::vector<std::string> TA2AVectors<FImpl, nBasis>::getInput(void)
{
    int Nl = par().Nl;
    std::string sub_string = "";
    if (Nl > 0) sub_string = "_subtract";
    std::vector<std::string> in = {par().solver + sub_string};

    int n = par().sources.size();

    for (unsigned int t = 0; t < n; t += 1)
    {
        in.push_back(par().sources[t]);
    }

    return in;
}

template <typename FImpl, int nBasis>
std::vector<std::string> TA2AVectors<FImpl, nBasis>::getReference(void)
{
    std::vector<std::string> ref = {par().action};

    if (!par().eigenPack.empty())
    {
        ref.push_back(par().eigenPack);
    }

    return ref;
}

template <typename FImpl, int nBasis>
std::vector<std::string> TA2AVectors<FImpl, nBasis>::getOutput(void)
{
    std::vector<std::string> out = {getName(), className_,
                                    getName() + "_w_high_5d", getName() + "_v_high_5d",
                                    getName() + "_w_high_4d", getName() + "_v_high_4d"};

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
    int Ls;

    std::string sub_string = "";
    if (Nl > 0) sub_string = "_subtract";
    auto &solver = envGet(Solver, par().solver + sub_string);
    Ls = env().getObjectLs(par().solver + sub_string);

    auto &action = envGet(FMat, par().action);

    envTmpLat(FermionField, "ferm_src", Ls);
    envTmpLat(FermionField, "unphys_ferm", Ls);
    envTmpLat(FermionField, "tmp");
    envTmpLat(FermionField, "tmp2");

    std::vector<FermionField> *evec;
    const std::vector<RealD> *eval;

    if (Nl > 0)
    {
        // Low modes
        auto &epack = envGet(EPack, par().eigenPack);

        LOG(Message) << "Creating a2a vectors " << getName() <<
                     " using eigenpack '" << par().eigenPack << "' ("
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

    int size_5d = 1;
    if (return_5d) size_5d = Nh;
    envCreate(std::vector<FermionField>, getName() + "_w_high_5d", Ls, size_5d, FermionField(env().getGrid(Ls)));
    envCreate(std::vector<FermionField>, getName() + "_v_high_5d", Ls, size_5d, FermionField(env().getGrid(Ls)));
    envCreate(std::vector<FermionField>, getName() + "_w_high_4d", 1, Nh, FermionField(env().getGrid(1)));
    envCreate(std::vector<FermionField>, getName() + "_v_high_4d", 1, Nh, FermionField(env().getGrid(1)));

    auto &w_high_5d = envGet(std::vector<FermionField>, getName() + "_w_high_5d");
    auto &v_high_5d = envGet(std::vector<FermionField>, getName() + "_v_high_5d");
    auto &w_high_4d = envGet(std::vector<FermionField>, getName() + "_w_high_4d");
    auto &v_high_4d = envGet(std::vector<FermionField>, getName() + "_v_high_4d");

    envCreate(A2ABase, className_, Ls,
                     evec, eval,
                     action,
                     solver,
                     w_high_5d, v_high_5d,
                     w_high_4d, v_high_4d,
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

    auto &a2areturn = envGet(A2ABase, className_);

    // High modes
    auto sources = par().sources;
    int Nsrc = par().sources.size();

    envGetTmp(FermionField, ferm_src);
    envGetTmp(FermionField, unphys_ferm);
    envGetTmp(FermionField, tmp);
    envGetTmp(FermionField, tmp2);

    int N_count = 0;
    for (unsigned int T = 0; T < Nsrc; T++)
    for (unsigned int s = 0; s < Ns; ++s)
    for (unsigned int c = 0; c < Nc; ++c)
        {
            auto &prop_src = envGet(PropagatorField, sources[T]);
            LOG(Message) << "A2A src for s = " << s << " , c = " << c << ", T = " << T << std::endl;
            // source conversion for 4D sources
            if (!env().isObject5d(sources[T]))
            {
                if (Ls_ == 1)
                {
                    PropToFerm<FImpl>(ferm_src, prop_src, s, c);
                    tmp = ferm_src;
                }
                else
                {
                    PropToFerm<FImpl>(tmp, prop_src, s, c);
                    action.ImportPhysicalFermionSource(tmp, ferm_src);
                    action.ImportUnphysicalFermion(tmp, unphys_ferm);
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
                    action.ExportPhysicalFermionSolution(ferm_src, tmp);
                    unphys_ferm = ferm_src;
                }
            }
            LOG(Message) << "a2areturn.high_modes Ncount = " << N_count << std::endl;
            a2areturn.high_modes(ferm_src, unphys_ferm, tmp, N_count);
            N_count++;
        }
}
END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSolver_A2AVectors_hpp_
