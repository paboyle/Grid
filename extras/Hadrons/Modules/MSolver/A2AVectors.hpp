#ifndef Hadrons_MSolver_A2AVectors_hpp_
#define Hadrons_MSolver_A2AVectors_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>
#include <Grid/Hadrons/Solver.hpp>
#include <Grid/Hadrons/EigenPack.hpp>
#include <Grid/Hadrons/AllToAllVectors.hpp>
#include <Grid/Hadrons/DilutedNoise.hpp>

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
                                  std::string, noise,
                                  std::string, action,
                                  std::string, eigenPack,
                                  std::string, solver);
};

template <typename FImpl, typename Pack>
class TA2AVectors : public Module<A2AVectorsPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    SOLVER_TYPE_ALIASES(FImpl,);
    typedef A2AModesSchurDiagTwo<FermionField, FMat, Solver> A2ABase;
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
};

MODULE_REGISTER_TMP(A2AVectors, 
    ARG(TA2AVectors<FIMPL, FermionEigenPack<FIMPL>>), MSolver);
MODULE_REGISTER_TMP(ZA2AVectors, 
    ARG(TA2AVectors<ZFIMPL, FermionEigenPack<ZFIMPL>>), MSolver);

/******************************************************************************
 *                 TA2AVectors implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
TA2AVectors<FImpl, Pack>::TA2AVectors(const std::string name)
: Module<A2AVectorsPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
std::vector<std::string> TA2AVectors<FImpl, Pack>::getInput(void)
{
    int Nl = par().Nl;
    std::string sub_string = "";
    if (Nl > 0) sub_string = "_subtract";

    std::vector<std::string> in = {par().solver + sub_string, par().noise};

    return in;
}

template <typename FImpl, typename Pack>
std::vector<std::string> TA2AVectors<FImpl, Pack>::getReference(void)
{
    std::vector<std::string> ref = {par().action};

    if (!par().eigenPack.empty())
    {
        ref.push_back(par().eigenPack);
    }

    return ref;
}

template <typename FImpl, typename Pack>
std::vector<std::string> TA2AVectors<FImpl, Pack>::getOutput(void)
{
    std::vector<std::string> out = {getName()};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
void TA2AVectors<FImpl, Pack>::setup(void)
{
    int Nl = par().Nl;
    bool return_5d = par().return_5d;   
    auto &noise  = envGet(DilutedNoise<FImpl>, par().noise);
    int Ls;

    std::string sub_string = "";
    if (Nl > 0) sub_string = "_subtract";
    auto &solver = envGet(Solver, par().solver + sub_string);
    Ls = env().getObjectLs(par().solver + sub_string);

    auto &action = envGet(FMat, par().action);

    envTmpLat(FermionField, "ferm_src", Ls);
    envTmpLat(FermionField, "unphys_ferm", Ls);
    envTmpLat(FermionField, "tmp");

    std::vector<FermionField> *evec;
    const std::vector<RealD> *eval;

    if (Nl > 0)
    {
        // Low modes
        auto &epack = envGet(Pack, par().eigenPack);

        LOG(Message) << "Creating a2a vectors " << getName() <<
                     " using eigenpack '" << par().eigenPack << "' ("
                     << epack.evec.size() << " modes)" <<
                     " and " << noise.size() << " high modes." << std::endl;
        evec = &epack.evec;
        eval = &epack.eval;
    }
    else
    {
        LOG(Message) << "Creating a2a vectors " << getName() <<
                     " using " << noise.size() << " high modes only." << std::endl;
    }

    envCreate(A2ABase, getName(), Ls, evec, eval, action, solver, Nl, noise.size(),
              return_5d);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
void TA2AVectors<FImpl, Pack>::execute(void)
{
    auto &action = envGet(FMat, par().action);
    auto &noise  = envGet(DilutedNoise<FImpl>, par().noise);

    int Ls;
    int Nl = par().Nl;

    std::string sub_string = "";
    if (Nl > 0) sub_string = "_subtract";
    Ls = env().getObjectLs(par().solver + sub_string);

    auto &a2areturn = envGet(A2ABase, getName());

    // High modes
    envGetTmp(FermionField, ferm_src);
    envGetTmp(FermionField, unphys_ferm);
    envGetTmp(FermionField, tmp);
    for (unsigned int i = 0; i < noise.size(); i++)
    {
        LOG(Message) << "A2A src for noise vector " << i << std::endl;
        // source conversion for 4D sources
        if (!env().isObject5d(par().noise))
        {
            if (Ls == 1)
            {
                ferm_src = noise[i];
                tmp = ferm_src;
            }
            else
            {
                tmp = noise[i];
                action.ImportPhysicalFermionSource(noise[i], ferm_src);
                action.ImportUnphysicalFermion(noise[i], unphys_ferm);
            }
        }
        // source conversion for 5D sources
        else
        {
            if (Ls != env().getObjectLs(par().noise))
            {
                HADRONS_ERROR(Size, "Ls mismatch between quark action and source");
            }
            else
            {
                ferm_src = noise[i];
                action.ExportPhysicalFermionSolution(ferm_src, tmp);
                unphys_ferm = ferm_src;
            }
        }
        LOG(Message) << "solveHighMode i = " << i << std::endl;
        a2areturn.high_modes(ferm_src, unphys_ferm, tmp, i);
    }
}
END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSolver_A2AVectors_hpp_
