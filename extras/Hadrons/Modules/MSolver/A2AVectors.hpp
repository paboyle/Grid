#ifndef Hadrons_MSolver_A2AVectors_hpp_
#define Hadrons_MSolver_A2AVectors_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>
#include <Grid/Hadrons/Solver.hpp>
#include <Grid/Hadrons/EigenPack.hpp>
#include <Grid/Hadrons/A2AVectors.hpp>
#include <Grid/Hadrons/DilutedNoise.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                     Create all-to-all vector class                         *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSolver)

class A2AVectorsPar: Serializable
{
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(A2AVectorsPar,
                                  bool, return_5d,
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
    typedef A2AVectorsSchurDiagTwo<FImpl> A2A;
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
    std::string  solverName_;
    unsigned int Nl_{0};
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
    std::string              sub_string;
    std::vector<std::string> in;

    if (!par().eigenPack.empty())
    {
        in.push_back(par().eigenPack);
        sub_string = (!par().eigenPack.empty()) ? "_subtract" : "";
    }
    in.push_back(par().solver + sub_string);
    in.push_back(par().noise);

    return in;
}

template <typename FImpl, typename Pack>
std::vector<std::string> TA2AVectors<FImpl, Pack>::getOutput(void)
{
    std::vector<std::string> out = {getName() + "_v", getName() + "_w"};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
void TA2AVectors<FImpl, Pack>::setup(void)
{
    bool        hasLowModes = (!par().eigenPack.empty());
    std::string sub_string  = (hasLowModes) ? "_subtract" : "";
    bool        return_5d   = par().return_5d;   
    auto        &noise      = envGet(DilutedNoise<FImpl>, par().noise);
    auto        &action     = envGet(FMat, par().action);
    auto        &solver     = envGet(Solver, par().solver + sub_string);
    int         Ls          = env().getObjectLs(par().action);

    LOG(Message) << "Creating all-to-all vectors ";
    if (hasLowModes)
    {
        auto &epack = envGet(Pack, par().eigenPack);
        
        Nl_ = epack.evec.size();
        std::cout << " using eigenpack '" << par().eigenPack << "' ("
                  << Nl_ << " low modes) and noise '"
                  << par().noise << "' (" << noise.size() 
                  << " noise vectors)" << std::endl;
    }
    else
    {
        std::cout << " using noise '" << par().noise << "' (" << noise.size() 
                  << " noise vectors)" << std::endl;
    }
    envCreate(std::vector<FermionField>, getName() + "_v", 1, 
              Nl_ + noise.size(), FermionField(env().getGrid()));
    envCreate(std::vector<FermionField>, getName() + "_w", 1, 
              Nl_ + noise.size(), FermionField(env().getGrid()));
    if (Ls > 1)
    {
        envTmpLat(FermionField, "f5", Ls);
    }
    envTmp(A2A, "a2a", 1, action, solver);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
void TA2AVectors<FImpl, Pack>::execute(void)
{
    std::string sub_string = (Nl_ > 0) ? "_subtract" : "";
    auto        &action    = envGet(FMat, par().action);
    auto        &solver    = envGet(Solver, par().solver + sub_string);
    auto        &noise     = envGet(DilutedNoise<FImpl>, par().noise);
    auto        &v         = envGet(std::vector<FermionField>, getName() + "_v");
    auto        &w         = envGet(std::vector<FermionField>, getName() + "_w");
    int         Ls         = env().getObjectLs(par().action);

    envGetTmp(A2A, a2a);
    // Low modes
    for (unsigned int il = 0; il < Nl_; il++)
    {
        auto &epack  = envGet(Pack, par().eigenPack);

        LOG(Message) << "V vector i = " << il << " (low mode)" << std::endl;
        if (Ls == 1)
        {
            a2a.makeLowModeV(v[il], epack.evec[il], epack.eval[il]);
        }
        else
        {
            envGetTmp(FermionField, f5);
            a2a.makeLowModeV5D(v[il], f5, epack.evec[il], epack.eval[il]);
        }
        LOG(Message) << "W vector i = " << il << " (low mode)" << std::endl;
        if (Ls == 1)
        {
            a2a.makeLowModeW(w[il], epack.evec[il], epack.eval[il]);
        }
        else
        {
            envGetTmp(FermionField, f5);
            a2a.makeLowModeW5D(w[il], f5, epack.evec[il], epack.eval[il]);
        }   
    }

    // High modes
    for (unsigned int ih = 0; ih < noise.size(); ih++)
    {
        LOG(Message) << "V vector i = " << Nl_ + ih
                     << " (" << ((Nl_ > 0) ? "high " : "") 
                     << "stochastic mode)" << std::endl;
        if (Ls == 1)
        {
            a2a.makeHighModeV(v[Nl_ + ih], noise[ih]);
        }
        else
        {
            envGetTmp(FermionField, f5);
            a2a.makeHighModeV5D(v[Nl_ + ih], f5, noise[ih]);
            std::cout << norm2(v[Nl_ + ih]) << std::endl;
        }
        LOG(Message) << "W vector i = " << Nl_ + ih
                     << " (" << ((Nl_ > 0) ? "high " : "") 
                     << "stochastic mode)" << std::endl;
        if (Ls == 1)
        {
            a2a.makeHighModeW(w[Nl_ + ih], noise[ih]);
        }
        else
        {
            envGetTmp(FermionField, f5);
            a2a.makeHighModeW5D(w[Nl_ + ih], f5, noise[ih]);
            std::cout << norm2(w[Nl_ + ih]) << std::endl;
        }
    }


        // // source conversion for 4D sources
        // if (!env().isObject5d(par().noise))
        // {
        //     if (Ls == 1)
        //     {
        //         ferm_src = noise[ih];
        //         tmp = ferm_src;
        //     }
        //     else
        //     {
        //         tmp = noise[ih];
        //         action.ImportPhysicalFermionSource(noise[ih], ferm_src);
        //         action.ImportUnphysicalFermion(noise[ih], unphys_ferm);
        //     }
        // }
        // // source conversion for 5D sources
        // else
        // {
        //     if (Ls != env().getObjectLs(par().noise))
        //     {
        //         HADRONS_ERROR(Size, "Ls mismatch between quark action and source");
        //     }
        //     else
        //     {
        //         ferm_src = noise[ih];
        //         action.ExportPhysicalFermionSolution(ferm_src, tmp);
        //         unphys_ferm = ferm_src;
        //     }
        // }
        // a2a.high_modes(ih, ferm_src, unphys_ferm, tmp, solver);
}
END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSolver_A2AVectors_hpp_
