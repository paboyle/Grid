#ifndef Hadrons_MContraction_A2AMeson_hpp_
#define Hadrons_MContraction_A2AMeson_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>
#include <Grid/Hadrons/AllToAllVectors.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         A2AMeson                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class MesonFieldPar : Serializable
{
  public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(MesonFieldPar,
                                    int, Nl,
                                    int, N,
                                    std::string, A2A1,
                                    std::string, A2A2,
                                    std::string, action,
                                    std::string, epack1,
                                    std::string, epack2,
                                    std::string, gamma,
                                    std::string, output);
};

template <typename FImpl>
class TMesonFieldGmu : public Module<MesonFieldPar>
{
  public:
    FERM_TYPE_ALIASES(FImpl, );

    typedef A2AModesSchurDiagTwo<typename FImpl::FermionField, FMat> A2ABase;

    class Result : Serializable
    {
      public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
                                        Gamma::Algebra, gamma_mu,
                                        std::vector<std::vector<std::vector<ComplexD>>>, MesonField);
    };

  public:
    // constructor
    TMesonFieldGmu(const std::string name);
    // destructor
    virtual ~TMesonFieldGmu(void){};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER(A2AMeson, ARG(TMesonFieldGmu<FIMPL>), MContraction);

/******************************************************************************
*                  TMesonFieldGmu implementation                             *
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TMesonFieldGmu<FImpl>::TMesonFieldGmu(const std::string name)
    : Module<MesonFieldPar>(name)
{
}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TMesonFieldGmu<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().A2A1, par().A2A2, par().action};
    in.push_back(par().A2A1 + "_ret");
    in.push_back(par().A2A2 + "_ret");
    int Nl = par().Nl;
    if (Nl > 0)
    {
        in.push_back(par().epack1);
        in.push_back(par().epack2);
    }

    return in;
}

template <typename FImpl>
std::vector<std::string> TMesonFieldGmu<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TMesonFieldGmu<FImpl>::setup(void)
{
    int nt = env().getDim(Tp);
    int N = par().N;

    int Ls_ = env().getObjectLs(par().A2A1 + "_ret");

    envTmpLat(FermionField, "w", Ls_);
    envTmpLat(FermionField, "v", Ls_);
    envTmpLat(FermionField, "tmpv_5d", Ls_);
    envTmpLat(FermionField, "tmpw_5d", Ls_);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TMesonFieldGmu<FImpl>::execute(void)
{
    LOG(Message) << "Computing A2A meson field for gamma_mu  = " << par().gamma << ", taking w from " << par().A2A1 << " and v from " << par().A2A2 << std::endl;

    Result result;
    std::istringstream sstr(par().gamma);
    Gamma::Algebra g_mu;
    sstr >> g_mu;
    Gamma gamma_mu(g_mu);

    int N = par().N;
    LOG(Message) << "N for A2A cont: " << N << std::endl;
    result.gamma_mu = g_mu;
    result.MesonField.resize(N);

    int nt = env().getDim(Tp);
    std::vector<ComplexD> MesonField_ij;
    MesonField_ij.resize(nt);

    auto &a2a1 = envGet(A2ABase, par().A2A1 + "_ret");
    auto &a2a2 = envGet(A2ABase, par().A2A2 + "_ret");

    envGetTmp(FermionField, w);
    envGetTmp(FermionField, v);
    envGetTmp(FermionField, tmpv_5d);
    envGetTmp(FermionField, tmpw_5d);

    for (unsigned int i = 0; i < N; i++)
    {
        a2a1.return_w(i, tmpw_5d, w);
        for (unsigned int j = 0; j < N; j++)
        {
            a2a2.return_v(j, tmpv_5d, v);
            v = gamma_mu*v;
            sliceInnerProductVector(MesonField_ij, w, v, Tp);
            result.MesonField[j][i] = MesonField_ij;
        }
        if (i % 10 == 0)
        {
            LOG(Message) << "MF for i = " << i << " of " << N << std::endl;
        }
    }

    saveResult(par().output, "meson", result);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_A2AMeson_hpp_
