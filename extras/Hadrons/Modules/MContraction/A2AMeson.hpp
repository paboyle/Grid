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

typedef std::pair<Gamma::Algebra, Gamma::Algebra> GammaPair;

class A2AMesonPar : Serializable
{
  public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2AMesonPar,
                                    int, Nl,
                                    int, N,
                                    std::string, A2A1,
                                    std::string, A2A2,
                                    std::string, gammas,
                                    std::string, output);
};

template <typename FImpl>
class TA2AMeson : public Module<A2AMesonPar>
{
  public:
    FERM_TYPE_ALIASES(FImpl, );
    SOLVER_TYPE_ALIASES(FImpl, );

    typedef A2AModesSchurDiagTwo<typename FImpl::FermionField, FMat, Solver> A2ABase;

    class Result : Serializable
    {
      public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
                                        Gamma::Algebra, gamma_snk,
                                        Gamma::Algebra, gamma_src,
                                        std::vector<Complex>, corr);
    };

  public:
    // constructor
    TA2AMeson(const std::string name);
    // destructor
    virtual ~TA2AMeson(void){};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual void parseGammaString(std::vector<GammaPair> &gammaList);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER(A2AMeson, ARG(TA2AMeson<FIMPL>), MContraction);
MODULE_REGISTER(ZA2AMeson, ARG(TA2AMeson<ZFIMPL>), MContraction);

/******************************************************************************
*                  TA2AMeson implementation                             *
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TA2AMeson<FImpl>::TA2AMeson(const std::string name)
    : Module<A2AMesonPar>(name)
{
}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TA2AMeson<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().A2A1, par().A2A2};
    in.push_back(par().A2A1 + "_class");
    in.push_back(par().A2A2 + "_class");

    return in;
}

template <typename FImpl>
std::vector<std::string> TA2AMeson<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {};

    return out;
}

template <typename FImpl>
void TA2AMeson<FImpl>::parseGammaString(std::vector<GammaPair> &gammaList)
{
    gammaList.clear();
    // Parse individual contractions from input string.
    gammaList = strToVec<GammaPair>(par().gammas);
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2AMeson<FImpl>::setup(void)
{
    int nt = env().getDim(Tp);
    int N = par().N;

    int Ls_ = env().getObjectLs(par().A2A1 + "_class");

    envTmp(std::vector<FermionField>, "w1", 1, N, FermionField(env().getGrid(1)));
    envTmp(std::vector<FermionField>, "v1", 1, N, FermionField(env().getGrid(1)));
    envTmpLat(FermionField, "tmpv_5d", Ls_);
    envTmpLat(FermionField, "tmpw_5d", Ls_);

    envTmp(std::vector<ComplexD>, "MF_x", 1, nt);
    envTmp(std::vector<ComplexD>, "MF_y", 1, nt);
    envTmp(std::vector<ComplexD>, "tmp", 1, nt);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2AMeson<FImpl>::execute(void)
{
    LOG(Message) << "Computing A2A meson contractions" << std::endl;

    Result result;
    Gamma g5(Gamma::Algebra::Gamma5);
    std::vector<GammaPair> gammaList;
    int nt = env().getDim(Tp);

    parseGammaString(gammaList);

    result.gamma_snk = gammaList[0].first;
    result.gamma_src = gammaList[0].second;
    result.corr.resize(nt);

    int Nl = par().Nl;
    int N = par().N;
    LOG(Message) << "N for A2A cont: " << N << std::endl;

    envGetTmp(std::vector<ComplexD>, MF_x);
    envGetTmp(std::vector<ComplexD>, MF_y);
    envGetTmp(std::vector<ComplexD>, tmp);

    for (unsigned int t = 0; t < nt; ++t)
    {
        tmp[t] = TensorRemove(MF_x[t] * MF_y[t] * 0.0);
    }

    Gamma gSnk(gammaList[0].first);
    Gamma gSrc(gammaList[0].second);

    auto &a2a1_fn = envGet(A2ABase, par().A2A1 + "_class");

    envGetTmp(std::vector<FermionField>, w1);
    envGetTmp(std::vector<FermionField>, v1);
    envGetTmp(FermionField, tmpv_5d);
    envGetTmp(FermionField, tmpw_5d);

    LOG(Message) << "Finding v and w vectors for N =  " << N << std::endl;
    for (int i = 0; i < N; i++)
    {
        a2a1_fn.return_v(i, tmpv_5d, v1[i]);
        a2a1_fn.return_w(i, tmpw_5d, w1[i]);
    }
    LOG(Message) << "Found v and w vectors for N =  " << N << std::endl;
    for (unsigned int i = 0; i < N; i++)
    {
        v1[i] = gSnk * v1[i];
    }
    int ty;
    for (unsigned int i = 0; i < N; i++)
    {
        for (unsigned int j = 0; j < N; j++)
        {
            mySliceInnerProductVector(MF_x, w1[i], v1[j], Tp);
            mySliceInnerProductVector(MF_y, w1[j], v1[i], Tp);
            for (unsigned int t = 0; t < nt; ++t)
            {
                for (unsigned int tx = 0; tx < nt; tx++)
                {
                    ty = (tx + t) % nt;
                    tmp[t] += TensorRemove((MF_x[tx]) * (MF_y[ty]));
                }
            }
        }
        if (i % 10 == 0)
        {
            LOG(Message) << "MF for i = " << i << " of " << N << std::endl;
        }
    }
    double NTinv = 1.0 / static_cast<double>(nt);
    for (unsigned int t = 0; t < nt; ++t)
    {
        result.corr[t] = NTinv * tmp[t];
    }

    saveResult(par().output, "meson", result);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_A2AMeson_hpp_
