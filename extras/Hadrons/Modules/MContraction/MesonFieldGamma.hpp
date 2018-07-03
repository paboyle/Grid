#ifndef Hadrons_MContraction_MesonFieldGamma_hpp_
#define Hadrons_MContraction_MesonFieldGamma_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>
#include <Grid/Hadrons/AllToAllVectors.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         MesonFieldGamma                                 *
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
                                    std::string, gammas,
                                    std::string, output);
};

template <typename FImpl>
class TMesonFieldGamma : public Module<MesonFieldPar>
{
  public:
    FERM_TYPE_ALIASES(FImpl, );
    SOLVER_TYPE_ALIASES(FImpl, );

    typedef A2AModesSchurDiagTwo<typename FImpl::FermionField, FMat, Solver> A2ABase;

    class Result : Serializable
    {
      public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
                                        Gamma::Algebra, gamma,
                                        std::vector<std::vector<std::vector<ComplexD>>>, MesonField);
    };

  public:
    // constructor
    TMesonFieldGamma(const std::string name);
    // destructor
    virtual ~TMesonFieldGamma(void){};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual void parseGammaString(std::vector<Gamma::Algebra> &gammaList);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER(MesonFieldGamma, ARG(TMesonFieldGamma<FIMPL>), MContraction);
MODULE_REGISTER(ZMesonFieldGamma, ARG(TMesonFieldGamma<ZFIMPL>), MContraction);

/******************************************************************************
*                  TMesonFieldGamma implementation                             *
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TMesonFieldGamma<FImpl>::TMesonFieldGamma(const std::string name)
    : Module<MesonFieldPar>(name)
{
}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TMesonFieldGamma<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().A2A1 + "_class", par().A2A2 + "_class"};
    return in;
}

template <typename FImpl>
std::vector<std::string> TMesonFieldGamma<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {};

    return out;
}

template <typename FImpl>
void TMesonFieldGamma<FImpl>::parseGammaString(std::vector<Gamma::Algebra> &gammaList)
{
    gammaList.clear();
    // Determine gamma matrices to insert at source/sink.
    if (par().gammas.compare("all") == 0)
    {
        // Do all contractions.
        for (unsigned int i = 1; i < Gamma::nGamma; i += 2)
        {
            gammaList.push_back(((Gamma::Algebra)i));
        }
    }
    else
    {
        // Parse individual contractions from input string.
        gammaList = strToVec<Gamma::Algebra>(par().gammas);
    }
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TMesonFieldGamma<FImpl>::setup(void)
{
    int nt = env().getDim(Tp);
    int N = par().N;

    int Ls_ = env().getObjectLs(par().A2A1 + "_class");

    envTmpLat(FermionField, "w", Ls_);
    envTmpLat(FermionField, "v", Ls_);
    envTmpLat(FermionField, "tmpv_5d", Ls_);
    envTmpLat(FermionField, "tmpw_5d", Ls_);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TMesonFieldGamma<FImpl>::execute(void)
{
    LOG(Message) << "Computing A2A meson field for gamma = " << par().gammas << ", taking w from " << par().A2A1 << " and v from " << par().A2A2 << std::endl;

    int N = par().N;
    int nt = env().getDim(Tp);

    std::vector<Result> result;
    std::vector<Gamma::Algebra> gammaResultList;
    std::vector<Gamma> gammaList;

    parseGammaString(gammaResultList);
    result.resize(gammaResultList.size());

    Gamma g5(Gamma::Algebra::Gamma5);
    gammaList.resize(gammaResultList.size(), g5);

    for (unsigned int i = 0; i < result.size(); ++i)
    {
        result[i].gamma = gammaResultList[i];
        result[i].MesonField.resize(N, std::vector<std::vector<ComplexD>>(N, std::vector<ComplexD>(nt)));

        Gamma gamma(gammaResultList[i]);
        gammaList[i] = gamma;
    }

    std::vector<ComplexD> MesonField_ij;
    MesonField_ij.resize(nt);

    auto &a2a1 = envGet(A2ABase, par().A2A1 + "_class");
    auto &a2a2 = envGet(A2ABase, par().A2A2 + "_class");

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
            for (unsigned int k = 0; k < result.size(); k++)
            {
                v = gammaList[k]*v;
                sliceInnerProductVector(MesonField_ij, w, v, Tp);
                result[k].MesonField[i][j] = MesonField_ij;
            }
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

#endif // Hadrons_MContraction_MesonFieldGm_hpp_
