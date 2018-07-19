#ifndef Hadrons_MContraction_MesonFieldGamma_hpp_
#define Hadrons_MContraction_MesonFieldGamma_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>
#include <Grid/Hadrons/AllToAllVectors.hpp>
#include <Grid/Hadrons/AllToAllReduction.hpp>
#include <Grid/Grid_Eigen_Dense.h>
#include <fstream>

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
                                    int, Nblock,
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
    virtual void vectorOfWs(std::vector<FermionField> &w, int i, int Nblock, FermionField &tmpw_5d, std::vector<FermionField> &vec_w);
    virtual void vectorOfVs(std::vector<FermionField> &v, int j, int Nblock, FermionField &tmpv_5d, std::vector<FermionField> &vec_v);
    virtual void gammaMult(std::vector<FermionField> &v, Gamma gamma);
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

template <typename FImpl>
void TMesonFieldGamma<FImpl>::vectorOfWs(std::vector<FermionField> &w, int i, int Nblock, FermionField &tmpw_5d, std::vector<FermionField> &vec_w)
{
    for (unsigned int ni = 0; ni < Nblock; ni++)
    {
        vec_w[ni] = w[i + ni];
    }
}

template <typename FImpl>
void TMesonFieldGamma<FImpl>::vectorOfVs(std::vector<FermionField> &v, int j, int Nblock, FermionField &tmpv_5d, std::vector<FermionField> &vec_v)
{
    for (unsigned int nj = 0; nj < Nblock; nj++)
    {
        vec_v[nj] = v[j+nj];
    }
}

template <typename FImpl>
void TMesonFieldGamma<FImpl>::gammaMult(std::vector<FermionField> &v, Gamma gamma)
{
    int Nblock = v.size();
    for (unsigned int nj = 0; nj < Nblock; nj++)
    {
        v[nj] = gamma * v[nj];
    }
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TMesonFieldGamma<FImpl>::setup(void)
{
    int nt = env().getDim(Tp);
    int N = par().N;
    int Nblock = par().Nblock;

    int Ls_ = env().getObjectLs(par().A2A1 + "_class");

    envTmpLat(FermionField, "tmpv_5d", Ls_);
    envTmpLat(FermionField, "tmpw_5d", Ls_);

    envTmp(std::vector<FermionField>, "w", 1, N, FermionField(env().getGrid(1)));
    envTmp(std::vector<FermionField>, "v", 1, N, FermionField(env().getGrid(1)));

    envTmp(Eigen::MatrixXcd, "MF", 1, Eigen::MatrixXcd::Zero(nt, N * N));

    envTmp(std::vector<FermionField>, "w_block", 1, Nblock, FermionField(env().getGrid(1)));
    envTmp(std::vector<FermionField>, "v_block", 1, Nblock, FermionField(env().getGrid(1)));
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TMesonFieldGamma<FImpl>::execute(void)
{
    LOG(Message) << "Computing A2A meson field for gamma = " << par().gammas << ", taking w from " << par().A2A1 << " and v from " << par().A2A2 << std::endl;

    int N = par().N;
    int nt = env().getDim(Tp);
    int Nblock = par().Nblock;

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

    auto &a2a1 = envGet(A2ABase, par().A2A1 + "_class");
    auto &a2a2 = envGet(A2ABase, par().A2A2 + "_class");

    envGetTmp(FermionField, tmpv_5d);
    envGetTmp(FermionField, tmpw_5d);

    envGetTmp(std::vector<FermionField>, v);
    envGetTmp(std::vector<FermionField>, w);
    LOG(Message) << "Finding v and w vectors for N =  " << N << std::endl;
    for (int i = 0; i < N; i++)
    {
        a2a2.return_v(i, tmpv_5d, v[i]);
        a2a1.return_w(i, tmpw_5d, w[i]);
    }
    LOG(Message) << "Found v and w vectors for N =  " << N << std::endl;

    std::vector<std::vector<ComplexD>> MesonField_ij;
    LOG(Message) << "Before blocked MFs, Nblock = " << Nblock << std::endl;
    envGetTmp(std::vector<FermionField>, v_block);
    envGetTmp(std::vector<FermionField>, w_block);
    MesonField_ij.resize(Nblock * Nblock, std::vector<ComplexD>(nt));

    envGetTmp(Eigen::MatrixXcd, MF);

    LOG(Message) << "Before blocked MFs, Nblock = " << Nblock << std::endl;
    for (unsigned int i = 0; i < N; i += Nblock)
    {
        vectorOfWs(w, i, Nblock, tmpw_5d, w_block);
        for (unsigned int j = 0; j < N; j += Nblock)
        {
            vectorOfVs(v, j, Nblock, tmpv_5d, v_block);
            for (unsigned int k = 0; k < result.size(); k++)
            {
                gammaMult(v_block, gammaList[k]);
                sliceInnerProductMesonField(MesonField_ij, w_block, v_block, Tp);
                for (unsigned int nj = 0; nj < Nblock; nj++)
                {
                    for (unsigned int ni = 0; ni < Nblock; ni++)
                    {
                        MF.col((i + ni) + (j + nj) * N) = Eigen::VectorXcd::Map(&MesonField_ij[nj * Nblock + ni][0], MesonField_ij[nj * Nblock + ni].size());
                    }
                }
            }
        }
        if (i % 10 == 0)
        {
            LOG(Message) << "MF for i = " << i << " of " << N << std::endl;
        }
    }
    LOG(Message) << "Before Global sum, Nblock = " << Nblock << std::endl;
    v_block[0]._grid->GlobalSumVector(MF.data(), MF.size());
    LOG(Message) << "After Global sum, Nblock = " << Nblock << std::endl;
    for (unsigned int i = 0; i < N; i++)
    {
        for (unsigned int j = 0; j < N; j++)
        {
            for (unsigned int k = 0; k < result.size(); k++)
            {
                for (unsigned int t = 0; t < nt; t++)
                {
                    result[k].MesonField[i][j][t] = MF.col(i + N * j)[t];
                }
            }
        }
    }
    saveResult(par().output, "meson", result);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_MesonFieldGm_hpp_
