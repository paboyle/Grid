#ifndef Hadrons_MContraction_A2AWeakHamiltonian_hpp_
#define Hadrons_MContraction_A2AWeakHamiltonian_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/A2AVectors.hpp>
#include <Hadrons/A2AMatrix.hpp>
#include <Hadrons/DiskVector.hpp>

template <typename FImpl>
class A2A_DVutils
{
  public:
    typedef typename FImpl::ComplexField ComplexField;
    typedef typename FImpl::FermionField FermionField;
    typedef typename FImpl::PropagatorField PropagatorField;

    typedef typename FImpl::SiteSpinor vobj;
    typedef typename vobj::scalar_object sobj;
    typedef typename vobj::scalar_type scalar_type;
    typedef typename vobj::vector_type vector_type;

    typedef Grid::iSpinMatrix<vector_type> SpinMatrix_v;
    typedef Grid::iSpinMatrix<scalar_type> SpinMatrix_s;
    typedef Grid::iSinglet<vector_type> Scalar_v;
    typedef Grid::iSinglet<scalar_type> Scalar_s;

    typedef Grid::iSpinColourMatrix<vector_type> SpinColourMatrix_v;

    template <typename TensorType> // input: rank 3 tensor, e.g. Eigen::Tensor<ComplexD, 3>, EigenDiskVector<ComplexD>
    static void ContractWWVV_dv(std::vector<PropagatorField> &WWVV,
                                const TensorType &WW_sd,
                                const FermionField *vs,
                                const FermionField *vd);
};
template <class FImpl>
template <typename TensorType>
void A2A_DVutils<FImpl>::ContractWWVV_dv(std::vector<PropagatorField> &WWVV,
                                         const TensorType &WW_sd,
                                         const FermionField *vs,
                                         const FermionField *vd)
{
    Grid::GridBase *grid = vs[0]._grid;

    int nd = grid->_ndimension;
    int Nsimd = grid->Nsimd();
    int N_t = WW_sd.dimensions()[0];
    int N_s = WW_sd.dimensions()[1];
    int N_d = WW_sd.dimensions()[2];

    int d_unroll = 32; // Empirical optimisation

    for (int t = 0; t < N_t; t++)
    {
        std::cout << " Contraction t = " << t << std::endl;
        Grid::Hadrons::EigenDiskVector<Grid::ComplexD>::Matrix buf = WW_sd[t];
        parallel_for(int ss = 0; ss < grid->oSites(); ss++)
        {
            for (int d_o = 0; d_o < N_d; d_o += d_unroll)
            {
                for (int s = 0; s < N_s; s++)
                {
                    auto tmp1 = vs[s]._odata[ss];
                    vobj tmp2 = tmp1 * 0.0; // would rather write zero

                    for (int d = d_o; d < MIN(d_o + d_unroll, N_d); d++)
                    {
                        Scalar_v coeff = buf(s, d);
                        mac(&tmp2, &coeff, &vd[d]._odata[ss]);
                    }

                    //////////////////////////
                    // Fast outer product of tmp1 with a sum of terms suppressed by d_unroll
                    //////////////////////////
                    tmp2 = conjugate(tmp2);
                    for (int s1 = 0; s1 < Grid::Ns; s1++)
                    {
                        for (int s2 = 0; s2 < Grid::Ns; s2++)
                        {
                            WWVV[t]._odata[ss]()(s1, s2)(0, 0) += tmp1()(s1)(0) * tmp2()(s2)(0);
                            WWVV[t]._odata[ss]()(s1, s2)(0, 1) += tmp1()(s1)(0) * tmp2()(s2)(1);
                            WWVV[t]._odata[ss]()(s1, s2)(0, 2) += tmp1()(s1)(0) * tmp2()(s2)(2);
                            WWVV[t]._odata[ss]()(s1, s2)(1, 0) += tmp1()(s1)(1) * tmp2()(s2)(0);
                            WWVV[t]._odata[ss]()(s1, s2)(1, 1) += tmp1()(s1)(1) * tmp2()(s2)(1);
                            WWVV[t]._odata[ss]()(s1, s2)(1, 2) += tmp1()(s1)(1) * tmp2()(s2)(2);
                            WWVV[t]._odata[ss]()(s1, s2)(2, 0) += tmp1()(s1)(2) * tmp2()(s2)(0);
                            WWVV[t]._odata[ss]()(s1, s2)(2, 1) += tmp1()(s1)(2) * tmp2()(s2)(1);
                            WWVV[t]._odata[ss]()(s1, s2)(2, 2) += tmp1()(s1)(2) * tmp2()(s2)(2);
                        }
                    }
                }
            }
        }
    }
}

BEGIN_HADRONS_NAMESPACE
/******************************************************************************
 *                         A2AWeakHamiltonian                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class A2AWeakHamiltonianPar : Serializable
{
  public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2AWeakHamiltonianPar,
                                    std::string, vStrange,
                                    std::string, vLight1,
                                    std::string, wLight1,
                                    std::string, vLight2,
                                    std::string, vCharm,
                                    std::string, wCharm,
                                    std::string, mfDirLL,
                                    std::string, mfDirLS,
                                    std::string, output,
                                    int, dt_min,
                                    int, dt_max);
};

template <typename FImpl>
class TA2AWeakHamiltonian : public Module<A2AWeakHamiltonianPar>
{
  public:
    FERM_TYPE_ALIASES(FImpl, );
    SOLVER_TYPE_ALIASES(FImpl, );
    class Result : Serializable
    {
      public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
                                        int, dt_min,
                                        int, dt_max,
                                        std::vector<Complex>, corr_wing,
                                        std::vector<Complex>, corr_conncected,
                                        std::vector<Complex>, corr_eye_light,
                                        std::vector<Complex>, corr_saucer_light,
                                        std::vector<Complex>, corr_eye_charm,
                                        std::vector<Complex>, corr_saucer_charm);
    };

  public:
    // constructor
    TA2AWeakHamiltonian(const std::string name);
    // destructor
    virtual ~TA2AWeakHamiltonian(void){};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(A2AWeakHamiltonian, TA2AWeakHamiltonian<FIMPL>, MContraction);

/******************************************************************************
 *                 TA2AWeakHamiltonian implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TA2AWeakHamiltonian<FImpl>::TA2AWeakHamiltonian(const std::string name)
    : Module<A2AWeakHamiltonianPar>(name)
{
}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TA2AWeakHamiltonian<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().vStrange,
                                   par().vLight1, par().wLight1,
                                   par().vLight2,
                                   par().vCharm, par().wCharm};

    return in;
}

template <typename FImpl>
std::vector<std::string> TA2AWeakHamiltonian<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2AWeakHamiltonian<FImpl>::setup(void)
{
    // preallocate memory for meson fields
    int nt = env().getDim(Tp);

    envTmp(std::vector<PropagatorField>, "WWVV_LL", 1, nt, PropagatorField(env().getGrid()));
    envTmp(std::vector<PropagatorField>, "WWVV_LS", 1, nt, PropagatorField(env().getGrid()));
    envTmp(std::vector<PropagatorField>, "WWVV_LLS", 1, nt, PropagatorField(env().getGrid()));
    envTmp(PropagatorField, "VW_U", 1, PropagatorField(env().getGrid()));
    envTmp(PropagatorField, "VW_C", 1, PropagatorField(env().getGrid()));
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2AWeakHamiltonian<FImpl>::execute(void)
{
    auto &vStrange = envGet(std::vector<FermionField>, par().vStrange);
    auto &vLight1 = envGet(std::vector<FermionField>, par().vLight1);
    auto &wLight1 = envGet(std::vector<FermionField>, par().wLight1);
    auto &vLight2 = envGet(std::vector<FermionField>, par().vLight2);
    auto &vCharm = envGet(std::vector<FermionField>, par().vCharm);
    auto &wCharm = envGet(std::vector<FermionField>, par().wCharm);

    int dt_min = par().dt_min;
    int dt_max = par().dt_max;

    int nt = env().getDim(Tp);

    int N_i = vLight1.size();
    int N_j = vStrange.size();
    GridBase *grid = vLight1[0]._grid;
    int nd = grid->_ndimension;
    double vol = 1.0;
    for (int dim = 0; dim < nd; dim++)
    {
        vol = vol * grid->GlobalDimensions()[dim];
    }
    auto G5 = Gamma(Gamma::Algebra::Gamma5);

    Result result;
    result.dt_min = dt_min;
    result.dt_max = dt_max;
    result.corr_wing.resize(nt);
    result.corr_conncected.resize(nt);
    result.corr_eye_light.resize(nt);
    result.corr_saucer_light.resize(nt);
    result.corr_eye_charm.resize(nt);
    result.corr_saucer_charm.resize(nt);

    LOG(Message) << "Computing weak hamiltonian non-eye contraction using A2A" << std::endl;
    LOG(Message) << "Vs: '" << par().vStrange << "'" << std::endl;
    LOG(Message) << "Wl_1: '" << par().wLight1 << "' Vl_1: '" << par().vLight1 << "'" << std::endl;
    LOG(Message) << "Vl_2: '" << par().vLight2 << "'" << std::endl;
    LOG(Message) << "Wc: '" << par().wCharm << "' Vc: '" << par().vCharm << "'" << std::endl;

    LOG(Message) << "Meson Field light/strange size: " << nt << "*" << N_i << "*" << N_j
                 << " (filesize " << sizeString(nt * N_i * N_j * sizeof(Complex))
                 << "/momentum/bilinear)" << std::endl;

    LOG(Message) << "Meson Field light/light size: " << nt << "*" << N_i << "*" << N_i
                 << " (filesize " << sizeString(nt * N_i * N_i * sizeof(Complex))
                 << "/momentum/bilinear)" << std::endl;

    std::string WWLL_dv_file = "eigendiskvector_LL";
    std::string WVLL_G5_dv_file = "eigendiskvector_LL_G5";
    std::string WWLS_dv_file = "eigendiskvector_LS";
    std::string LLS_dv_file = "eigendiskvector_LLS";

    EigenDiskVector<ComplexD> pionFieldWW_LL_dv(WWLL_dv_file, nt, 1, true, grid);
    EigenDiskVector<ComplexD> pionFieldWV_LL_G5_dv(WVLL_G5_dv_file, nt, 1, true, grid);
    EigenDiskVector<ComplexD> pionFieldWW_LS_dv(WWLS_dv_file, nt, 1, true, grid);
    EigenDiskVector<ComplexD> pionField_LLS_dv(LLS_dv_file, nt, 1, true, grid);

    int traj = vm().getTrajectory();

    std::string mfDirLL = par().mfDirLL;
    std::string mfDirLS = par().mfDirLS;
    std::string WWLL_file = mfDirLL + "/mf_ll_ww_12." + std::to_string(traj) + "/Identity_0_0_0.h5";
    std::string WWLL_dataset = "Identity_0_0_0";
    std::string WVLL_G5_file = mfDirLL + "/mf_ll_wv_21." + std::to_string(traj) + "/Gamma5_0_0_0.h5";
    std::string WVLL_G5_dataset = "Gamma5_0_0_0";
    std::string WWLS_file = mfDirLS + "/mf_ls." + std::to_string(traj) + "/Identity_0_0_0.h5";
    std::string WWLS_dataset = "Identity_0_0_0";
    double t;
    std::cout << "-- Loading '" << WWLL_file << "'..." << std::endl;
    A2AMatrixIo<HADRONS_A2AM_IO_TYPE> WWLL_a2aIo(WWLL_file, WWLL_dataset, nt, N_i, N_i);
    WWLL_a2aIo.load(pionFieldWW_LL_dv, &t, grid);
    std::cout << "Read " << WWLL_a2aIo.getSize() << " bytes in " << t << " usec, " << WWLL_a2aIo.getSize() / t * 1.0e6 / 1024 / 1024 << " MB/s" << std::endl;

    std::cout << "-- Loading '" << WVLL_G5_file << "'..." << std::endl;
    A2AMatrixIo<HADRONS_A2AM_IO_TYPE> WVLL_a2aIo(WVLL_G5_file, WVLL_G5_dataset, nt, N_i, N_i);
    WVLL_a2aIo.load(pionFieldWV_LL_G5_dv, &t, grid);
    std::cout << "Read " << WVLL_a2aIo.getSize() << " bytes in " << t << " usec, " << WVLL_a2aIo.getSize() / t * 1.0e6 / 1024 / 1024 << " MB/s" << std::endl;

    std::cout << "-- Loading '" << WWLS_file << "'..." << std::endl;
    A2AMatrixIo<HADRONS_A2AM_IO_TYPE> WWLS_a2aIo(WWLS_file, WWLS_dataset, nt, N_i, N_j);
    WWLS_a2aIo.load(pionFieldWW_LS_dv, &t, grid);
    std::cout << "Read " << WWLS_a2aIo.getSize() << " bytes in " << t << " usec, " << WWLS_a2aIo.getSize() / t * 1.0e6 / 1024 / 1024 << " MB/s" << std::endl;

    envGetTmp(std::vector<PropagatorField>, WWVV_LL);
    envGetTmp(std::vector<PropagatorField>, WWVV_LS);
    envGetTmp(std::vector<PropagatorField>, WWVV_LLS);
    envGetTmp(PropagatorField, VW_U);
    envGetTmp(PropagatorField, VW_C);

    Eigen::Matrix<ComplexD, Eigen::Dynamic, Eigen::Dynamic> PI_WVLL_G5_t1(N_i, N_i); // <ComplexD, Eigen::Dynamic, Eigen::Dynamic, Eigen ::RowMajor>
    Eigen::Matrix<ComplexD, Eigen::Dynamic, Eigen::Dynamic> PI_WW_LS_t0(N_i, N_j);

    std::cout << GridLogMessage << " dt " << dt_min << "..." << dt_max << std::endl;

    double matmult = -usecond();
    for (int t0 = 0; t0 < nt; t0++)
    {
        Eigen::Matrix<ComplexD, Eigen::Dynamic, Eigen::Dynamic> PIik;

        PIik = Eigen::Matrix<ComplexD, Eigen::Dynamic, Eigen::Dynamic>::Zero(N_i, N_j);
        std::cout << GridLogMessage << " t0 " << t0 << std::endl;
        // for (int dt = dt_min; dt < dt_max; dt++)
        {
            int dt = dt_min;
            int t1 = (t0 + dt) % nt;
            std::cout << GridLogMessage << " t1 " << t1 << std::endl;

            PI_WW_LS_t0 = pionFieldWW_LS_dv[t0];
            PI_WVLL_G5_t1 = pionFieldWV_LL_G5_dv[t1];

            PIik += PI_WVLL_G5_t1 * PI_WW_LS_t0;
        }
        pionField_LLS_dv[t0] = PIik;
    }
    matmult += usecond();

    double million = 1.0e6;
    std::cout << "Performing " << nt * (dt_max - dt_min) << " matrix multiplications: " << matmult / million << " s " << std::endl;

    for (int t = 0; t < nt; t++)
    {
        WWVV_LL[t] = zero;
        WWVV_LS[t] = zero;
        WWVV_LLS[t] = zero;
    }

    std::cout << "Computing contr1 " << std::endl;
    double contr1 = -usecond();
    A2A_DVutils<FImpl>::ContractWWVV_dv(WWVV_LL, pionFieldWW_LL_dv, &vLight1[0], &vLight2[0]);
    contr1 += usecond();

    std::cout << "Computing contr2 " << std::endl;
    double contr2 = -usecond();
    A2A_DVutils<FImpl>::ContractWWVV_dv(WWVV_LS, pionFieldWW_LS_dv, &vLight1[0], &vStrange[0]);
    contr2 += usecond();

    std::cout << "Computing contr3 " << std::endl;
    double contr3 = -usecond();
    A2A_DVutils<FImpl>::ContractWWVV_dv(WWVV_LLS, pionField_LLS_dv, &vLight2[0], &vStrange[0]);
    contr3 += usecond();
    A2Autils<FImpl>::ContractVW(VW_U, vLight1, wLight1);
    A2Autils<FImpl>::ContractVW(VW_C, vCharm, wCharm);

    std::cout << "Computing contr1    " << contr1 / million << " s " << std::endl;
    std::cout << "Computing contr2    " << contr2 / million << " s " << std::endl;
    std::cout << "Computing contr3    " << contr2 / million << " s " << std::endl;
    //////////////////////////////
    // Implicit gamma-5
    //////////////////////////////
    for (int t = 0; t < nt; t++)
    {
        WWVV_LL[t] = WWVV_LL[t] * G5;
        WWVV_LS[t] = WWVV_LS[t] * G5;
        WWVV_LLS[t] = WWVV_LLS[t] * G5;
    }

    ////////////////////////////////////////////////////////
    // Contraction
    ////////////////////////////////////////////////////////

    ComplexField c_field(grid);
    c_field = zero;
    ComplexField w_field(grid);
    w_field = zero;

    ComplexField c_field_VV(grid);
    c_field_VV = zero;
    ComplexField w_field_VV(grid);
    w_field_VV = zero;

    ComplexField c_field_AA(grid);
    c_field_AA = zero;
    ComplexField w_field_AA(grid);
    w_field_AA = zero;

    ComplexField s_field(grid);
    s_field = zero;
    ComplexField e_field(grid);
    e_field = zero;

    ComplexField s_field_VV(grid);
    s_field_VV = zero;
    ComplexField e_field_VV(grid);
    e_field_VV = zero;

    ComplexField s_field_AA(grid);
    s_field_AA = zero;
    ComplexField e_field_AA(grid);
    e_field_AA = zero;

    ComplexField s_c_field(grid);
    s_c_field = zero;
    ComplexField e_c_field(grid);
    e_c_field = zero;

    ComplexField s_c_field_VV(grid);
    s_c_field_VV = zero;
    ComplexField e_c_field_VV(grid);
    e_c_field_VV = zero;

    ComplexField s_c_field_AA(grid);
    s_c_field_AA = zero;
    ComplexField e_c_field_AA(grid);
    e_c_field_AA = zero;

    //////////////////////////////////////////////////
    // Used to store appropriate correlation funcs
    //////////////////////////////////////////////////
    std::vector<TComplex> C1;

    double t_contr = -usecond();

    // Wick contraction
    auto VX = Gamma(Gamma::Algebra::GammaX);
    auto VY = Gamma(Gamma::Algebra::GammaY);
    auto VZ = Gamma(Gamma::Algebra::GammaZ);
    auto VT = Gamma(Gamma::Algebra::GammaT);

    auto AX = Gamma(Gamma::Algebra::GammaXGamma5);
    auto AY = Gamma(Gamma::Algebra::GammaYGamma5);
    auto AZ = Gamma(Gamma::Algebra::GammaZGamma5);
    auto AT = Gamma(Gamma::Algebra::GammaTGamma5);

    std::vector<Gamma> VV({VX, VY, VZ, VT});
    std::vector<Gamma> AA({AX, AY, AZ, AT});

    std::cout << GridLogMessage << " dt " << dt_min << "..." << dt_max << std::endl;
    for (int t0 = 0; t0 < nt; t0++)
    {
        std::cout << GridLogMessage << " t0 " << t0 << std::endl;
        // for (int dt = dt_min; dt < dt_max; dt++)
        {
            int dt = dt_min;

            int t1 = (t0 + dt) % nt;
            std::cout << GridLogMessage << " t1 " << t1 << std::endl;

            A2Autils<FImpl>::ContractFourQuarkColourDiagonal(WWVV_LS[t0], WWVV_LL[t1], VV, VV, c_field_VV, w_field_VV); // VV
            A2Autils<FImpl>::ContractFourQuarkColourDiagonal(WWVV_LS[t0], WWVV_LL[t1], AA, AA, c_field_AA, w_field_AA); // AA

            c_field = c_field_VV - c_field_AA;
            w_field = w_field_VV - w_field_AA;

            sliceSum(c_field, C1, Tp);

            for (int t = 0; t < nt; t++)
            {
                result.corr_conncected[t] += C1[(t + t0) % nt]()()() / vol;
            }

            sliceSum(w_field, C1, Tp);

            for (int t = 0; t < nt; t++)
            {
                result.corr_wing[t] += C1[(t + t0) % nt]()()() / vol;
            }
        }

        A2Autils<FImpl>::ContractFourQuarkColourDiagonal(WWVV_LLS[t0], VW_U, VV, VV, s_field_VV, e_field_VV); // VV
        A2Autils<FImpl>::ContractFourQuarkColourDiagonal(WWVV_LLS[t0], VW_U, AA, AA, s_field_AA, e_field_AA); // AA

        s_field = s_field_VV - s_field_AA;
        e_field = e_field_VV - e_field_AA;

        sliceSum(s_field, C1, Tp);

        for (int t = 0; t < nt; t++)
        {
            result.corr_saucer_light[t] += C1[(t + t0) % nt]()()() / vol;
        }

        sliceSum(e_field, C1, Tp);

        for (int t = 0; t < nt; t++)
        {
            result.corr_eye_light[t] += C1[(t + t0) % nt]()()() / vol;
        }

        A2Autils<FImpl>::ContractFourQuarkColourDiagonal(WWVV_LLS[t0], VW_C, VV, VV, s_c_field_VV, e_c_field_VV); // VV
        A2Autils<FImpl>::ContractFourQuarkColourDiagonal(WWVV_LLS[t0], VW_C, AA, AA, s_c_field_AA, e_c_field_AA); // AA

        s_c_field = s_c_field_VV - s_c_field_AA;
        e_c_field = e_c_field_VV - e_c_field_AA;

        sliceSum(s_c_field, C1, Tp);

        for (int t = 0; t < nt; t++)
        {
            result.corr_saucer_charm[t] += C1[(t + t0) % nt]()()() / vol;
        }

        sliceSum(e_c_field, C1, Tp);

        for (int t = 0; t < nt; t++)
        {
            result.corr_eye_charm[t] += C1[(t + t0) % nt]()()() / vol;
        }
    }

    LOG(Message) << "Writing results to " << par().output << "." << std::endl;
    saveResult(par().output, "results", result);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_A2AWeakHamiltonian_hpp_
