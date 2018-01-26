#include <Grid/Hadrons/Modules/MScalar/ChargedProp.hpp>
#include <Grid/Hadrons/Modules/MScalar/ScalarVP.hpp>
#include <Grid/Hadrons/Modules/MScalar/Scalar.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MScalar;

/*
 * Scalar QED vacuum polarisation up to O(alpha)
 *
 * Conserved vector 2-point function diagram notation:
 *        _______
 *       /       \
 * U_nu *         * U_mu
 *       \_______/
 *
 *                (   adj(S(a\hat{nu}|x)) U_mu(x) S(0|x+a\hat{mu}) U_nu(0)    )
 *          = 2 Re(                             -                             )
 *                ( adj(S(a\hat{nu}|x+a\hat{mu})) adj(U_mu(x)) S(0|x) U_nu(0) )
 *  
 *
 *            _______
 *           /       \
 * free = 1 *         * 1
 *           \_______/
 *
 *
 *
 *             _______
 *            /       \
 * S = iA_nu *         * iA_mu
 *            \_______/
 *
 *
 *         Delta_1
 *         ___*___
 *        /       \
 * X = 1 *         * 1
 *        \___*___/
 *         Delta_1
 *
 *          Delta_1                     Delta_1
 *          ___*___                     ___*___
 *         /       \                   /       \
 *      1 *         * iA_mu  +  iA_nu *         * 1
 *         \_______/                   \_______/
 * 4C =        _______                     _______
 *            /       \                   /       \
 *      +  1 *         * iA_mu  +  iA_nu *         * 1
 *            \___*___/                   \___*___/
 *             Delta_1                     Delta_1
 *
 *     Delta_1   Delta_1
 *          _*___*_             _______
 *         /       \           /       \
 * 2E = 1 *         * 1  +  1 *         * 1
 *         \_______/           \_*___*_/
 *                         Delta_1   Delta_1
 *
 *          Delta_2
 *          ___*___             _______
 *         /       \           /       \
 * 2T = 1 *         * 1  +  1 *         * 1
 *         \_______/           \___*___/
 *                              Delta_2
 *
 *
 *                    _______
 *                   /       \
 * srcT = -A_nu^2/2 *         * 1
 *                   \_______/
 *
 *
 *
 *            _______
 *           /       \
 * snkT = 1 *         * -A_mu^2/2
 *           \_______/
 *
 * Full VP to O(alpha) = free + q^2*(S+X+4C+2E+2T+srcT+snkT)
 */

/******************************************************************************
*                  TScalarVP implementation                             *
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
TScalarVP::TScalarVP(const std::string name)
: Module<ScalarVPPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
std::vector<std::string> TScalarVP::getInput(void)
{
    prop0Name_ = par().scalarProp + "_0";
    propQName_ = par().scalarProp + "_Q";
    propSunName_ = par().scalarProp + "_Sun";
    propTadName_ = par().scalarProp + "_Tad";

	std::vector<std::string> in = {par().emField, prop0Name_, propQName_,
                                   propSunName_, propTadName_};
    
    return in;
}

std::vector<std::string> TScalarVP::getOutput(void)
{
    std::vector<std::string> out;
    
    for (unsigned int mu = 0; mu < env().getNd(); ++mu)
    {
        // out.push_back(getName() + "_propQ_" + std::to_string(mu));

        for (unsigned int nu = 0; nu < env().getNd(); ++nu)
        {
            out.push_back(getName() + "_" + std::to_string(mu)
                          + "_" + std::to_string(nu));
        }
    }

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
void TScalarVP::setup(void)
{
	freeMomPropName_ = FREEMOMPROP(static_cast<TChargedProp *>(vm().getModule(par().scalarProp))->par().mass);
	GFSrcName_ = par().scalarProp + "_DinvSrc";
    fftName_   = par().scalarProp + "_fft";
	phaseName_.clear();
	muPropQName_.clear();
    vpTensorName_.clear();
    momPhaseName_.clear();
	for (unsigned int mu = 0; mu < env().getNd(); ++mu)
    {
        phaseName_.push_back("_shiftphase_" + std::to_string(mu));
        muPropQName_.push_back(getName() + "_propQ_" + std::to_string(mu));

        std::vector<std::string> vpTensorName_mu;
        for (unsigned int nu = 0; nu < env().getNd(); ++nu)
        {
            vpTensorName_mu.push_back(getName() + "_" + std::to_string(mu)
                                      + "_" + std::to_string(nu));
        }
        vpTensorName_.push_back(vpTensorName_mu);
    }
    if (!par().output.empty())
    {
        for (unsigned int i_p = 0; i_p < par().outputMom.size(); ++i_p)
        {
            momPhaseName_.push_back("_momentumphase_" + std::to_string(i_p));
        }
    }

    for (unsigned int mu = 0; mu < env().getNd(); ++mu)
	{
	    envCreateLat(ScalarField, muPropQName_[mu]);

        for (unsigned int nu = 0; nu < env().getNd(); ++nu)
        {
            envCreateLat(ScalarField, vpTensorName_[mu][nu]);
        }
	}
    if (!par().output.empty())
    {
        momPhasesDone_ = env().hasCreatedObject(momPhaseName_[0]);
        for (unsigned int i_p = 0; i_p < par().outputMom.size(); ++i_p)
        {
            envCacheLat(ScalarField, momPhaseName_[i_p]);
        }
    }
    envTmpLat(ScalarField, "buf");
    envTmpLat(ScalarField, "result");
    envTmpLat(ScalarField, "Amu");
    envTmpLat(ScalarField, "Usnk");
    envTmpLat(ScalarField, "tmpProp");
}

// execution ///////////////////////////////////////////////////////////////////
void TScalarVP::execute(void)
{
    // CACHING ANALYTIC EXPRESSIONS
    makeCaches();

    Complex ci(0.0,1.0);
    Real    q        = static_cast<TChargedProp *>(vm().getModule(par().scalarProp))->par().charge;
    auto    &prop0   = envGet(ScalarField, prop0Name_);
    auto    &propQ   = envGet(ScalarField, propQName_);
    auto    &propSun = envGet(ScalarField, propSunName_);
    auto    &propTad = envGet(ScalarField, propTadName_);
    auto    &GFSrc   = envGet(ScalarField, GFSrcName_);
    auto    &G       = envGet(ScalarField, freeMomPropName_);
    auto    &fft     = envGet(FFT, fftName_);
    phase_.clear();
    for (unsigned int mu = 0; mu < env().getNd(); ++mu)
    {
        auto &phmu = envGet(ScalarField, phaseName_[mu]);
        phase_.push_back(&phmu);
    }
    
    // PROPAGATORS FROM SHIFTED SOURCES
    LOG(Message) << "Computing O(q) charged scalar propagators..."
                 << std::endl;
    std::vector<ScalarField *> muPropQ;
    for (unsigned int mu = 0; mu < env().getNd(); ++mu)
    {
        auto &propmu = envGet(ScalarField, muPropQName_[mu]);

        // -G*momD1*G*F*tau_mu*Src (momD1 = F*D1*Finv)
        propmu = adj(*phase_[mu])*GFSrc;
        momD1(propmu, fft);
        propmu = -G*propmu;
        fft.FFT_all_dim(propmu, propmu, FFT::backward);

        muPropQ.push_back(&propmu);
    }

    // CONTRACTIONS
    auto        &A = envGet(EmField, par().emField);
    envGetTmp(ScalarField, buf);
    envGetTmp(ScalarField, result);
    envGetTmp(ScalarField, Amu);
    envGetTmp(ScalarField, Usnk);
    envGetTmp(ScalarField, tmpProp);
    TComplex    Anu0, Usrc;
    std::vector<int> coor0 = {0, 0, 0, 0};
    std::vector<std::vector<ScalarField *> > vpTensor;
    for (unsigned int mu = 0; mu < env().getNd(); ++mu)
    {
        std::vector<ScalarField *> vpTensor_mu;
        for (unsigned int nu = 0; nu < env().getNd(); ++nu)
        {
            auto &vpmunu = envGet(ScalarField, vpTensorName_[mu][nu]);
            vpTensor_mu.push_back(&vpmunu);
        }
        vpTensor.push_back(vpTensor_mu);
    }

    // Open output files if necessary
    std::vector<ResultWriter *> writer;
    if (!par().output.empty())
    {
        LOG(Message) << "Preparing output files..." << std::endl;
        for (unsigned int i_p = 0; i_p < par().outputMom.size(); ++i_p)
        {
            std::vector<int> mom = strToVec<int>(par().outputMom[i_p]);
            std::string filename = par().output + "_"
                                   + std::to_string(mom[0])
                                   + std::to_string(mom[1])
                                   + std::to_string(mom[2]);
            ResultWriter *writer_i = new ResultWriter(RESULT_FILE_NAME(filename));
            writer.push_back(writer_i);
            write(*writer_i, "charge", q);
            write(*writer_i, "mass", static_cast<TChargedProp *>(vm().getModule(par().scalarProp))->par().mass);
        }
    }

    // Do contractions
    for (unsigned int nu = 0; nu < env().getNd(); ++nu)
    {
        peekSite(Anu0, peekLorentz(A, nu), coor0);

        for (unsigned int mu = 0; mu < env().getNd(); ++mu)
        {
            LOG(Message) << "Computing Pi[" << mu << "][" << nu << "]..."
                         << std::endl;
            Amu = peekLorentz(A, mu);

            // free
            tmpProp = Cshift(prop0, nu, -1);     // S_0(0|x-a\hat{\nu})
                                                 // = S_0(a\hat{\nu}|x)
            Usrc    = Complex(1.0,0.0);
            vpContraction(buf, prop0, tmpProp, Usrc, mu);
            *vpTensor[mu][nu] = buf;
            // Output if necessary
            if (!par().output.empty())
            {
                writeVP(writer, buf,
                        "Pi_free_"+std::to_string(mu)+"_"+std::to_string(nu));
            }

            // srcT
            result = buf * (-0.5)*q*q*Anu0*Anu0;
            *vpTensor[mu][nu] += result;
            // Output if necessary
            if (!par().output.empty())
            {
                writeVP(writer, result,
                        "Pi_srcT_"+std::to_string(mu)+"_"+std::to_string(nu));
            }

            // snkT
            result = buf * (-0.5)*q*q*Amu*Amu;
            *vpTensor[mu][nu] += result;
            // Output if necessary
            if (!par().output.empty())
            {
                writeVP(writer, result,
                        "Pi_snkT_"+std::to_string(mu)+"_"+std::to_string(nu));
            }

            // S
            tmpProp = Cshift(prop0, nu, -1);     // S_0(a\hat{\nu}|x)
            Usrc    = ci*Anu0;
            Usnk    = ci*Amu;
            vpContraction(result, prop0, tmpProp, Usrc, Usnk, mu);
            result = q*q*result;
            *vpTensor[mu][nu] += result;
            // Output if necessary
            if (!par().output.empty())
            {
                writeVP(writer, result,
                        "Pi_S_"+std::to_string(mu)+"_"+std::to_string(nu));
            }

            // 4C
            tmpProp = Cshift(prop0, nu, -1);     // S_0(a\hat{\nu}|x)
            Usrc    = Complex(1.0,0.0);
            Usnk    = ci*Amu;
            vpContraction(result, propQ, tmpProp, Usrc, Usnk, mu);
            Usrc    = ci*Anu0;
            vpContraction(buf, propQ, tmpProp, Usrc, mu);
            result += buf;
            vpContraction(buf, prop0, *muPropQ[nu], Usrc, mu);
            result += buf;
            Usrc = Complex(1.0,0.0);
            Usnk = ci*Amu;
            vpContraction(buf, prop0, *muPropQ[nu], Usrc, Usnk, mu);
            result += buf;
            result = q*q*result;
            *vpTensor[mu][nu] += result;
            // Output if necessary
            if (!par().output.empty())
            {
                writeVP(writer, result,
                        "Pi_4C_"+std::to_string(mu)+"_"+std::to_string(nu));
            }

            // X
            Usrc = Complex(1.0,0.0);
            vpContraction(result, propQ, *muPropQ[nu], Usrc, mu);
            result = q*q*result;
            *vpTensor[mu][nu] += result;
            // Output if necessary
            if (!par().output.empty())
            {
                writeVP(writer, result,
                        "Pi_X_"+std::to_string(mu)+"_"+std::to_string(nu));
            }

            // 2E
            tmpProp = Cshift(prop0, nu, -1);     // S_0(a\hat{\nu}|x)
            Usrc    = Complex(1.0,0.0);
            vpContraction(result, propSun, tmpProp, Usrc, mu);
            tmpProp = Cshift(propSun, nu, -1);     // S_\Sigma(0|x-a\hat{\nu})
                               //(Note: <S(0|x-a\hat{\nu})> = <S(a\hat{\nu}|x)>)
            vpContraction(buf, prop0, tmpProp, Usrc, mu);
            result += buf;
            result = q*q*result;
            *vpTensor[mu][nu] += result;
            // Output if necessary
            if (!par().output.empty())
            {
                writeVP(writer, result,
                        "Pi_2E_"+std::to_string(mu)+"_"+std::to_string(nu));
            }

            // 2T
            tmpProp = Cshift(prop0, nu, -1);     // S_0(a\hat{\nu}|x)
            Usrc    = Complex(1.0,0.0);
            vpContraction(result, propTad, tmpProp, Usrc, mu);
            tmpProp = Cshift(propTad, nu, -1);     // S_T(0|x-a\hat{\nu})
            vpContraction(buf, prop0, tmpProp, Usrc, mu);
            result += buf;
            result = q*q*result;
            *vpTensor[mu][nu] += result;
            // Output if necessary
            if (!par().output.empty())
            {
                writeVP(writer, result,
                        "Pi_2T_"+std::to_string(mu)+"_"+std::to_string(nu));
            }

            // Output full VP if necessary
            if (!par().output.empty())
            {
                writeVP(writer, *vpTensor[mu][nu],
                        "Pi_"+std::to_string(mu)+"_"+std::to_string(nu));
            }
        }
    }
    if (!par().output.empty())
    {
        for (unsigned int i_p = 0; i_p < par().outputMom.size(); ++i_p)
        {
            delete writer[i_p];
        }
    }
}

void TScalarVP::makeCaches(void)
{
    envGetTmp(ScalarField, buf);

    if ( (!par().output.empty()) && (!momPhasesDone_) )
    {
        LOG(Message) << "Caching phases for momentum projections..."
                     << std::endl;
        std::vector<int> &l = env().getGrid()->_fdimensions;
        Complex          ci(0.0,1.0);

        // Calculate phase factors
        for (unsigned int i_p = 0; i_p < par().outputMom.size(); ++i_p)
        {
            std::vector<int> mom = strToVec<int>(par().outputMom[i_p]);
            auto &momph_ip = envGet(ScalarField, momPhaseName_[i_p]);
            momph_ip = zero;
            for (unsigned int j = 0; j < env().getNd()-1; ++j)
            {
                Real twoPiL = M_PI*2./l[j];
                LatticeCoordinate(buf, j);
                buf = mom[j]*twoPiL*buf;
                momph_ip = momph_ip + buf;
            }
            momph_ip = exp(-ci*momph_ip);
            momPhase_.push_back(&momph_ip);
        }
    }
}

void TScalarVP::vpContraction(ScalarField &vp,
                   ScalarField &prop_0_x, ScalarField &prop_nu_x,
                   TComplex u_src, ScalarField &u_snk, int mu)
{
    // Note: this function assumes a point source is used.
    vp = adj(prop_nu_x) * u_snk * Cshift(prop_0_x, mu, 1) * u_src;
    vp -= Cshift(adj(prop_nu_x), mu, 1) * adj(u_snk) * prop_0_x * u_src;
    vp = 2.0*real(vp);
}

void TScalarVP::vpContraction(ScalarField &vp,
                   ScalarField &prop_0_x, ScalarField &prop_nu_x,
                   TComplex u_src, int mu)
{
    // Note: this function assumes a point source is used.
    vp = adj(prop_nu_x) * Cshift(prop_0_x, mu, 1) * u_src;
    vp -= Cshift(adj(prop_nu_x), mu, 1) * prop_0_x * u_src;
    vp = 2.0*real(vp);
}

void TScalarVP::writeVP(std::vector<ResultWriter *> &writer, const ScalarField &vp, std::string dsetName)
{
    std::vector<TComplex>   vecBuf;
    std::vector<Complex>    result;
    envGetTmp(ScalarField, buf);

    for (unsigned int i_p = 0; i_p < par().outputMom.size(); ++i_p)
    {
        buf = vp*(*momPhase_[i_p]);
        sliceSum(buf, vecBuf, Tp);
        result.resize(vecBuf.size());
        for (unsigned int t = 0; t < vecBuf.size(); ++t)
        {
            result[t] = TensorRemove(vecBuf[t]);
        }
        write(*writer[i_p], dsetName, result);
    }
}

void TScalarVP::momD1(ScalarField &s, FFT &fft)
{
    auto        &A = envGet(EmField, par().emField);
    Complex     ci(0.0,1.0);

    envGetTmp(ScalarField, buf);
    envGetTmp(ScalarField, result);
    envGetTmp(ScalarField, Amu);

    result = zero;
    for (unsigned int mu = 0; mu < env().getNd(); ++mu)
    {
        Amu = peekLorentz(A, mu);
        buf = (*phase_[mu])*s;
        fft.FFT_all_dim(buf, buf, FFT::backward);
        buf = Amu*buf;
        fft.FFT_all_dim(buf, buf, FFT::forward);
        result = result - ci*buf;
    }
    fft.FFT_all_dim(s, s, FFT::backward);
    for (unsigned int mu = 0; mu < env().getNd(); ++mu)
    {
        Amu = peekLorentz(A, mu);
        buf = Amu*s;
        fft.FFT_all_dim(buf, buf, FFT::forward);
        result = result + ci*adj(*phase_[mu])*buf;
    }

    s = result;
}
