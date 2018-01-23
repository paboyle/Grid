#include <Grid/Hadrons/Modules/MScalar/ChargedProp.hpp>
#include <Grid/Hadrons/Modules/MScalar/ScalarVP.hpp>
#include <Grid/Hadrons/Modules/MScalar/Scalar.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MScalar;

/*
 * Scalar QED vacuum polarisation up to O(alpha)
 *
 *
 *                          _______
 *                         /       \              (   adj(S(a\hat{nu}|x)) U_mu(x) S(0|x+a\hat{mu}) U_nu(0)    )
 * Diagram notation: U_nu *         * U_mu  = 2 Re(                             -                             )
 *                         \_______/              ( adj(S(a\hat{nu}|x+a\hat{mu})) adj(U_mu(x)) S(0|x) U_nu(0) )
 *                          
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
    propQName_ = par().scalarProp + "_Q";
    propSunName_ = par().scalarProp + "_Sun";
    propTadName_ = par().scalarProp + "_Tad";

	std::vector<std::string> in = {par().emField, propQName_, propSunName_,
                                   propTadName_};
    
    return in;
}

std::vector<std::string> TScalarVP::getOutput(void)
{
    std::vector<std::string> out;
    
    for (unsigned int mu = 0; mu < env().getNd(); ++mu)
    {
        out.push_back(getName() + "_propQ_" + std::to_string(mu));

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
	freeMomPropName_ = FREEMOMPROP(static_cast<TChargedProp *>(env().getModule(par().scalarProp))->par().mass);
	GFSrcName_ = "_" + par().scalarProp + "_DinvSrc";
    prop0Name_ = par().scalarProp + "_0";

	phaseName_.clear();
	muPropQName_.clear();
    vpTensorName_.clear();
    
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

    for (unsigned int mu = 0; mu < env().getNd(); ++mu)
	{
	    env().registerLattice<ScalarField>(muPropQName_[mu]);

        for (unsigned int nu = 0; nu < env().getNd(); ++nu)
        {
            env().registerLattice<ScalarField>(vpTensorName_[mu][nu]);
        }
	}
}

// execution ///////////////////////////////////////////////////////////////////
void TScalarVP::execute(void)
{
    // Get objects cached by ChargedProp module
    Complex     ci(0.0,1.0);
    FFT         fft(env().getGrid());
    Real        q = static_cast<TChargedProp *>(env().getModule(par().scalarProp))->par().charge;

    freeMomProp_ = env().getObject<ScalarField>(freeMomPropName_);
    for (unsigned int mu = 0; mu < env().getNd(); ++mu)
    {
        phase_.push_back(env().getObject<ScalarField>(phaseName_[mu]));
    }
    GFSrc_ = env().getObject<ScalarField>(GFSrcName_);
    prop0_ = env().getObject<ScalarField>(prop0Name_);

    // Propagator from unshifted source
	ScalarField &propQ   = *env().getObject<ScalarField>(propQName_);
    ScalarField &propSun   = *env().getObject<ScalarField>(propSunName_);
    ScalarField &propTad   = *env().getObject<ScalarField>(propTadName_);

    // Propagators from shifted sources
    LOG(Message) << "Computing O(q) charged scalar propagators..."
                 << std::endl;
    std::vector<ScalarField> muPropQ;
    for (unsigned int mu = 0; mu < env().getNd(); ++mu)
    {
        muPropQ.push_back(*env().createLattice<ScalarField>(muPropQName_[mu]));

        // -G*momD1*G*F*tau_mu*Src (momD1 = F*D1*Finv)
        muPropQ[mu] = adj(*phase_[mu])*(*GFSrc_);
        momD1(muPropQ[mu], fft);
        muPropQ[mu] = -(*freeMomProp_)*muPropQ[mu];
        fft.FFT_all_dim(muPropQ[mu], muPropQ[mu], FFT::backward);
    }

    // CONTRACTIONS
    ScalarField prop1(env().getGrid()), prop2(env().getGrid());
    EmField     &A = *env().getObject<EmField>(par().emField);
    ScalarField Amu(env().getGrid()), U_snk(env().getGrid());
    ScalarField tmp_vp1(env().getGrid()), tmp_vp2(env().getGrid());
    TComplex    Anu0, U_src;
    std::vector<int> coor0 = {0, 0, 0, 0};
    std::vector<std::vector<ScalarField> > vpTensor;
    for (unsigned int mu = 0; mu < env().getNd(); ++mu)
    {
        std::vector<ScalarField> vpTensor_mu;
        for (unsigned int nu = 0; nu < env().getNd(); ++nu)
        {
            vpTensor_mu.push_back(*env().createLattice<ScalarField>(vpTensorName_[mu][nu]));
        }
        vpTensor.push_back(vpTensor_mu);
    }

    // Open output files if necessary
    std::vector<CorrWriter *> writer;
    std::vector<ScalarField> momphases;
    if (!par().output.empty())
    {
        LOG(Message) << "Preparing output files..." << std::endl;
        for (unsigned int i_p = 0; i_p < par().outputMom.size(); ++i_p)
        {
            std::vector<int> mom = strToVec<int>(par().outputMom[i_p]);

            // Open output files
            std::string filename = par().output + "_"
                                   + std::to_string(mom[0])
                                   + std::to_string(mom[1])
                                   + std::to_string(mom[2]) + "."
                                   + std::to_string(env().getTrajectory());

            if (env().getGrid()->IsBoss())
            {
                CorrWriter *writer_i = new CorrWriter(filename);
                writer.push_back(writer_i);

                write(*writer[i_p], "charge", q);
                write(*writer[i_p], "mass", static_cast<TChargedProp *>(env().getModule(par().scalarProp))->par().mass);
            }

            // Calculate phase factors
            tmp_vp1 = Complex(1.0,0.0);
            for (unsigned int j = 0; j < env().getNd()-1; ++j)
            {
                for (unsigned int momcount = 0; momcount < mom[j]; ++momcount)
                {
                    tmp_vp1 = tmp_vp1*(*phase_[j]);
                }
            }
            tmp_vp1 = adj(tmp_vp1);
            momphases.push_back(tmp_vp1);
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
            prop1 = *prop0_;                     // S_0(0|x)
            prop2 = Cshift(*prop0_, nu, -1);     // S_0(0|x-a\hat{\nu})
                                                 // = S_0(a\hat{\nu}|x)
            U_src = Complex(1.0,0.0);
            vpContraction(tmp_vp1, prop1, prop2, U_src, mu);
            vpTensor[mu][nu] = tmp_vp1;
            // Output if necessary
            if (!par().output.empty())
            {
                writeVP(writer, tmp_vp1, momphases,
                        "Pi_free_"+std::to_string(mu)+"_"+std::to_string(nu));
            }

            // srcT
            tmp_vp2 = tmp_vp1 * (-0.5)*q*q*Anu0*Anu0;
            vpTensor[mu][nu] += tmp_vp2;
            // Output if necessary
            if (!par().output.empty())
            {
                writeVP(writer, tmp_vp2, momphases,
                        "Pi_srcT_"+std::to_string(mu)+"_"+std::to_string(nu));
            }

            // snkT
            tmp_vp2 = tmp_vp1 * (-0.5)*q*q*Amu*Amu;
            vpTensor[mu][nu] += tmp_vp2;
            // Output if necessary
            if (!par().output.empty())
            {
                writeVP(writer, tmp_vp2, momphases,
                        "Pi_snkT_"+std::to_string(mu)+"_"+std::to_string(nu));
            }

            // S
            prop1 = *prop0_;                     // S_0(0|x)
            prop2 = Cshift(*prop0_, nu, -1);     // S_0(a\hat{\nu}|x)
            U_src = ci*q*Anu0;
            U_snk = ci*q*Amu;
            vpContraction(tmp_vp1, prop1, prop2, U_src, U_snk, mu);
            vpTensor[mu][nu] += tmp_vp1;
            // Output if necessary
            if (!par().output.empty())
            {
                writeVP(writer, tmp_vp1, momphases,
                        "Pi_S_"+std::to_string(mu)+"_"+std::to_string(nu));
            }

            // 4C
            prop1 = q*propQ;                     // q*S_1(0|x)
            prop2 = Cshift(*prop0_, nu, -1);     // S_0(a\hat{\nu}|x)
            U_src = Complex(1.0,0.0);
            U_snk = ci*q*Amu;
            vpContraction(tmp_vp1, prop1, prop2, U_src, U_snk, mu);
            U_src = ci*q*Anu0;
            vpContraction(tmp_vp2, prop1, prop2, U_src, mu);
            tmp_vp1 += tmp_vp2;
            prop1 = *prop0_;                     // S_0(0|x)
            prop2 = q*muPropQ[nu];               // q*S_1(a\hat{\nu}|x)
            vpContraction(tmp_vp2, prop1, prop2, U_src, mu);
            tmp_vp1 += tmp_vp2;
            U_src = Complex(1.0,0.0);
            U_snk = ci*q*Amu;
            vpContraction(tmp_vp2, prop1, prop2, U_src, U_snk, mu);
            tmp_vp1 += tmp_vp2;
            vpTensor[mu][nu] += tmp_vp1;
            // Output if necessary
            if (!par().output.empty())
            {
                writeVP(writer, tmp_vp1, momphases,
                        "Pi_4C_"+std::to_string(mu)+"_"+std::to_string(nu));
            }

            // X
            prop1 = q*propQ;                     // q*S_1(0|x)
            prop2 = q*muPropQ[nu];               // q*S_1(a\hat{\nu}|x)
            U_src = Complex(1.0,0.0);
            vpContraction(tmp_vp1, prop1, prop2, U_src, mu);
            vpTensor[mu][nu] += tmp_vp1;
            // Output if necessary
            if (!par().output.empty())
            {
                writeVP(writer, tmp_vp1, momphases,
                        "Pi_X_"+std::to_string(mu)+"_"+std::to_string(nu));
            }

            // 2E
            prop1 = q*q*propSun;                 // q^2*S_\Sigma(0|x)
            prop2 = Cshift(*prop0_, nu, -1);     // S_0(a\hat{\nu}|x)
            U_src = Complex(1.0,0.0);
            vpContraction(tmp_vp1, prop1, prop2, U_src, mu);
            prop1 = *prop0_;                     // S_0(0|x)
            prop2 = q*q*Cshift(propSun, nu, -1); // q^2*S_\Sigma(0|x-a\hat{\nu})
                               //(Note: <S(0|x-a\hat{\nu})> = <S(a\hat{\nu}|x)>)
            vpContraction(tmp_vp2, prop1, prop2, U_src, mu);
            tmp_vp1 += tmp_vp2;
            vpTensor[mu][nu] += tmp_vp1;
            // Output if necessary
            if (!par().output.empty())
            {
                writeVP(writer, tmp_vp1, momphases,
                        "Pi_2E_"+std::to_string(mu)+"_"+std::to_string(nu));
            }

            // 2T
            prop1 = q*q*propTad;                 // q^2*S_T(0|x)
            prop2 = Cshift(*prop0_, nu, -1);     // S_0(a\hat{\nu}|x)
            U_src = Complex(1.0,0.0);
            vpContraction(tmp_vp1, prop1, prop2, U_src, mu);
            prop1 = *prop0_;                     // S_0(0|x)
            prop2 = q*q*Cshift(propTad, nu, -1); // q^2*S_T(0|x-a\hat{\nu})
            vpContraction(tmp_vp2, prop1, prop2, U_src, mu);
            tmp_vp1 += tmp_vp2;
            vpTensor[mu][nu] += tmp_vp1;
            // Output if necessary
            if (!par().output.empty())
            {
                writeVP(writer, tmp_vp1, momphases,
                        "Pi_2T_"+std::to_string(mu)+"_"+std::to_string(nu));
            }

            // Output full VP if necessary
            if (!par().output.empty())
            {
                writeVP(writer, vpTensor[mu][nu], momphases,
                        "Pi_"+std::to_string(mu)+"_"+std::to_string(nu));
            }
        }
    }
    if (!par().output.empty())
    {
        for (unsigned int i_p = 0; i_p < par().outputMom.size(); ++i_p)
        {
            if (env().getGrid()->IsBoss())
            {
                delete writer[i_p];
            }
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

void TScalarVP::writeVP(const std::vector<CorrWriter *> &writers, const ScalarField &vp,
             const std::vector<ScalarField> &momphases, std::string dsetName)
{
    std::vector<TComplex>   vecBuf;
    std::vector<Complex>    result;
    ScalarField vpPhase(env().getGrid());

    for (unsigned int i_p = 0; i_p < momphases.size(); ++i_p)
    {
        vpPhase = vp*momphases[i_p];
        sliceSum(vpPhase, vecBuf, Tp);
        result.resize(vecBuf.size());
        for (unsigned int t = 0; t < vecBuf.size(); ++t)
        {
            result[t] = TensorRemove(vecBuf[t]);
        }
        if (env().getGrid()->IsBoss())
        {
            write(*writers[i_p], dsetName, result);
        }
    }
}

void TScalarVP::momD1(ScalarField &s, FFT &fft)
{
    EmField     &A = *env().getObject<EmField>(par().emField);
    ScalarField buf(env().getGrid()), result(env().getGrid()),
                Amu(env().getGrid());
    Complex     ci(0.0,1.0);

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
