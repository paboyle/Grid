#include <Grid/Hadrons/Modules/MScalar/ChargedProp.hpp>
#include <Grid/Hadrons/Modules/MScalar/ScalarVP.hpp>
#include <Grid/Hadrons/Modules/MScalar/Scalar.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MScalar;

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
            out.push_back(getName() + "_" + std::to_string(mu) + "_" + std::to_string(nu));
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
    ScalarField Amu(env().getGrid()), tmp_vp(env().getGrid());
    TComplex    Anu0;
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
    std::vector<TComplex>   vecBuf;
    std::vector<Complex>    result;
    ScalarField vpPhase(env().getGrid());
    std::vector<CorrWriter *> writer, writer0, writerD;
    std::vector<ScalarField> momphases;
    if (!par().output.empty())
    {
        LOG(Message) << "Preparing output files..." << std::endl;
        for (unsigned int i_p = 0; i_p < par().outputMom.size(); ++i_p)
        {
            std::vector<int> mom = strToVec<int>(par().outputMom[i_p]);

            // Open output files
            std::string           filename = par().output + "_" + std::to_string(mom[0])
                                                          + std::to_string(mom[1])
                                                          + std::to_string(mom[2])
                                                          + "." +
                                             std::to_string(env().getTrajectory());
            std::string           filename0 = par().output + "_" + std::to_string(mom[0])
                                                           + std::to_string(mom[1])
                                                           + std::to_string(mom[2])
                                                           + "_free." +
                                              std::to_string(env().getTrajectory());
            std::string           filenameD = par().output + "_" + std::to_string(mom[0])
                                                           + std::to_string(mom[1])
                                                           + std::to_string(mom[2])
                                                           + "_diagrams." +
                                              std::to_string(env().getTrajectory());

            if (env().getGrid()->IsBoss())
            {
                CorrWriter *writer_i = new CorrWriter(filename);
                writer.push_back(writer_i);
                CorrWriter *writer0_i = new CorrWriter(filename0);
                writer0.push_back(writer0_i);
                CorrWriter *writerD_i = new CorrWriter(filenameD);
                writerD.push_back(writerD_i);

                write(*writer[i_p], "charge", q);
                write(*writer[i_p], "mass", static_cast<TChargedProp *>(env().getModule(par().scalarProp))->par().mass);
                write(*writer0[i_p], "charge", 0.0);
                write(*writer0[i_p], "mass", static_cast<TChargedProp *>(env().getModule(par().scalarProp))->par().mass);
                write(*writerD[i_p], "charge", q);
                write(*writerD[i_p], "mass", static_cast<TChargedProp *>(env().getModule(par().scalarProp))->par().mass);
            }

            // Calculate phase factors
            vpPhase = Complex(1.0,0.0);
            for (unsigned int j = 0; j < env().getNd()-1; ++j)
            {
                for (unsigned int momcount = 0; momcount < mom[j]; ++momcount)
                {
                    vpPhase = vpPhase*(*phase_[j]);
                }
            }
            vpPhase = adj(vpPhase);
            momphases.push_back(vpPhase);
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

            // Free VP
            prop1 = *prop0_;
            prop2 = Cshift(*prop0_, nu, -1);
            tmp_vp = adj(prop2) * Cshift(prop1, mu, 1);
            tmp_vp -= Cshift(adj(prop2), mu, 1) * prop1;
            tmp_vp = 2.0*real(tmp_vp);
            vpTensor[mu][nu] = tmp_vp;

            // Output if necessary
            if (!par().output.empty())
            {
                for (unsigned int i_p = 0; i_p < par().outputMom.size(); ++i_p)
                {
                    vpPhase = tmp_vp*momphases[i_p];
                    sliceSum(vpPhase, vecBuf, Tp);
                    result.resize(vecBuf.size());
                    for (unsigned int t = 0; t < vecBuf.size(); ++t)
                    {
                        result[t] = TensorRemove(vecBuf[t]);
                    }
                    if (env().getGrid()->IsBoss())
                    {
                        write(*writer0[i_p],
                              "Pi_"+std::to_string(mu)+"_"+std::to_string(nu),
                              result);
                    }
                }
            }

            // S
            // X
            // 4C

            // "Exchange" terms
            prop1 += q*propQ;
            prop2 += q*muPropQ[nu];
            tmp_vp = adj(prop2) * (1.0 + ci*q*Amu)
                     * Cshift(prop1, mu, 1) * (1.0 + ci*q*Anu0);
            tmp_vp -= Cshift(adj(prop2), mu, 1) * (1.0 - ci*q*Amu)
                      * prop1 * (1.0 + ci*q*Anu0);
            tmp_vp = 2.0*real(tmp_vp);
            vpTensor[mu][nu] = tmp_vp;

            // Output if necessary
            if (!par().output.empty())
            {
                for (unsigned int i_p = 0; i_p < par().outputMom.size(); ++i_p)
                {
                    vpPhase = tmp_vp*momphases[i_p];
                    sliceSum(vpPhase, vecBuf, Tp);
                    result.resize(vecBuf.size());
                    for (unsigned int t = 0; t < vecBuf.size(); ++t)
                    {
                        result[t] = TensorRemove(vecBuf[t]);
                    }
                    if (env().getGrid()->IsBoss())
                    {
                        write(*writerD[i_p],
                              "Pi_exchange_"+std::to_string(mu)+"_"+std::to_string(nu),
                              result);
                    }
                }
            }

            // Subtract O(alpha^2) term
            prop1 = q*propQ;
            prop2 = q*muPropQ[nu];
            tmp_vp = Cshift(adj(prop2), mu, 1) * (-ci)*q*Amu
                     * prop1 * ci*q*Anu0;
            tmp_vp -= adj(prop2) * ci*q*Amu
                      * Cshift(prop1, mu, 1) * ci*q*Anu0;
            tmp_vp = 2.0*real(tmp_vp);
            vpTensor[mu][nu] += tmp_vp;

            // Output if necessary
            if (!par().output.empty())
            {
                for (unsigned int i_p = 0; i_p < par().outputMom.size(); ++i_p)
                {
                    vpPhase = tmp_vp*momphases[i_p];
                    sliceSum(vpPhase, vecBuf, Tp);
                    result.resize(vecBuf.size());
                    for (unsigned int t = 0; t < vecBuf.size(); ++t)
                    {
                        result[t] = TensorRemove(vecBuf[t]);
                    }
                    if (env().getGrid()->IsBoss())
                    {
                        write(*writerD[i_p],
                              "Pi_alpha2_"+std::to_string(mu)+"_"+std::to_string(nu),
                              result);
                    }
                }
            }

            // Sunset from unshifted source
            prop1 = q*q*propSun;
            prop2 = Cshift(*prop0_, nu, -1);
            tmp_vp = adj(prop2) * Cshift(prop1, mu, 1);
            tmp_vp -= Cshift(adj(prop2), mu, 1) * prop1;
            tmp_vp = 2.0*real(tmp_vp);
            vpTensor[mu][nu] += tmp_vp;

            // Output if necessary
            if (!par().output.empty())
            {
                for (unsigned int i_p = 0; i_p < par().outputMom.size(); ++i_p)
                {
                    vpPhase = tmp_vp*momphases[i_p];
                    sliceSum(vpPhase, vecBuf, Tp);
                    result.resize(vecBuf.size());
                    for (unsigned int t = 0; t < vecBuf.size(); ++t)
                    {
                        result[t] = TensorRemove(vecBuf[t]);
                    }
                    if (env().getGrid()->IsBoss())
                    {
                        write(*writerD[i_p],
                              "Pi_sunset_unshifted_"+std::to_string(mu)+"_"+std::to_string(nu),
                              result);
                    }
                }
            }

            // Sunset from shifted source
            prop1 = Cshift(prop1, nu, -1);
            tmp_vp = Cshift(adj(*prop0_), mu, 1) * prop1;
            tmp_vp -= adj(*prop0_) * Cshift(prop1, mu, 1);
            tmp_vp = 2.0*real(tmp_vp);
            vpTensor[mu][nu] += tmp_vp;

            // Output if necessary
            if (!par().output.empty())
            {
                for (unsigned int i_p = 0; i_p < par().outputMom.size(); ++i_p)
                {
                    vpPhase = tmp_vp*momphases[i_p];
                    sliceSum(vpPhase, vecBuf, Tp);
                    result.resize(vecBuf.size());
                    for (unsigned int t = 0; t < vecBuf.size(); ++t)
                    {
                        result[t] = TensorRemove(vecBuf[t]);
                    }
                    if (env().getGrid()->IsBoss())
                    {
                        write(*writerD[i_p],
                              "Pi_sunset_shifted_"+std::to_string(mu)+"_"+std::to_string(nu),
                              result);
                    }
                }
            }

            // Tadpole from unshifted source
            prop1 = q*q*propTad;
            prop2 = Cshift(*prop0_, nu, -1);
            tmp_vp = adj(prop2) * Cshift(prop1, mu, 1);
            tmp_vp -= Cshift(adj(prop2), mu, 1) * prop1;
            tmp_vp = 2.0*real(tmp_vp);
            vpTensor[mu][nu] += tmp_vp;

            // Output if necessary
            if (!par().output.empty())
            {
                for (unsigned int i_p = 0; i_p < par().outputMom.size(); ++i_p)
                {
                    vpPhase = tmp_vp*momphases[i_p];
                    sliceSum(vpPhase, vecBuf, Tp);
                    result.resize(vecBuf.size());
                    for (unsigned int t = 0; t < vecBuf.size(); ++t)
                    {
                        result[t] = TensorRemove(vecBuf[t]);
                    }
                    if (env().getGrid()->IsBoss())
                    {
                        write(*writerD[i_p],
                              "Pi_tadpole_unshifted_"+std::to_string(mu)+"_"+std::to_string(nu),
                              result);
                    }
                }
            }

            // Tadpole from shifted source
            prop1 = Cshift(prop1, nu, -1);
            tmp_vp = Cshift(adj(*prop0_), mu, 1) * prop1;
            tmp_vp -= adj(*prop0_) * Cshift(prop1, mu, 1);
            tmp_vp = 2.0*real(tmp_vp);
            vpTensor[mu][nu] += tmp_vp;

            // Output if necessary
            if (!par().output.empty())
            {
                for (unsigned int i_p = 0; i_p < par().outputMom.size(); ++i_p)
                {
                    vpPhase = tmp_vp*momphases[i_p];
                    sliceSum(vpPhase, vecBuf, Tp);
                    result.resize(vecBuf.size());
                    for (unsigned int t = 0; t < vecBuf.size(); ++t)
                    {
                        result[t] = TensorRemove(vecBuf[t]);
                    }
                    if (env().getGrid()->IsBoss())
                    {
                        write(*writerD[i_p],
                              "Pi_tadpole_shifted_"+std::to_string(mu)+"_"+std::to_string(nu),
                              result);
                    }
                }
            }

            // Source tadpole
            prop1 = *prop0_;
            tmp_vp = adj(prop2)
                     * Cshift(prop1, mu, 1)
                     * (-0.5)*q*q*Anu0*Anu0;
            tmp_vp -= Cshift(adj(prop2), mu, 1)
                      * prop1
                      * (-0.5)*q*q*Anu0*Anu0;
            tmp_vp = 2.0*real(tmp_vp);
            vpTensor[mu][nu] += tmp_vp;

            // Output if necessary
            if (!par().output.empty())
            {
                for (unsigned int i_p = 0; i_p < par().outputMom.size(); ++i_p)
                {
                    vpPhase = tmp_vp*momphases[i_p];
                    sliceSum(vpPhase, vecBuf, Tp);
                    result.resize(vecBuf.size());
                    for (unsigned int t = 0; t < vecBuf.size(); ++t)
                    {
                        result[t] = TensorRemove(vecBuf[t]);
                    }
                    if (env().getGrid()->IsBoss())
                    {
                        write(*writerD[i_p],
                              "Pi_sourcetadpole_"+std::to_string(mu)+"_"+std::to_string(nu),
                              result);
                    }
                }
            }

            // Sink tadpole
            tmp_vp = adj(prop2)
                     * (-0.5)*q*q*Amu*Amu
                     * Cshift(prop1, mu, 1);
            tmp_vp -= Cshift(adj(prop2), mu, 1)
                      * (-0.5)*q*q*Amu*Amu
                      * prop1;
            tmp_vp = 2.0*real(tmp_vp);
            vpTensor[mu][nu] += tmp_vp;

            // Output if necessary
            if (!par().output.empty())
            {
                for (unsigned int i_p = 0; i_p < par().outputMom.size(); ++i_p)
                {
                    vpPhase = tmp_vp*momphases[i_p];
                    sliceSum(vpPhase, vecBuf, Tp);
                    result.resize(vecBuf.size());
                    for (unsigned int t = 0; t < vecBuf.size(); ++t)
                    {
                        result[t] = TensorRemove(vecBuf[t]);
                    }
                    if (env().getGrid()->IsBoss())
                    {
                        write(*writerD[i_p],
                              "Pi_sinktadpole_"+std::to_string(mu)+"_"+std::to_string(nu),
                              result);
                    }
                }
            }

            // Output if necessary
            if (!par().output.empty())
            {
                for (unsigned int i_p = 0; i_p < par().outputMom.size(); ++i_p)
                {
                    vpPhase = vpTensor[mu][nu]*momphases[i_p];
                    sliceSum(vpPhase, vecBuf, Tp);
                    result.resize(vecBuf.size());
                    for (unsigned int t = 0; t < vecBuf.size(); ++t)
                    {
                        result[t] = TensorRemove(vecBuf[t]);
                    }
                    if (env().getGrid()->IsBoss())
                    {
                        write(*writer[i_p],
                              "Pi_"+std::to_string(mu)+"_"+std::to_string(nu),
                              result);
                    }
                }
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
                delete writer0[i_p];
                delete writerD[i_p];
            }
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
