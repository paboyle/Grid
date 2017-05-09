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
	std::vector<std::string> in = {par().source, par().emField};
    
    return in;
}

std::vector<std::string> TScalarVP::getOutput(void)
{
    std::vector<std::string> out = {getName()+"_propQ",
                                    getName()+"_propSun",
                                    getName()+"_propTad"};
    
    for (unsigned int mu = 0; mu < env().getNd(); ++mu)
    {
        out.push_back(getName() + "_propQ_" + std::to_string(mu));

        for (unsigned int nu = 0; nu < env().getNd(); ++nu)
        {
            out.push_back(getName() + "_" + std::to_string(mu) + "_" + std::to_string(nu));
            out.push_back(getName() + "_free_" + std::to_string(mu) + "_" + std::to_string(nu));
        }
    }

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
void TScalarVP::setup(void)
{
	freeMomPropName_ = FREEMOMPROP(par().mass);
	GFSrcName_ = "_" + getName() + "_DinvSrc";
    prop0Name_ = getName() + "_prop0";
    propQName_ = getName() + "_propQ";
    propSunName_ = getName() + "_propSun";
    propTadName_ = getName() + "_propTad";

	phaseName_.clear();
	muPropQName_.clear();
    vpTensorName_.clear();
    freeVpTensorName_.clear();
    
	for (unsigned int mu = 0; mu < env().getNd(); ++mu)
    {
        phaseName_.push_back("_shiftphase_" + std::to_string(mu));
        muPropQName_.push_back(getName() + "_propQ_" + std::to_string(mu));

        std::vector<std::string> vpTensorName_mu;
        std::vector<std::string> freeVpTensorName_mu;
        for (unsigned int nu = 0; nu < env().getNd(); ++nu)
        {
            vpTensorName_mu.push_back(getName() + "_" + std::to_string(mu)
                                      + "_" + std::to_string(nu));
            freeVpTensorName_mu.push_back(getName() + "_free_" + std::to_string(mu)
                                       + "_" + std::to_string(nu));
        }
        vpTensorName_.push_back(vpTensorName_mu);
        freeVpTensorName_.push_back(freeVpTensorName_mu);
    }

    if (!env().hasRegisteredObject(freeMomPropName_))
    {
        env().registerLattice<ScalarField>(freeMomPropName_);
    }
    if (!env().hasRegisteredObject(phaseName_[0]))
    {
        for (unsigned int mu = 0; mu < env().getNd(); ++mu)
        {
            env().registerLattice<ScalarField>(phaseName_[mu]);
        }
    }
    if (!env().hasRegisteredObject(GFSrcName_))
    {
        env().registerLattice<ScalarField>(GFSrcName_);
    }
    if (!env().hasRegisteredObject(prop0Name_))
    {
        env().registerLattice<ScalarField>(prop0Name_);
    }
    env().registerLattice<ScalarField>(propQName_);
    for (unsigned int mu = 0; mu < env().getNd(); ++mu)
	{
	    env().registerLattice<ScalarField>(muPropQName_[mu]);
	}
    env().registerLattice<ScalarField>(propSunName_);
    env().registerLattice<ScalarField>(propTadName_);
    for (unsigned int mu = 0; mu < env().getNd(); ++mu)
    {
        for (unsigned int nu = 0; nu < env().getNd(); ++nu)
        {
            env().registerLattice<ScalarField>(vpTensorName_[mu][nu]);
            env().registerLattice<ScalarField>(freeVpTensorName_[mu][nu]);
        }
    }
}

// execution ///////////////////////////////////////////////////////////////////
void TScalarVP::execute(void)
{
	// CACHING ANALYTIC EXPRESSIONS
    ScalarField &source = *env().getObject<ScalarField>(par().source);
    Complex     ci(0.0,1.0);
    FFT         fft(env().getGrid());
    Real        q = par().charge;

    // cache momentum-space free scalar propagator
    if (!env().hasCreatedObject(freeMomPropName_))
    {
        LOG(Message) << "Caching momentum space free scalar propagator"
                     << " (mass= " << par().mass << ")..." << std::endl;
        freeMomProp_ = env().createLattice<ScalarField>(freeMomPropName_);
        Scalar<SIMPL>::MomentumSpacePropagator(*freeMomProp_, par().mass);
    }
    else
    {
        freeMomProp_ = env().getObject<ScalarField>(freeMomPropName_);
    }
    // cache phases
    if (!env().hasCreatedObject(phaseName_[0]))
    {
        std::vector<int> &l = env().getGrid()->_fdimensions;
        
        LOG(Message) << "Caching shift phases..." << std::endl;
        for (unsigned int mu = 0; mu < env().getNd(); ++mu)
        {
            Real    twoPiL = M_PI*2./l[mu];
            
            phase_.push_back(env().createLattice<ScalarField>(phaseName_[mu]));
            LatticeCoordinate(*(phase_[mu]), mu);
            *(phase_[mu]) = exp(ci*twoPiL*(*(phase_[mu])));
        }
    }
    else
    {
        for (unsigned int mu = 0; mu < env().getNd(); ++mu)
        {
            phase_.push_back(env().getObject<ScalarField>(phaseName_[mu]));
        }
    }
    // cache G*F*src
    if (!env().hasCreatedObject(GFSrcName_))
    {
        GFSrc_ = env().createLattice<ScalarField>(GFSrcName_);
        fft.FFT_all_dim(*GFSrc_, source, FFT::forward);
        *GFSrc_ = (*freeMomProp_)*(*GFSrc_);
    }
    else
    {
        GFSrc_ = env().getObject<ScalarField>(GFSrcName_);
    }
    // cache position-space free scalar propagators
    if (!env().hasCreatedObject(prop0Name_))
    {
        prop0_ = env().createLattice<ScalarField>(prop0Name_);
        fft.FFT_all_dim(*prop0_, *GFSrc_, FFT::backward);
    }
    else
    {
        prop0_ = env().getObject<ScalarField>(prop0Name_);
    }

    // PROPAGATOR CALCULATION
    // Propagator from unshifted source
    LOG(Message) << "Computing O(alpha) charged scalar propagator"
                 << " (mass= " << par().mass
                 << ", charge= " << q << ")..."
                 << std::endl;
	ScalarField &propQ   = *env().createLattice<ScalarField>(propQName_);
    ScalarField &propSun   = *env().createLattice<ScalarField>(propSunName_);
    ScalarField &propTad   = *env().createLattice<ScalarField>(propTadName_);
    chargedProp(propQ, propSun, propTad, *GFSrc_, fft);

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
    ScalarField Amu(env().getGrid());
    TComplex    Anu0;
    std::vector<int> coor0 = {0, 0, 0, 0};
    std::vector<std::vector<ScalarField> > vpTensor, freeVpTensor;
    for (unsigned int mu = 0; mu < env().getNd(); ++mu)
    {
        std::vector<ScalarField> vpTensor_mu;
        std::vector<ScalarField> freeVpTensor_mu;
        for (unsigned int nu = 0; nu < env().getNd(); ++nu)
        {
            vpTensor_mu.push_back(*env().createLattice<ScalarField>(vpTensorName_[mu][nu]));
            freeVpTensor_mu.push_back(*env().createLattice<ScalarField>(freeVpTensorName_[mu][nu]));
        }
        vpTensor.push_back(vpTensor_mu);
        freeVpTensor.push_back(freeVpTensor_mu);
    }

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
            freeVpTensor[mu][nu] = adj(prop2) * Cshift(prop1, mu, 1);
            freeVpTensor[mu][nu] -= Cshift(adj(prop2), mu, 1) * prop1;

            // "Exchange" terms
            prop1 += q*propQ;
            prop2 += q*muPropQ[nu];
            vpTensor[mu][nu] = adj(prop2) * (1.0 + ci*q*Amu)
                               * Cshift(prop1, mu, 1) * (1.0 + ci*q*Anu0);
            vpTensor[mu][nu] -= Cshift(adj(prop2), mu, 1) * (1.0 - ci*q*Amu)
                                * prop1 * (1.0 + ci*q*Anu0);

            // Subtract O(alpha^2) term
            prop1 = q*propQ;
            prop2 = q*muPropQ[nu];
            vpTensor[mu][nu] -= adj(prop2) * ci*q*Amu
                                * Cshift(prop1, mu, 1) * ci*q*Anu0;
            vpTensor[mu][nu] += Cshift(adj(prop2), mu, 1) * (-ci)*q*Amu
                                * prop1 * ci*q*Anu0;

            // Sunset+tadpole from source
            prop1 = q*q*(propSun + propTad);
            prop2 = Cshift(*prop0_, nu, -1);
            vpTensor[mu][nu] += adj(prop2) * Cshift(prop1, mu, 1);
            vpTensor[mu][nu] -= Cshift(adj(prop2), mu, 1) * prop1;

            // Sunset+tadpole from shifted source
            prop1 = Cshift(prop1, nu, -1);
            vpTensor[mu][nu] += Cshift(adj(*prop0_), mu, 1) * prop1;
            vpTensor[mu][nu] -= adj(*prop0_) * Cshift(prop1, mu, 1);

            // Source tadpole
            prop1 = *prop0_;
            vpTensor[mu][nu] += adj(prop2)
                                * Cshift(prop1, mu, 1)
                                * (-0.5)*q*q*Anu0*Anu0;
            vpTensor[mu][nu] -= Cshift(adj(prop2), mu, 1)
                                * prop1
                                * (-0.5)*q*q*Anu0*Anu0;

            // Sink tadpole
            vpTensor[mu][nu] += adj(prop2)
                                * (-0.5)*q*q*Amu*Amu
                                * Cshift(prop1, mu, 1);
            vpTensor[mu][nu] -= Cshift(adj(prop2), mu, 1)
                                * (-0.5)*q*q*Amu*Amu
                                * prop1;

            vpTensor[mu][nu] = 2.0*real(vpTensor[mu][nu]);
        }
    }

    // OUTPUT IF NECESSARY
    if (!par().output.empty())
    {
        std::string           filename = par().output + "." +
                                         std::to_string(env().getTrajectory());
        
        LOG(Message) << "Saving zero-momentum projection to '"
                     << filename << "'..." << std::endl;
        
        CorrWriter              writer(filename);
        std::vector<TComplex>   vecBuf;
        std::vector<Complex>    result;
        
        write(writer, "charge", q);
        write(writer, "mass", par().mass);

        for (unsigned int mu = 0; mu < env().getNd(); ++mu)
        {
            for (unsigned int nu = 0; nu < env().getNd(); ++nu)
            {
                sliceSum(vpTensor[mu][nu], vecBuf, Tp);
                result.resize(vecBuf.size());
                for (unsigned int t = 0; t < vecBuf.size(); ++t)
                {
                    result[t] = TensorRemove(vecBuf[t]);
                }
                write(writer, "Pi_"+std::to_string(mu)+"_"+std::to_string(nu),
                      result);

                sliceSum(freeVpTensor[mu][nu], vecBuf, Tp);
                result.resize(vecBuf.size());
                for (unsigned int t = 0; t < vecBuf.size(); ++t)
                {
                    result[t] = TensorRemove(vecBuf[t]);
                }
                write(writer,
                      "Pi_"+std::to_string(mu)+"_"+std::to_string(nu)+"_free",
                      result);
            }
        }
    }
}

// Calculate O(q) and O(q^2) terms of position-space charged propagator
void TScalarVP::chargedProp(ScalarField &prop_q, ScalarField &prop_sun,
							ScalarField &prop_tad, ScalarField &GFSrc,
							FFT &fft)
{
	Complex     ci(0.0,1.0);
	ScalarField &G = *freeMomProp_;
    ScalarField buf(env().getGrid());
	
	// -G*momD1*G*F*Src (momD1 = F*D1*Finv)
    buf = GFSrc;
    momD1(buf, fft);
    buf = G*buf;
    prop_q = -buf;
    fft.FFT_all_dim(prop_q, prop_q, FFT::backward);

    // G*momD1*G*momD1*G*F*Src
    momD1(buf, fft);
    prop_sun = G*buf;
    fft.FFT_all_dim(prop_sun, prop_sun, FFT::backward);

    // -G*momD2*G*F*Src (momD2 = F*D2*Finv)
    buf = GFSrc;
    momD2(buf, fft);
    prop_tad = -G*buf;
    fft.FFT_all_dim(prop_tad, prop_tad, FFT::backward);
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

void TScalarVP::momD2(ScalarField &s, FFT &fft)
{
    EmField     &A = *env().getObject<EmField>(par().emField);
    ScalarField buf(env().getGrid()), result(env().getGrid()),
                Amu(env().getGrid());

    result = zero;
    
    for (unsigned int mu = 0; mu < env().getNd(); ++mu)
    {
        Amu = peekLorentz(A, mu);
        buf = (*phase_[mu])*s;
        fft.FFT_all_dim(buf, buf, FFT::backward);
        buf = Amu*Amu*buf;
        fft.FFT_all_dim(buf, buf, FFT::forward);
        result = result + .5*buf;
    }
    fft.FFT_all_dim(s, s, FFT::backward);
    for (unsigned int mu = 0; mu < env().getNd(); ++mu)
    {
        Amu = peekLorentz(A, mu);        
        buf = Amu*Amu*s;
        fft.FFT_all_dim(buf, buf, FFT::forward);
        result = result + .5*adj(*phase_[mu])*buf;
    }

    s = result;
}
