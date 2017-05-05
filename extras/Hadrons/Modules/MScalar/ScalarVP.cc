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
    std::vector<std::string> out = {getName(), getName()+"_propQ",
                                    getName()+"_propSun",
                                    getName()+"_propTad"};
    
    for (unsigned int mu = 0; mu < env().getNd(); ++mu)
    {
        out.push_back(getName() + "_propQ_" + std::to_string(mu));
        out.push_back(getName() + "_propSun_" + std::to_string(mu));
        out.push_back(getName() + "_propTad_" + std::to_string(mu));

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
	freeMomPropName_ = FREEMOMPROP(par().mass);
	GFSrcName_ = "_" + getName() + "_DinvSrc";
    prop0Name_ = getName() + "_prop0";
    propQName_ = getName() + "_propQ";
    propSunName_ = getName() + "_propSun";
    propTadName_ = getName() + "_propTad";

	phaseName_.clear();
	muPropQName_.clear();
	muPropSunName_.clear();
	muPropTadName_.clear();
    vpTensorName_.clear();
    
	for (unsigned int mu = 0; mu < env().getNd(); ++mu)
    {
        phaseName_.push_back("_shiftphase_" + std::to_string(mu));
        muPropQName_.push_back(getName() + "_propQ_" + std::to_string(mu));
        muPropSunName_.push_back(getName() + "_propSun_" + std::to_string(mu));
        muPropTadName_.push_back(getName() + "_propTad_" + std::to_string(mu));

        std::vector<std::string> vpTensorName_mu;
        for (unsigned int nu = 0; nu < env().getNd(); ++nu)
        {
            vpTensorName_mu.push_back(getName() + "_" + std::to_string(mu)
                                      + "_" + std::to_string(nu));
        }
        vpTensorName_.push_back(vpTensorName_mu);
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
    for (unsigned int mu = 0; mu < env().getNd(); ++mu)
    {
		env().registerLattice<ScalarField>(muPropSunName_[mu]);
	}
    env().registerLattice<ScalarField>(propTadName_);
    for (unsigned int mu = 0; mu < env().getNd(); ++mu)
    {
		env().registerLattice<ScalarField>(muPropTadName_[mu]);
	}
    for (unsigned int mu = 0; mu < env().getNd(); ++mu)
    {
        for (unsigned int nu = 0; nu < env().getNd(); ++nu)
        {
            env().registerLattice<ScalarField>(vpTensorName_[mu][nu]);
        }
    }
    env().registerLattice<ScalarField>(getName());
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
	ScalarField &propQ   = *env().createLattice<ScalarField>(propQName_);
    ScalarField &propSun   = *env().createLattice<ScalarField>(propSunName_);
    ScalarField &propTad   = *env().createLattice<ScalarField>(propTadName_);
    chargedProp(propQ, propSun, propTad, *GFSrc_, fft);

    // Propagators from shifted sources
    std::vector<ScalarField *> muPropQ_, muPropSun_, muPropTad_;
    ScalarField buf(env().getGrid());
    for (unsigned int mu = 0; mu < env().getNd(); ++mu)
    {
        muPropQ_.push_back(env().createLattice<ScalarField>(muPropQName_[mu]));
        muPropSun_.push_back(env().createLattice<ScalarField>(muPropSunName_[mu]));
        muPropTad_.push_back(env().createLattice<ScalarField>(muPropTadName_[mu]));

        buf = adj(*phase_[mu])*(*GFSrc_);
        chargedProp(*(muPropQ_[mu]), *(muPropSun_[mu]), *(muPropTad_[mu]),
        			buf, fft);
    }

    // CONTRACTIONS
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
    ScalarField prop1(env().getGrid()), prop2(env().getGrid());
    EmField     &A = *env().getObject<EmField>(par().emField);
    ScalarField Amu(env().getGrid());
    TComplex    Anu0;
    std::vector<int> coor0 = {0, 0, 0, 0};

    prop1 = *GFSrc_ + q*propQ + q*q*propSun + q*q*propTad;
    fft.FFT_all_dim(prop1, prop1, FFT::backward);
    for (unsigned int nu = 0; nu < env().getNd(); ++nu)
    {
        peekSite(Anu0, peekLorentz(A, nu), coor0);
        prop2 = adj(*phase_[nu])*(*GFSrc_) + q*(*(muPropQ_[nu]))
                + q*q*(*(muPropSun_[nu]) + *(muPropTad_[nu]));
        fft.FFT_all_dim(prop2, prop2, FFT::backward);

        std::vector<ScalarField> pi_nu;
        for (unsigned int mu = 0; mu < env().getNd(); ++mu)
        {
            LOG(Message) << "Computing Pi[" << mu << "][" << nu << "]..."
                         << std::endl;
            Amu = peekLorentz(A, mu);
            vpTensor[mu][nu] = adj(prop2)
                               * (1.0 + ci*q*Amu - 0.5*q*q*Amu*Amu)
                               * Cshift(prop1, mu, 1)
                               * (1.0 + ci*q*Anu0 - 0.5*q*q*Anu0*Anu0);
            vpTensor[mu][nu] -= Cshift(adj(prop2), mu, 1)
                                * (1.0 - ci*q*Amu - 0.5*q*q*Amu*Amu)
                                * prop1
                                * (1.0 + ci*q*Anu0 - 0.5*q*q*Anu0*Anu0);
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
            }
        }
    }
}

// Calculate O(q) and O(q^2) terms of momentum-space charged propagator
void TScalarVP::chargedProp(ScalarField &prop_q, ScalarField &prop_sun,
							ScalarField &prop_tad, ScalarField &GFSrc,
							FFT &fft)
{
	Complex     ci(0.0,1.0);
	ScalarField &G = *freeMomProp_;
    ScalarField buf(env().getGrid());

    LOG(Message) << "Computing charged scalar propagator"
                 << " (mass= " << par().mass
                 << ", charge= " << par().charge << ")..."
                 << std::endl;
	
	// -G*momD1*G*F*Src (momD1 = F*D1*Finv)
    buf = GFSrc;
    momD1(buf, fft);
    buf = G*buf;
    prop_q = -buf;

    // G*momD1*G*momD1*G*F*Src
    momD1(buf, fft);
    prop_sun = G*buf;

    // -G*momD2*G*F*Src (momD2 = F*D2*Finv)
    buf = GFSrc;
    momD2(buf, fft);
    prop_tad = -G*buf;
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
