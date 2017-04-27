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
	for (unsigned int mu = 0; mu < env().getNd(); ++mu)
    {
        phaseName_.push_back("_shiftphase_" + std::to_string(mu));
        muPropQName_.push_back(getName() + "_propQ_" + std::to_string(mu));
        muPropSunName_.push_back(getName() + "_propSun_" + std::to_string(mu));
        muPropTadName_.push_back(getName() + "_propTad_" + std::to_string(mu));
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
    env().registerLattice<ScalarField>(getName());
}

// execution ///////////////////////////////////////////////////////////////////
void TScalarVP::execute(void)
{
	// CACHING ANALYTIC EXPRESSIONS
    ScalarField &source = *env().getObject<ScalarField>(par().source);
    Complex     ci(0.0,1.0);
    FFT         fft(env().getGrid());
    double      q = par().charge;

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
