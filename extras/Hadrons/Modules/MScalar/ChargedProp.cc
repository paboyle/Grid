#include <Grid/Hadrons/Modules/MScalar/ChargedProp.hpp>
#include <Grid/Hadrons/Modules/MScalar/Scalar.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MScalar;

/******************************************************************************
*                     TChargedProp implementation                             *
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
TChargedProp::TChargedProp(const std::string name)
: Module<ChargedPropPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
std::vector<std::string> TChargedProp::getInput(void)
{
    std::vector<std::string> in = {par().source, par().emField};
    
    return in;
}

std::vector<std::string> TChargedProp::getOutput(void)
{
    std::vector<std::string> out = {getName(), getName()+"_Q",
                                    getName()+"_Sun", getName()+"_Tad"};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
void TChargedProp::setup(void)
{
    freeMomPropName_ = FREEMOMPROP(par().mass);
    phaseName_.clear();
    for (unsigned int mu = 0; mu < env().getNd(); ++mu)
    {
        phaseName_.push_back("_shiftphase_" + std::to_string(mu));
    }
    GFSrcName_ = "_" + getName() + "_DinvSrc";
    prop0Name_ = getName() + "_0";
    propQName_ = getName() + "_Q";
    propSunName_ = getName() + "_Sun";
    propTadName_ = getName() + "_Tad";
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
    env().registerLattice<ScalarField>(getName());
    env().registerLattice<ScalarField>(propQName_);
    env().registerLattice<ScalarField>(propSunName_);
    env().registerLattice<ScalarField>(propTadName_);
}

// execution ///////////////////////////////////////////////////////////////////
void TChargedProp::execute(void)
{
    // CACHING ANALYTIC EXPRESSIONS
    ScalarField &source = *env().getObject<ScalarField>(par().source);
    Complex     ci(0.0,1.0);
    FFT         fft(env().getGrid());
    
    // cache momentum-space free scalar propagator
    if (!env().hasCreatedObject(freeMomPropName_))
    {
        LOG(Message) << "Caching momentum space free scalar propagator"
                     << " (mass= " << par().mass << ")..." << std::endl;
        freeMomProp_ = env().createLattice<ScalarField>(freeMomPropName_);
        SIMPL::MomentumSpacePropagator(*freeMomProp_, par().mass);
    }
    else
    {
        freeMomProp_ = env().getObject<ScalarField>(freeMomPropName_);
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
    // cache position-space free scalar propagator
    if (!env().hasCreatedObject(prop0Name_))
    {
        prop0_ = env().createLattice<ScalarField>(prop0Name_);
        *prop0_ = *GFSrc_;
        fft.FFT_all_dim(*prop0_, *prop0_, FFT::backward);
    }
    else
    {
        prop0_ = env().getObject<ScalarField>(prop0Name_);
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

    // PROPAGATOR CALCULATION
    LOG(Message) << "Computing charged scalar propagator"
                 << " (mass= " << par().mass
                 << ", charge= " << par().charge << ")..." << std::endl;
    
    ScalarField &prop   = *env().createLattice<ScalarField>(getName());
    ScalarField &propQ   = *env().createLattice<ScalarField>(propQName_);
    ScalarField &propSun   = *env().createLattice<ScalarField>(propSunName_);
    ScalarField &propTad   = *env().createLattice<ScalarField>(propTadName_);
    ScalarField buf(env().getGrid());
    ScalarField &GFSrc = *GFSrc_, &G = *freeMomProp_;
    double      q = par().charge;
    
    // -G*momD1*G*F*Src (momD1 = F*D1*Finv)
    buf = GFSrc;
    momD1(buf, fft);
    buf = -G*buf;
    fft.FFT_all_dim(propQ, buf, FFT::backward);

    // G*momD1*G*momD1*G*F*Src (here buf = G*momD1*G*F*Src)
    buf = -buf;
    momD1(buf, fft);
    propSun = G*buf;
    fft.FFT_all_dim(propSun, propSun, FFT::backward);

    // -G*momD2*G*F*Src (momD2 = F*D2*Finv)
    buf = GFSrc;
    momD2(buf, fft);
    buf = -G*buf;
    fft.FFT_all_dim(propTad, buf, FFT::backward);

    // full charged scalar propagator
    prop = (*prop0_) + q*propQ + q*q*propSun + q*q*propTad;
    
    // OUTPUT IF NECESSARY
    if (!par().output.empty())
    {
        for (unsigned int i_p = 0; i_p < par().outputMom.size(); ++i_p)
        {
            std::vector<int> mom = strToVec<int>(par().outputMom[i_p]);
            std::string           filename = par().output + "_" + std::to_string(mom[0])
                                                          + std::to_string(mom[1])
                                                          + std::to_string(mom[2])
                                                          + "." +
                                         std::to_string(env().getTrajectory());

            LOG(Message) << "Saving (" << par().outputMom[i_p] << ") momentum projection to '"
                     << filename << "'..." << std::endl;

            CorrWriter            writer(filename);
            std::vector<TComplex> vecBuf;
            std::vector<Complex>  result;

            write(writer, "charge", q);
            write(writer, "mass", par().mass);

            // Write full propagator
            buf = prop;
            for (unsigned int j = 0; j < env().getNd()-1; ++j)
            {
                for (unsigned int momcount = 0; momcount < mom[j]; ++momcount)
                {
                    buf = buf*adj(*phase_[j]);
                }
            }
            sliceSum(buf, vecBuf, Tp);
            result.resize(vecBuf.size());
            for (unsigned int t = 0; t < vecBuf.size(); ++t)
            {
                result[t] = TensorRemove(vecBuf[t]);
            }
            write(writer, "prop", result);

            // Write free propagator
            buf = *prop0_;
            for (unsigned int j = 0; j < env().getNd()-1; ++j)
            {
                for (unsigned int momcount = 0; momcount < mom[j]; ++momcount)
                {
                    buf = buf*adj(*phase_[j]);
                }
            }
            sliceSum(buf, vecBuf, Tp);
            for (unsigned int t = 0; t < vecBuf.size(); ++t)
            {
                result[t] = TensorRemove(vecBuf[t]);
            }
            write(writer, "prop_0", result);

            // Write propagator O(q) term
            buf = propQ;
            for (unsigned int j = 0; j < env().getNd()-1; ++j)
            {
                for (unsigned int momcount = 0; momcount < mom[j]; ++momcount)
                {
                    buf = buf*adj(*phase_[j]);
                }
            }
            sliceSum(buf, vecBuf, Tp);
            for (unsigned int t = 0; t < vecBuf.size(); ++t)
            {
                result[t] = TensorRemove(vecBuf[t]);
            }
            write(writer, "prop_Q", result);

            // Write propagator sunset term
            buf = propSun;
            for (unsigned int j = 0; j < env().getNd()-1; ++j)
            {
                for (unsigned int momcount = 0; momcount < mom[j]; ++momcount)
                {
                    buf = buf*adj(*phase_[j]);
                }
            }
            sliceSum(buf, vecBuf, Tp);
            for (unsigned int t = 0; t < vecBuf.size(); ++t)
            {
                result[t] = TensorRemove(vecBuf[t]);
            }
            write(writer, "prop_Sun", result);

            // Write propagator tadpole term
            buf = propTad;
            for (unsigned int j = 0; j < env().getNd()-1; ++j)
            {
                for (unsigned int momcount = 0; momcount < mom[j]; ++momcount)
                {
                    buf = buf*adj(*phase_[j]);
                }
            }
            sliceSum(buf, vecBuf, Tp);
            for (unsigned int t = 0; t < vecBuf.size(); ++t)
            {
                result[t] = TensorRemove(vecBuf[t]);
            }
            write(writer, "prop_Tad", result);
        }
    }
}

void TChargedProp::momD1(ScalarField &s, FFT &fft)
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

void TChargedProp::momD2(ScalarField &s, FFT &fft)
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
