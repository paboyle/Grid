/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MScalar/ChargedProp.cc

Copyright (C) 2015-2018

Author: Antonin Portelli <antonin.portelli@me.com>
Author: James Harrison <jch1g10@soton.ac.uk>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

See the full license in the file "LICENSE" in the top level distribution directory
*************************************************************************************/
/*  END LEGAL */
#include <Hadrons/Modules/MScalar/ChargedProp.hpp>
#include <Hadrons/Modules/MScalar/Scalar.hpp>

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
    std::vector<std::string> out = {getName(), getName()+"_0", getName()+"_Q",
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
    GFSrcName_ = getName() + "_DinvSrc";
	prop0Name_ = getName() + "_0";
    propQName_ = getName() + "_Q";
    propSunName_ = getName() + "_Sun";
    propTadName_ = getName() + "_Tad";
    fftName_   = getName() + "_fft";

    freeMomPropDone_ = env().hasCreatedObject(freeMomPropName_);
    GFSrcDone_       = env().hasCreatedObject(GFSrcName_);
    phasesDone_      = env().hasCreatedObject(phaseName_[0]);
	prop0Done_		 = env().hasCreatedObject(prop0Name_);
    envCacheLat(ScalarField, freeMomPropName_);
    for (unsigned int mu = 0; mu < env().getNd(); ++mu)
    {
        envCacheLat(ScalarField, phaseName_[mu]);
    }
    envCacheLat(ScalarField, GFSrcName_);
	envCacheLat(ScalarField, prop0Name_);
    envCreateLat(ScalarField, getName());
    envCreateLat(ScalarField, propQName_);
    envCreateLat(ScalarField, propSunName_);
    envCreateLat(ScalarField, propTadName_);
    envTmpLat(ScalarField, "buf");
    envTmpLat(ScalarField, "result");
    envTmpLat(ScalarField, "Amu");
    envCache(FFT, fftName_, 1, env().getGrid());
}

// execution ///////////////////////////////////////////////////////////////////
void TChargedProp::execute(void)
{
    // CACHING ANALYTIC EXPRESSIONS
    makeCaches();

    // PROPAGATOR CALCULATION
    LOG(Message) << "Computing charged scalar propagator"
                 << " (mass= " << par().mass
                 << ", charge= " << par().charge << ")..." << std::endl;
    
    auto   &prop    = envGet(ScalarField, getName());
	auto   &prop0   = envGet(ScalarField, prop0Name_);
	auto   &propQ   = envGet(ScalarField, propQName_);
	auto   &propSun = envGet(ScalarField, propSunName_);
	auto   &propTad = envGet(ScalarField, propTadName_);
    auto   &GFSrc   = envGet(ScalarField, GFSrcName_);
    auto   &G       = envGet(ScalarField, freeMomPropName_);
    auto   &fft     = envGet(FFT, fftName_);
    double q        = par().charge;
    envGetTmp(ScalarField, buf);

    // -G*momD1*G*F*Src (momD1 = F*D1*Finv)
    propQ = GFSrc;
    momD1(propQ, fft);
    propQ = -G*propQ;
    propSun = -propQ;
    fft.FFT_dim(propQ, propQ, env().getNd()-1, FFT::backward);

    // G*momD1*G*momD1*G*F*Src (here buf = G*momD1*G*F*Src)
    momD1(propSun, fft);
    propSun = G*propSun;
    fft.FFT_dim(propSun, propSun, env().getNd()-1, FFT::backward);

    // -G*momD2*G*F*Src (momD2 = F*D2*Finv)
    propTad = GFSrc;
    momD2(propTad, fft);
    propTad = -G*propTad;
    fft.FFT_dim(propTad, propTad, env().getNd()-1, FFT::backward);
    
    // full charged scalar propagator
    fft.FFT_dim(buf, GFSrc, env().getNd()-1, FFT::backward);
    prop = buf + q*propQ + q*q*propSun + q*q*propTad;

    // OUTPUT IF NECESSARY
    if (!par().output.empty())
    {
        Result result;
        TComplex            site;
        std::vector<int>    siteCoor;

        LOG(Message) << "Saving momentum-projected propagator to '"
                     << RESULT_FILE_NAME(par().output, vm().getTrajectory()) << "'..."
                     << std::endl;
        result.projection.resize(par().outputMom.size());
        result.lattice_size = env().getGrid()->_fdimensions;
        result.mass = par().mass;
        result.charge = q;
        siteCoor.resize(env().getNd());
        for (unsigned int i_p = 0; i_p < par().outputMom.size(); ++i_p)
        {
            result.projection[i_p].momentum = strToVec<int>(par().outputMom[i_p]);

            LOG(Message) << "Calculating (" << par().outputMom[i_p]
                         << ") momentum projection" << std::endl;

            result.projection[i_p].corr_0.resize(env().getGrid()->_fdimensions[env().getNd()-1]);
            result.projection[i_p].corr.resize(env().getGrid()->_fdimensions[env().getNd()-1]);
            result.projection[i_p].corr_Q.resize(env().getGrid()->_fdimensions[env().getNd()-1]);
            result.projection[i_p].corr_Sun.resize(env().getGrid()->_fdimensions[env().getNd()-1]);
            result.projection[i_p].corr_Tad.resize(env().getGrid()->_fdimensions[env().getNd()-1]);

            for (unsigned int j = 0; j < env().getNd()-1; ++j)
            {
                siteCoor[j] = result.projection[i_p].momentum[j];
            }

            for (unsigned int t = 0; t < result.projection[i_p].corr.size(); ++t)
            {
                siteCoor[env().getNd()-1] = t;
                peekSite(site, prop, siteCoor);
                result.projection[i_p].corr[t]=TensorRemove(site);
                peekSite(site, buf, siteCoor);
                result.projection[i_p].corr_0[t]=TensorRemove(site);
                peekSite(site, propQ, siteCoor);
                result.projection[i_p].corr_Q[t]=TensorRemove(site);
                peekSite(site, propSun, siteCoor);
                result.projection[i_p].corr_Sun[t]=TensorRemove(site);
                peekSite(site, propTad, siteCoor);
                result.projection[i_p].corr_Tad[t]=TensorRemove(site);
            }
        }
        saveResult(par().output, "prop", result);
    }

    std::vector<int> mask(env().getNd(),1);
    mask[env().getNd()-1] = 0;
    fft.FFT_dim_mask(prop, prop, mask, FFT::backward);
    fft.FFT_dim_mask(propQ, propQ, mask, FFT::backward);
    fft.FFT_dim_mask(propSun, propSun, mask, FFT::backward);
    fft.FFT_dim_mask(propTad, propTad, mask, FFT::backward);
}

void TChargedProp::makeCaches(void)
{
    auto &freeMomProp = envGet(ScalarField, freeMomPropName_);
    auto &GFSrc       = envGet(ScalarField, GFSrcName_);
	auto &prop0		  = envGet(ScalarField, prop0Name_);
    auto &fft         = envGet(FFT, fftName_);

    if (!freeMomPropDone_)
    {
        LOG(Message) << "Caching momentum-space free scalar propagator"
                     << " (mass= " << par().mass << ")..." << std::endl;
        SIMPL::MomentumSpacePropagator(freeMomProp, par().mass);
    }
    if (!GFSrcDone_)
    {   
        auto &source = envGet(ScalarField, par().source);
        
        LOG(Message) << "Caching G*F*src..." << std::endl;
        fft.FFT_all_dim(GFSrc, source, FFT::forward);
        GFSrc = freeMomProp*GFSrc;
    }
	if (!prop0Done_)
	{
		LOG(Message) << "Caching position-space free scalar propagator..."
                     << std::endl;
		fft.FFT_all_dim(prop0, GFSrc, FFT::backward);
	}
    if (!phasesDone_)
    {
        std::vector<int> &l = env().getGrid()->_fdimensions;
        Complex          ci(0.0,1.0);
        
        LOG(Message) << "Caching shift phases..." << std::endl;
        for (unsigned int mu = 0; mu < env().getNd(); ++mu)
        {
            Real twoPiL = M_PI*2./l[mu];
            auto &phmu  = envGet(ScalarField, phaseName_[mu]);
            
            LatticeCoordinate(phmu, mu);
            phmu = exp(ci*twoPiL*phmu);
            phase_.push_back(&phmu);
        }
    }
    else
    {
        phase_.clear();
        for (unsigned int mu = 0; mu < env().getNd(); ++mu)
        {
            phase_.push_back(env().getObject<ScalarField>(phaseName_[mu]));
        }
    }
}

void TChargedProp::momD1(ScalarField &s, FFT &fft)
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

void TChargedProp::momD2(ScalarField &s, FFT &fft)
{
    auto &A = envGet(EmField, par().emField);

    envGetTmp(ScalarField, buf);
    envGetTmp(ScalarField, result);
    envGetTmp(ScalarField, Amu);

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
