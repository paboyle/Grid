/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: extras/Hadrons/Modules/MScalar/ChargedProp.cc

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
    std::vector<std::string> out = {getName()};
    
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
    fftName_   = getName() + "_fft";

    freeMomPropDone_ = env().hasCreatedObject(freeMomPropName_);
    GFSrcDone_       = env().hasCreatedObject(GFSrcName_);
    phasesDone_      = env().hasCreatedObject(phaseName_[0]);
    envCacheLat(ScalarField, freeMomPropName_);
    for (unsigned int mu = 0; mu < env().getNd(); ++mu)
    {
        envCacheLat(ScalarField, phaseName_[mu]);
    }
    envCacheLat(ScalarField, GFSrcName_);
    envCreateLat(ScalarField, getName());
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
    
    auto   &prop  = envGet(ScalarField, getName());
    auto   &GFSrc = envGet(ScalarField, GFSrcName_);
    auto   &G     = envGet(ScalarField, freeMomPropName_);
    auto   &fft   = envGet(FFT, fftName_);
    double q      = par().charge;
    envGetTmp(ScalarField, result); 
    envGetTmp(ScalarField, buf); 

    // G*F*Src
    prop = GFSrc;

    // - q*G*momD1*G*F*Src (momD1 = F*D1*Finv)
    buf = GFSrc;
    momD1(buf, fft);
    buf = G*buf;
    prop = prop - q*buf;

    // + q^2*G*momD1*G*momD1*G*F*Src (here buf = G*momD1*G*F*Src)
    momD1(buf, fft);
    prop = prop + q*q*G*buf;

    // - q^2*G*momD2*G*F*Src (momD2 = F*D2*Finv)
    buf = GFSrc;
    momD2(buf, fft);
    prop = prop - q*q*G*buf;

    // final FT
    fft.FFT_all_dim(prop, prop, FFT::backward);
    
    // OUTPUT IF NECESSARY
    if (!par().output.empty())
    {
        std::string           filename = par().output + "." +
                                         std::to_string(vm().getTrajectory());
        
        LOG(Message) << "Saving zero-momentum projection to '"
                     << filename << "'..." << std::endl;
        
        CorrWriter            writer(filename);
        std::vector<TComplex> vecBuf;
        std::vector<Complex>  result;
        
        sliceSum(prop, vecBuf, Tp);
        result.resize(vecBuf.size());
        for (unsigned int t = 0; t < vecBuf.size(); ++t)
        {
            result[t] = TensorRemove(vecBuf[t]);
        }
        write(writer, "charge", q);
        write(writer, "prop", result);
    }
}

void TChargedProp::makeCaches(void)
{
    auto &freeMomProp = envGet(ScalarField, freeMomPropName_);
    auto &GFSrc       = envGet(ScalarField, GFSrcName_);
    auto &fft         = envGet(FFT, fftName_);

    if (!freeMomPropDone_)
    {
        LOG(Message) << "Caching momentum space free scalar propagator"
                     << " (mass= " << par().mass << ")..." << std::endl;
        SIMPL::MomentumSpacePropagator(freeMomProp, par().mass);
    }
    if (!GFSrcDone_)
    {   
        FFT  fft(env().getGrid());
        auto &source = envGet(ScalarField, par().source);

        LOG(Message) << "Caching G*F*src..." << std::endl;
        fft.FFT_all_dim(GFSrc, source, FFT::forward);
        GFSrc = freeMomProp*GFSrc;
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
