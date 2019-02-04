/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MScalar/VPCounterTerms.cc

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
#include <Hadrons/Modules/MScalar/VPCounterTerms.hpp>
#include <Hadrons/Modules/MScalar/Scalar.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MScalar;

/******************************************************************************
*                  TVPCounterTerms implementation                             *
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
TVPCounterTerms::TVPCounterTerms(const std::string name)
: Module<VPCounterTermsPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
std::vector<std::string> TVPCounterTerms::getInput(void)
{
    std::vector<std::string> in = {par().source};
    
    return in;
}

std::vector<std::string> TVPCounterTerms::getOutput(void)
{
    std::vector<std::string> out;
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
void TVPCounterTerms::setup(void)
{
	freeMomPropName_ = FREEMOMPROP(par().mass);
    phaseName_.clear();
    for (unsigned int mu = 0; mu < env().getNd(); ++mu)
    {
        phaseName_.push_back("_shiftphase_" + std::to_string(mu));
    }
    GFSrcName_ = getName() + "_DinvSrc";
    phatsqName_ = getName() + "_pHatSquared";
    prop0Name_ = getName() + "_freeProp";
    twoscalarName_ = getName() + "_2scalarProp";
    psquaredName_ = getName() + "_psquaredProp";
    if (!par().output.empty())
    {
        for (unsigned int i_p = 0; i_p < par().outputMom.size(); ++i_p)
        {
            momPhaseName_.push_back("_momentumphase_" + std::to_string(i_p));
        }
    }

    envCreateLat(ScalarField, freeMomPropName_);
    for (unsigned int mu = 0; mu < env().getNd(); ++mu)
    {
        envCreateLat(ScalarField, phaseName_[mu]);
    }
    envCreateLat(ScalarField, phatsqName_);
    envCreateLat(ScalarField, GFSrcName_);
    envCreateLat(ScalarField, prop0Name_);
    envCreateLat(ScalarField, twoscalarName_);
    envCreateLat(ScalarField, psquaredName_);
    if (!par().output.empty())
    {
        for (unsigned int i_p = 0; i_p < par().outputMom.size(); ++i_p)
        {
            envCacheLat(ScalarField, momPhaseName_[i_p]);
        }
    }
    envTmpLat(ScalarField, "buf");
    envTmpLat(ScalarField, "tmp_vp");
    envTmpLat(ScalarField, "vpPhase");
}

// execution ///////////////////////////////////////////////////////////////////
void TVPCounterTerms::execute(void)
{
	auto &source = envGet(ScalarField, par().source);
    Complex     ci(0.0,1.0);
    FFT         fft(env().getGrid());
    envGetTmp(ScalarField, buf);
    envGetTmp(ScalarField, tmp_vp);
    
    // Momentum-space free scalar propagator
    auto &G = envGet(ScalarField, freeMomPropName_);
    SIMPL::MomentumSpacePropagator(G, par().mass);

    // Phases and hat{p}^2
    auto &phatsq = envGet(ScalarField, phatsqName_);
    std::vector<int> &l = env().getGrid()->_fdimensions;
    
    LOG(Message) << "Calculating shift phases..." << std::endl;
    phatsq = zero;
    for (unsigned int mu = 0; mu < env().getNd(); ++mu)
    {
        Real    twoPiL = M_PI*2./l[mu];
        auto &phmu  = envGet(ScalarField, phaseName_[mu]);

        LatticeCoordinate(buf, mu);
        phmu = exp(ci*twoPiL*buf);
        phase_.push_back(&phmu);
        buf = 2.*sin(.5*twoPiL*buf);
		phatsq = phatsq + buf*buf;
    }

    // G*F*src
    auto &GFSrc       = envGet(ScalarField, GFSrcName_);
    fft.FFT_all_dim(GFSrc, source, FFT::forward);
    GFSrc = G*GFSrc;

    // Position-space free scalar propagator
    auto &prop0       = envGet(ScalarField, prop0Name_);
    prop0 = GFSrc;
    fft.FFT_all_dim(prop0, prop0, FFT::backward);

    // Propagators for counter-terms
    auto &twoscalarProp        = envGet(ScalarField, twoscalarName_);
    auto &psquaredProp         = envGet(ScalarField, psquaredName_);

    twoscalarProp = G*GFSrc;
    fft.FFT_all_dim(twoscalarProp, twoscalarProp, FFT::backward);

    psquaredProp = G*phatsq*GFSrc;
    fft.FFT_all_dim(psquaredProp, psquaredProp, FFT::backward);

    // Prepare output data structure if necessary
    Result outputData;
    if (!par().output.empty())
    {
        outputData.projection.resize(par().outputMom.size());
        outputData.lattice_size = env().getGrid()->_fdimensions;
        outputData.mass = par().mass;
        for (unsigned int i_p = 0; i_p < par().outputMom.size(); ++i_p)
        {
            outputData.projection[i_p].momentum = strToVec<int>(par().outputMom[i_p]);
            outputData.projection[i_p].twoScalar.resize(env().getNd());
            outputData.projection[i_p].threeScalar.resize(env().getNd());
            outputData.projection[i_p].pSquaredInsertion.resize(env().getNd());
            for (unsigned int nu = 0; nu < env().getNd(); ++nu)
            {
                outputData.projection[i_p].twoScalar[nu].resize(env().getNd());
                outputData.projection[i_p].threeScalar[nu].resize(env().getNd());
                outputData.projection[i_p].pSquaredInsertion[nu].resize(env().getNd());
            }
            // Calculate phase factors
            auto &momph_ip = envGet(ScalarField, momPhaseName_[i_p]);
            momph_ip = zero;
            for (unsigned int j = 0; j < env().getNd()-1; ++j)
            {
                Real twoPiL = M_PI*2./l[j];
                LatticeCoordinate(buf, j);
                buf = outputData.projection[i_p].momentum[j]*twoPiL*buf;
                momph_ip = momph_ip + buf;
            }
            momph_ip = exp(-ci*momph_ip);
            momPhase_.push_back(&momph_ip);
        }
    }

    // Contractions
    for (unsigned int nu = 0; nu < env().getNd(); ++nu)
    {
    	buf = adj(Cshift(prop0, nu, -1));
        for (unsigned int mu = 0; mu < env().getNd(); ++mu)
        {
            // Two-scalar loop
            tmp_vp = buf * Cshift(prop0, mu, 1);
            tmp_vp -= Cshift(buf, mu, 1) * prop0;
            tmp_vp = 2.0*real(tmp_vp);
            // Output if necessary
            if (!par().output.empty())
            {
                for (unsigned int i_p = 0; i_p < par().outputMom.size(); ++i_p)
                {
                    project(outputData.projection[i_p].twoScalar[mu][nu],
                            tmp_vp, i_p);
                }
            }

        	// Three-scalar loop (no vertex)
    		tmp_vp = buf * Cshift(twoscalarProp, mu, 1);
            tmp_vp -= Cshift(buf, mu, 1) * twoscalarProp;
            tmp_vp = 2.0*real(tmp_vp);
            // Output if necessary
            if (!par().output.empty())
            {
                for (unsigned int i_p = 0; i_p < par().outputMom.size(); ++i_p)
                {
                    project(outputData.projection[i_p].threeScalar[mu][nu],
                            tmp_vp, i_p);
                }
            }

            // Three-scalar loop (hat{p}^2 insertion)
    		tmp_vp = buf * Cshift(psquaredProp, mu, 1);
            tmp_vp -= Cshift(buf, mu, 1) * psquaredProp;
            tmp_vp = 2.0*real(tmp_vp);
            // Output if necessary
            if (!par().output.empty())
            {
                for (unsigned int i_p = 0; i_p < par().outputMom.size(); ++i_p)
                {
                    project(outputData.projection[i_p].pSquaredInsertion[mu][nu],
                            tmp_vp, i_p);
                }
            }
        }
    }

    // OUTPUT IF NECESSARY
    if (!par().output.empty())
    {
        LOG(Message) << "Saving momentum-projected correlators to '"
                     << RESULT_FILE_NAME(par().output, vm().getTrajectory()) << "'..."
                     << std::endl;
        saveResult(par().output, "scalar_loops", outputData);
    }
}

void TVPCounterTerms::project(std::vector<Complex> &projection, const ScalarField &vp, int i_p)
{
    std::vector<TComplex>   vecBuf;
    envGetTmp(ScalarField, vpPhase);

    vpPhase = vp*(*momPhase_[i_p]);
    sliceSum(vpPhase, vecBuf, Tp);
    projection.resize(vecBuf.size());
    for (unsigned int t = 0; t < vecBuf.size(); ++t)
    {
        projection[t] = TensorRemove(vecBuf[t]);
    }
}
