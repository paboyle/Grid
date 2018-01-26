#include <Grid/Hadrons/Modules/MScalar/VPCounterTerms.hpp>
#include <Grid/Hadrons/Modules/MScalar/Scalar.hpp>

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
    twoscalarVertexName_ = getName() + "_2scalarProp_withvertex";
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
    envCreateLat(ScalarField, twoscalarVertexName_);
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
    auto &twoscalarVertexProp  = envGet(ScalarField, twoscalarVertexName_);
    auto &psquaredProp         = envGet(ScalarField, psquaredName_);

    twoscalarProp = G*GFSrc;
    fft.FFT_all_dim(twoscalarProp, twoscalarProp, FFT::backward);

    twoscalarVertexProp = zero;
    for (unsigned int mu = 0; mu < env().getNd(); ++mu)
    {
        buf = GFSrc;
        twoscalarVertexProp = twoscalarVertexProp + .5*((*phase_[mu]) + adj(*phase_[mu]))*buf;
    }
    twoscalarVertexProp = G*twoscalarVertexProp;
    fft.FFT_all_dim(twoscalarVertexProp, twoscalarVertexProp, FFT::backward);

    psquaredProp = G*phatsq*GFSrc;
    fft.FFT_all_dim(psquaredProp, psquaredProp, FFT::backward);

    // Prepare output files if necessary
    if (!par().output.empty())
    {
        LOG(Message) << "Preparing output files..." << std::endl;
        for (unsigned int i_p = 0; i_p < par().outputMom.size(); ++i_p)
        {
            std::vector<int> mom = strToVec<int>(par().outputMom[i_p]);

            // Open output files
            std::string           filename = par().output + "_" + std::to_string(mom[0])
                                                          + std::to_string(mom[1])
                                                          + std::to_string(mom[2]);
            saveResult(filename, "mass", par().mass);

            // Calculate phase factors
            auto &momph_ip = envGet(ScalarField, momPhaseName_[i_p]);
            momph_ip = Complex(1.0,0.0);
            for (unsigned int j = 0; j < env().getNd()-1; ++j)
            {
                for (unsigned int momcount = 0; momcount < mom[j]; ++momcount)
                {
                    momph_ip = momph_ip*(*phase_[j]);
                }
            }
            momph_ip = adj(momph_ip);
            momPhase_.push_back(&momph_ip);
        }
    }

    // Contractions
    for (unsigned int nu = 0; nu < env().getNd(); ++nu)
    {
    	buf = adj(Cshift(prop0, nu, -1));
        for (unsigned int mu = 0; mu < env().getNd(); ++mu)
        {
        	// Three-scalar loop (no vertex)
    		tmp_vp = buf * Cshift(twoscalarProp, mu, 1);
            tmp_vp -= Cshift(buf, mu, 1) * twoscalarProp;
            tmp_vp = 2.0*real(tmp_vp);

            // Output if necessary
            if (!par().output.empty())
            {
                writeVP(tmp_vp, "NoVertex_"+std::to_string(mu)+"_"+std::to_string(nu));
            }

            // Three-scalar loop (tadpole vertex)
    		tmp_vp = buf * Cshift(twoscalarVertexProp, mu, 1);
            tmp_vp -= Cshift(buf, mu, 1) * twoscalarVertexProp;
            tmp_vp = 2.0*real(tmp_vp);

            // Output if necessary
            if (!par().output.empty())
            {
                writeVP(tmp_vp, "TadVertex_"+std::to_string(mu)+"_"+std::to_string(nu));
            }

            // Three-scalar loop (hat{p}^2 insertion)
    		tmp_vp = buf * Cshift(psquaredProp, mu, 1);
            tmp_vp -= Cshift(buf, mu, 1) * psquaredProp;
            tmp_vp = 2.0*real(tmp_vp);

            // Output if necessary
            if (!par().output.empty())
            {
                writeVP(tmp_vp, "pSquaredInsertion_"+std::to_string(mu)+"_"+std::to_string(nu));
            }
        }
    }
}

void TVPCounterTerms::writeVP(const ScalarField &vp, std::string dsetName)
{
    std::vector<TComplex>   vecBuf;
    std::vector<Complex>    result;
    envGetTmp(ScalarField, vpPhase);

    for (unsigned int i_p = 0; i_p < par().outputMom.size(); ++i_p)
    {
        std::vector<int> mom = strToVec<int>(par().outputMom[i_p]);
        std::string filename = par().output + "_"
                               + std::to_string(mom[0])
                               + std::to_string(mom[1])
                               + std::to_string(mom[2]);
        vpPhase = vp*(*momPhase_[i_p]);
        sliceSum(vpPhase, vecBuf, Tp);
        result.resize(vecBuf.size());
        for (unsigned int t = 0; t < vecBuf.size(); ++t)
        {
            result[t] = TensorRemove(vecBuf[t]);
        }
        saveResult(filename, dsetName, result);
    }
}
