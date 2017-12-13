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
    GFSrcName_ = "_" + getName() + "_DinvSrc";
    phatsqName_ = "_" + getName() + "_pHatSquared";
    prop0Name_ = getName() + "_freeProp";
    twoscalarName_ = getName() + "_2scalarProp";
    twoscalarVertexName_ = getName() + "_2scalarProp_withvertex";
    psquaredName_ = getName() + "_psquaredProp";
    env().registerLattice<ScalarField>(freeMomPropName_);
    for (unsigned int mu = 0; mu < env().getNd(); ++mu)
    {
        env().registerLattice<ScalarField>(phaseName_[mu]);
    }
    env().registerLattice<ScalarField>(phatsqName_);
    env().registerLattice<ScalarField>(GFSrcName_);
    env().registerLattice<ScalarField>(prop0Name_);
    env().registerLattice<ScalarField>(twoscalarName_);
    env().registerLattice<ScalarField>(twoscalarVertexName_);
    env().registerLattice<ScalarField>(psquaredName_);
}

// execution ///////////////////////////////////////////////////////////////////
void TVPCounterTerms::execute(void)
{
	ScalarField &source = *env().getObject<ScalarField>(par().source);
    Complex     ci(0.0,1.0);
    FFT         fft(env().getGrid());
    ScalarField     buf(env().getGrid()), tmp_vp(env().getGrid());
    
    // Momentum-space free scalar propagator
    ScalarField &G = *env().createLattice<ScalarField>(freeMomPropName_);
    SIMPL::MomentumSpacePropagator(G, par().mass);

    // Phases and hat{p}^2
    ScalarField &phatsq = *env().createLattice<ScalarField>(phatsqName_);
    std::vector<int> &l = env().getGrid()->_fdimensions;
    
    LOG(Message) << "Calculating shift phases..." << std::endl;
    phatsq = zero;
    for (unsigned int mu = 0; mu < env().getNd(); ++mu)
    {
        Real    twoPiL = M_PI*2./l[mu];
        
        phase_.push_back(env().createLattice<ScalarField>(phaseName_[mu]));
        LatticeCoordinate(buf, mu);
        *(phase_[mu]) = exp(ci*twoPiL*buf);
        buf = 2.*sin(.5*twoPiL*buf);
		phatsq = phatsq + buf*buf;
    }

    // G*F*src
    ScalarField &GFSrc = *env().createLattice<ScalarField>(GFSrcName_);
    fft.FFT_all_dim(GFSrc, source, FFT::forward);
    GFSrc = G*GFSrc;

    // Position-space free scalar propagator
    ScalarField &prop0 = *env().createLattice<ScalarField>(prop0Name_);
    prop0 = GFSrc;
    fft.FFT_all_dim(prop0, prop0, FFT::backward);

    // Propagators for counter-terms
    ScalarField &twoscalarProp   = *env().createLattice<ScalarField>(twoscalarName_);
    ScalarField &twoscalarVertexProp  = *env().createLattice<ScalarField>(twoscalarVertexName_);
    ScalarField &psquaredProp   = *env().createLattice<ScalarField>(psquaredName_);

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

    // Open output files if necessary
    std::vector<TComplex>   vecBuf;
    std::vector<Complex>    result;
    ScalarField vpPhase(env().getGrid());
    std::vector<CorrWriter *> writer;
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

            if (env().getGrid()->IsBoss())
            {
                CorrWriter *writer_i = new CorrWriter(filename);
                writer.push_back(writer_i);

                write(*writer[i_p], "mass", par().mass);
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
                        write(*writer[i_p],
                              "NoVertex_"+std::to_string(mu)+"_"+std::to_string(nu),
                              result);
                    }
                }
            }

            // Three-scalar loop (tadpole vertex)
    		tmp_vp = buf * Cshift(twoscalarVertexProp, mu, 1);
            tmp_vp -= Cshift(buf, mu, 1) * twoscalarVertexProp;
            tmp_vp = 2.0*real(tmp_vp);

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
                        write(*writer[i_p],
                              "TadVertex_"+std::to_string(mu)+"_"+std::to_string(nu),
                              result);
                    }
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
                    vpPhase = tmp_vp*momphases[i_p];
                    sliceSum(vpPhase, vecBuf, Tp);
                    result.resize(vecBuf.size());
                    for (unsigned int t = 0; t < vecBuf.size(); ++t)
                    {
                        result[t] = TensorRemove(vecBuf[t]);
                    }
                    if (env().getGrid()->IsBoss())
                    {
                        write(*writer[i_p],
                              "pSquaredInsertion_"+std::to_string(mu)+"_"+std::to_string(nu),
                              result);
                    }
                }
            }
        }
    }

    // Close output files if necessary
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
