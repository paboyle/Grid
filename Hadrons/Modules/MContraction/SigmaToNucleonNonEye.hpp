/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MContraction/SigmaToNucleonNonEye.hpp

Copyright (C) 2015-2019

Author: Antonin Portelli <antonin.portelli@me.com>
Author: Felix Erben <felix.erben@ed.ac.uk>

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

#ifndef Hadrons_MContraction_SigmaToNucleonNonEye_hpp_
#define Hadrons_MContraction_SigmaToNucleonNonEye_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Grid/qcd/utils/BaryonUtils.h>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                               SigmaToNucleonNonEye                                       *
 ******************************************************************************/
/*
 * Sigma-to-Nucleon 3-pt diagrams, non-eye topologies.
 * 
 * Schematic:     
 *            qsTi          quTf           |            qsTi            qdTf
 *          /-->--¬       /-->--¬          |          /-->--¬         /-->--¬       
 *         /       \     /       \         |         /       \       /       \      
 *        /         \   /         \        |        /         \     /         \     
 *       /           \ /           \       |       /           \   /           \    
 *      *             * G           *      |       *           G * * G          * 
 *     |\             * G           |      |      |\           /   \           /|
 *     | \           / \           /|      |      | \         /     \         / |   
 *     |  \         /   \         / |      |      |  \       /       \       /  |
 *     |   \       /     \       /  |      |      |   \-->--/         \-->--/   |   
 *      \   \-->--/       \-->--/  /       |       \   quTi            quTf    /
 *       \    quTi          qdTf  /        |        \                         /
 *        \                      /         |         \                       /
 *         \--------->----------/          |          \--------->-----------/
 *                 quSpec                  |                  quSpec
 *
 *
 * analogously to the rare-kaon naming, the left diagram is named 'one-trace' and
 * the diagram on the right 'two-trace'
 * 
 * Propagators:
 *  * quTi,   source at ti 
 *  * quTf,   source at tf
 *  * quSpec, source at ti
 *  * qdTf,   source at tf 
 *  * qsTi,   source at ti
 */
BEGIN_MODULE_NAMESPACE(MContraction)

class SigmaToNucleonNonEyePar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(SigmaToNucleonNonEyePar,
                                    std::string, quTi,
                                    std::string, quTf,
                                    std::string, quSpec,
                                    std::string, qdTf,
                                    std::string, qsTi,
                                    unsigned int,   tf,
                                    std::string, sink,
                                    std::string, output);
};

template <typename FImpl>
class TSigmaToNucleonNonEye: public Module<SigmaToNucleonNonEyePar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    BASIC_TYPE_ALIASES(ScalarImplCR, Scalar);
    SINK_TYPE_ALIASES(Scalar);
    typedef typename SpinMatrixField::vector_object::scalar_object SpinMatrix;
    class Metadata: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Metadata,
                                        Gamma::Algebra, gammaH,
                                        Gamma::Algebra, gammaASigma,
                                        Gamma::Algebra, gammaBSigma,
                                        Gamma::Algebra, gammaANucl,
                                        Gamma::Algebra, gammaBNucl,
                                        int, trace);
    };
    typedef Correlator<Metadata, SpinMatrix> Result;
public:
    // constructor
    TSigmaToNucleonNonEye(const std::string name);
    // destructor
    virtual ~TSigmaToNucleonNonEye(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
    // Which gamma algebra was specified
    Gamma::Algebra  al;
};

MODULE_REGISTER_TMP(SigmaToNucleonNonEye, ARG(TSigmaToNucleonNonEye<FIMPL>), MContraction);

/******************************************************************************
 *                         TSigmaToNucleonNonEye implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TSigmaToNucleonNonEye<FImpl>::TSigmaToNucleonNonEye(const std::string name)
: Module<SigmaToNucleonNonEyePar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TSigmaToNucleonNonEye<FImpl>::getInput(void)
{
    std::vector<std::string> input = {par().quTi, par().quTf, par().quSpec, par().qdTf, par().qsTi, par().sink};
    
    return input;
}

template <typename FImpl>
std::vector<std::string> TSigmaToNucleonNonEye<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TSigmaToNucleonNonEye<FImpl>::setup(void)
{
    envTmpLat(SpinMatrixField, "c");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TSigmaToNucleonNonEye<FImpl>::execute(void)
{
    const Gamma GammaB(Gamma::Algebra::SigmaXZ); // C*gamma_5
    const Gamma Id(Gamma::Algebra::Identity); // C*gamma_5

    LOG(Message) << "Computing sigma-to-nucleon contractions '" << getName() << "'" << std::endl;
    LOG(Message) << "' with (Gamma^A,Gamma^B)_sigma = ( Identity, C*gamma_5 ) and (Gamma^A,Gamma^B)_nucl = ( Identity, C*gamma_5 )" << std::endl; 
    LOG(Message) << " using sink " << par().sink << "." << std::endl;
        
    envGetTmp(SpinMatrixField, c);
    std::vector<SpinMatrix> buf;

    std::vector<Result> result;
    Result              r;
    r.info.gammaASigma = Id.g;
    r.info.gammaBSigma = GammaB.g;
    r.info.gammaANucl  = Id.g;
    r.info.gammaBNucl  = GammaB.g;

    auto &quTi      = envGet(PropagatorField, par().quTi);
    auto &quTf      = envGet(PropagatorField, par().quTf);
    auto &quSpec    = envGet(SlicedPropagator, par().quSpec);
    auto &qdTf      = envGet(PropagatorField, par().qdTf);
    auto &qsTi      = envGet(PropagatorField, par().qsTi);
    auto qut         = quSpec[par().tf];
    for (auto &G: Gamma::gall)
    {
      r.info.gammaH = G.g;
      //Operator Q1, equivalent to the two-trace case in the rare-kaons module
      c=Zero();
      BaryonUtils<FIMPL>::Sigma_to_Nucleon_NonEye(quTi,quTf,qut,qdTf,qsTi,G,GammaB,GammaB,"Q1",c);
      sliceSum(c,buf,Tp);
      r.corr.clear();
      for (unsigned int t = 0; t < buf.size(); ++t)
      {
          r.corr.push_back(buf[t]);
      }
      r.info.trace = 2;
      result.push_back(r);
      //Operator Q2, equivalent to the one-trace case in the rare-kaons module
      c=Zero();
      BaryonUtils<FIMPL>::Sigma_to_Nucleon_NonEye(quTi,quTf,qut,qdTf,qsTi,G,GammaB,GammaB,"Q2",c);
      sliceSum(c,buf,Tp);
      r.corr.clear();
      for (unsigned int t = 0; t < buf.size(); ++t)
      {
          r.corr.push_back(buf[t]);
      }
      r.info.trace = 1;
      result.push_back(r);
    }

    saveResult(par().output, "stnNonEye", result);

}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_SigmaToNucleonNonEye_hpp_
