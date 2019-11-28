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
 *           qs_xi         qd_xf           |            qs_xi          qd_xf
 *          /--<--¬       /--<--¬          |          /--<--¬         /--<--¬       
 *         /       \     /       \         |         /       \       /       \      
 *        /         \   /         \        |        /         \     /         \     
 *       /           \ /           \       |       /           \   /           \    
 *      *             * G           *      |       *           G * * G          * 
 *     |\             * G           |      |      |\           /   \           /|
 *     | \           / \           /|      |      | \         /     \         / |   
 *     |  \         /   \         / |      |      |  \       /       \       /  |
 *     |   \       /     \       /  |      |      |   \-->--/         \-->--/   |   
 *      \   \-->--/       \-->--/  /       |       \   qu_xi           qu_xf   /
 *       \    qu_xf         qu_xf /        |        \                         /
 *        \                      /         |         \                       /
 *         \--------->----------/          |          \--------->-----------/
 *                 qu_spec                 |                  qu_spec
 *
 *
 * analogously to the rare-kaon naming, the left diagram is named 'one-trace' and
 * the diagram on the right 'two-trace'
 * 
 * Propagators:
 *  * qu_xi,   source at xi 
 *  * qu_xf,   source at xf
 *  * qu_spec, source at xi
 *  * qd_xf,   source at xf 
 *  * qs_xi,   source at xi
 */
BEGIN_MODULE_NAMESPACE(MContraction)

class SigmaToNucleonNonEyePar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(SigmaToNucleonNonEyePar,
                                    std::string, qu_xi,
                                    std::string, qu_xf,
                                    std::string, qu_spec,
                                    std::string, qd_xf,
                                    std::string, qs_xi,
                                    unsigned int,   xf,
                                    std::string, parity,
                                    std::string, sink,
                                    std::string, output);
};

template <typename FImpl1, typename FImpl2, typename FImpl3, typename FImpl4>
class TSigmaToNucleonNonEye: public Module<SigmaToNucleonNonEyePar>
{
public:
    FERM_TYPE_ALIASES(FImpl1, 1);
    FERM_TYPE_ALIASES(FImpl2, 2);
    FERM_TYPE_ALIASES(FImpl3, 3);
    FERM_TYPE_ALIASES(FImpl3, 4);
    BASIC_TYPE_ALIASES(ScalarImplCR, Scalar);
    SINK_TYPE_ALIASES(Scalar);
    typedef typename SpinMatrixField1::vector_object::scalar_object SpinMatrix;
    class Metadata: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Metadata,
                                        Gamma::Algebra, gamma_H,
                                        Gamma::Algebra, gammaA_sigma,
                                        Gamma::Algebra, gammaB_sigma,
                                        Gamma::Algebra, gammaA_nucl,
                                        Gamma::Algebra, gammaB_nucl,
                                        int, trace,
                                        int, parity);
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

MODULE_REGISTER_TMP(SigmaToNucleonNonEye, ARG(TSigmaToNucleonNonEye<FIMPL, FIMPL, FIMPL, FIMPL>), MContraction);

/******************************************************************************
 *                         TSigmaToNucleonNonEye implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2, typename FImpl3, typename FImpl4>
TSigmaToNucleonNonEye<FImpl1, FImpl2, FImpl3, FImpl4>::TSigmaToNucleonNonEye(const std::string name)
: Module<SigmaToNucleonNonEyePar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2, typename FImpl3, typename FImpl4>
std::vector<std::string> TSigmaToNucleonNonEye<FImpl1, FImpl2, FImpl3, FImpl4>::getInput(void)
{
    std::vector<std::string> input = {par().qu_xi, par().qu_xf, par().qu_spec, par().qd_xf, par().qs_xi, par().sink};
    
    return input;
}

template <typename FImpl1, typename FImpl2, typename FImpl3, typename FImpl4>
std::vector<std::string> TSigmaToNucleonNonEye<FImpl1, FImpl2, FImpl3, FImpl4>::getOutput(void)
{
    std::vector<std::string> out = {};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2, typename FImpl3, typename FImpl4>
void TSigmaToNucleonNonEye<FImpl1, FImpl2, FImpl3, FImpl4>::setup(void)
{
    envTmpLat(SpinMatrixField1, "c");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2, typename FImpl3, typename FImpl4>
void TSigmaToNucleonNonEye<FImpl1, FImpl2, FImpl3, FImpl4>::execute(void)
{
    const int  parity {par().parity.size()>0 ? std::stoi(par().parity) : 1};
    const Gamma GammaB(Gamma::Algebra::SigmaXZ); // C*gamma_5
    const Gamma Id(Gamma::Algebra::Identity); // C*gamma_5

    LOG(Message) << "Computing sigma-to-nucleon contractions '" << getName() << "'" << std::endl;
    LOG(Message) << "' with (Gamma^A,Gamma^B)_sigma = ( Identity, C*gamma_5 ) and (Gamma^A,Gamma^B)_nucl = ( Identity, C*gamma_5 )" << std::endl; 
    LOG(Message) << "and parity " << parity << " using sink " << par().sink << "." << std::endl;
        
    envGetTmp(SpinMatrixField1, c);
    std::vector<SpinMatrix> buf;

    std::vector<Result> result;
    Result              r;
    r.info.parity       = parity;
    r.info.gammaA_sigma = Id.g;
    r.info.gammaB_sigma = GammaB.g;
    r.info.gammaA_nucl  = Id.g;
    r.info.gammaB_nucl  = GammaB.g;

    auto &qu_xi      = envGet(PropagatorField1, par().qu_xi);
    auto &qu_xf      = envGet(PropagatorField1, par().qu_xf);
    auto &qu_spec    = envGet(SlicedPropagator2, par().qu_spec);
    auto &qd_xf      = envGet(PropagatorField3, par().qd_xf);
    auto &qs_xi      = envGet(PropagatorField4, par().qs_xi);
    auto qut         = qu_spec[par().xf];
    for (auto &G: Gamma::gall)
    {
      r.info.gamma_H = G.g;
      //Operator Q1, equivalent to the two-trace case in the rare-kaons module
      c=Zero();
      BaryonUtils<FIMPL>::Sigma_to_Nucleon_NonEye(qu_xi,qu_xf,qut,qd_xf,qs_xi,G,GammaB,GammaB,parity,"Q1",c);
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
      BaryonUtils<FIMPL>::Sigma_to_Nucleon_NonEye(qu_xi,qu_xf,qut,qd_xf,qs_xi,G,GammaB,GammaB,parity,"Q2",c);
      sliceSum(c,buf,Tp);
      r.corr.clear();
      for (unsigned int t = 0; t < buf.size(); ++t)
      {
          r.corr.push_back(buf[t]);
      }
      r.info.trace = 1;
      result.push_back(r);
    }

    saveResult(par().output, "StN_NonEye", result);

}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_SigmaToNucleonNonEye_hpp_
