/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MContraction/Gamma3pt.hpp

Copyright (C) 2015-2019

Author: Antonin Portelli <antonin.portelli@me.com>
Author: Lanny91 <andrew.lawson@gmail.com>

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

#ifndef Hadrons_MContraction_Gamma3pt_hpp_
#define Hadrons_MContraction_Gamma3pt_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/*
 * 3pt contraction with gamma matrix insertion.
 *
 * Schematic:
 *
 *             q2           q3
 *        /----<------*------<----¬
 *       /          gamma          \
 *      /                           \
 *   i *                            * f
 *      \                          /
 *       \                        /
 *        \----------->----------/
 *                   q1
 *
 *      trace(g5*q1*adj(q2)*g5*gamma*q3)
 * 
 *  options:
 *   - q1: sink smeared propagator, source at i
 *   - q2: propagator, source at i
 *   - q3: propagator, source at f
 *   - gammas: gamma matrices to insert
 *             (space-separated strings e.g. "GammaT GammaX GammaY") 
 *   - tSnk: sink position for propagator q1.
 *
 */

/******************************************************************************
 *                               Gamma3pt                                     *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class Gamma3ptPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(Gamma3ptPar,
                                    std::string,  q1,
                                    std::string,  q2,
                                    std::string,  q3,
                                    std::string,  gamma,
                                    unsigned int, tSnk,
                                    std::string,  output);
};

template <typename FImpl1, typename FImpl2, typename FImpl3>
class TGamma3pt: public Module<Gamma3ptPar>
{
    FERM_TYPE_ALIASES(FImpl1, 1);
    FERM_TYPE_ALIASES(FImpl2, 2);
    FERM_TYPE_ALIASES(FImpl3, 3);
    class Result: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
                                        Gamma::Algebra, gamma,
                                        std::vector<Complex>, corr);
    };
public:
    // constructor
    TGamma3pt(const std::string name);
    // destructor
    virtual ~TGamma3pt(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual void parseGammaString(std::vector<Gamma::Algebra> &gammaList);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(Gamma3pt, ARG(TGamma3pt<FIMPL, FIMPL, FIMPL>), MContraction);

/******************************************************************************
 *                       TGamma3pt implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2, typename FImpl3>
TGamma3pt<FImpl1, FImpl2, FImpl3>::TGamma3pt(const std::string name)
: Module<Gamma3ptPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2, typename FImpl3>
std::vector<std::string> TGamma3pt<FImpl1, FImpl2, FImpl3>::getInput(void)
{
    std::vector<std::string> in = {par().q1, par().q2, par().q3};
    
    return in;
}

template <typename FImpl1, typename FImpl2, typename FImpl3>
std::vector<std::string> TGamma3pt<FImpl1, FImpl2, FImpl3>::getOutput(void)
{
    std::vector<std::string> out = {};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2, typename FImpl3>
void TGamma3pt<FImpl1, FImpl2, FImpl3>::setup(void)
{
    envTmpLat(LatticeComplex, "c");
}

template <typename FImpl1, typename FImpl2, typename FImpl3>
void TGamma3pt<FImpl1, FImpl2, FImpl3>::parseGammaString(std::vector<Gamma::Algebra> &gammaList)
{
    gammaList.clear();
    // Determine gamma matrices to insert at source/sink.
    if (par().gamma.compare("all") == 0)
    {
        // Do all contractions.
        for (unsigned int i = 1; i < Gamma::nGamma; i += 2)
        {
            gammaList.push_back((Gamma::Algebra)i);
        }
    }
    else
    {
        // Parse individual contractions from input string.
        gammaList = strToVec<Gamma::Algebra>(par().gamma);
    } 
}
// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2, typename FImpl3>
void TGamma3pt<FImpl1, FImpl2, FImpl3>::execute(void)
{
    LOG(Message) << "Computing 3pt contractions '" << getName() << "' using"
                 << " quarks '" << par().q1 << "', '" << par().q2 << "' and '"
                 << par().q3 << "', with " << par().gamma << " insertions." 
                 << std::endl;

    // Initialise variables. q2 and q3 are normal propagators, q1 may be 
    // sink smeared.
    auto                        &q1 = envGet(SlicedPropagator1, par().q1);
    auto                        &q2 = envGet(PropagatorField2, par().q2);
    auto                        &q3 = envGet(PropagatorField2, par().q3);
    Gamma                       g5(Gamma::Algebra::Gamma5);
    std::vector<Gamma::Algebra> gammaList;
    std::vector<TComplex>       buf;
    std::vector<Result>         result;
    int                         nt = env().getDim(Tp);


    parseGammaString(gammaList);
    result.resize(gammaList.size());
    for (unsigned int i = 0; i < result.size(); ++i)
    {
        result[i].gamma = gammaList[i];
        result[i].corr.resize(nt);
    }
    
    // Extract relevant timeslice of sinked propagator q1, then contract &
    // sum over all spacial positions of gamma insertion.
    SitePropagator1 q1Snk = q1[par().tSnk];
    envGetTmp(LatticeComplex, c);
    for (unsigned int i = 0; i < result.size(); ++i)
    {
        Gamma gamma(gammaList[i]);
        c = trace(g5*q1Snk*adj(q2)*(g5*gamma)*q3);
        sliceSum(c, buf, Tp);
        for (unsigned int t = 0; t < buf.size(); ++t)
        {
            result[i].corr[t] = TensorRemove(buf[t]);
        }
    }
    saveResult(par().output, "gamma3pt", result);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_Gamma3pt_hpp_
