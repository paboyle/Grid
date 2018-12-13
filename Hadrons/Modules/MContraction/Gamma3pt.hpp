/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MContraction/Gamma3pt.hpp

Copyright (C) 2015-2018

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
 *   - gamma: gamma matrix to insert
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
                                    std::string,    q1,
                                    std::string,    q2,
                                    std::string,    q3,
                                    Gamma::Algebra, gamma,
                                    unsigned int,   tSnk,
                                    std::string,    output);
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

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2, typename FImpl3>
void TGamma3pt<FImpl1, FImpl2, FImpl3>::execute(void)
{
    LOG(Message) << "Computing 3pt contractions '" << getName() << "' using"
                 << " quarks '" << par().q1 << "', '" << par().q2 << "' and '"
                 << par().q3 << "', with " << par().gamma << " insertion." 
                 << std::endl;

    // Initialise variables. q2 and q3 are normal propagators, q1 may be 
    // sink smeared.
    auto                  &q1 = envGet(SlicedPropagator1, par().q1);
    auto                  &q2 = envGet(PropagatorField2, par().q2);
    auto                  &q3 = envGet(PropagatorField2, par().q3);
    Gamma                 g5(Gamma::Algebra::Gamma5);
    Gamma                 gamma(par().gamma);
    std::vector<TComplex> buf;
    Result                result;
    
    // Extract relevant timeslice of sinked propagator q1, then contract &
    // sum over all spacial positions of gamma insertion.
    SitePropagator1 q1Snk = q1[par().tSnk];
    envGetTmp(LatticeComplex, c);
    c = trace(g5*q1Snk*adj(q2)*(g5*gamma)*q3);
    sliceSum(c, buf, Tp);

    result.gamma = par().gamma;
    result.corr.resize(buf.size());
    for (unsigned int t = 0; t < buf.size(); ++t)
    {
        result.corr[t] = TensorRemove(buf[t]);
    }
    saveResult(par().output, "gamma3pt", result);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_Gamma3pt_hpp_
