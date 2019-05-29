/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MContraction/WeakMesonDecayKl2.hpp

Copyright (C) 2015-2018

Author: Antonin Portelli <antonin.portelli@me.com>
Author: Vera Guelpers <Vera.Guelpers@ed.ac.uk>


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

#ifndef Hadrons_MContraction_WeakMesonDecayKl2_hpp_
#define Hadrons_MContraction_WeakMesonDecayKl2_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/*
* Kl2 contraction
* -----------------------------
*
* contraction for Kl2 decay, including the lepton
*
* 	trace(q1*adj(q2)*g5*gL[mu]) * (gL[mu] * lepton)_{a,b}
*
* with open spinor indices (a,b) for the lepton part
*
*             q1                  lepton
*        /------------\       /------------
*       /              \     /
*      /                \H_W/
* g_5 *                  * * 
*      \                /
*       \              / 
*        \____________/
*             q2
*
* * options:
* - q1: input propagator 1 (string)
* - q2: input propagator 2 (string)
* - lepton: input lepton (string)
*/

/******************************************************************************
 *                                TWeakMesonDecayKl2                             *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class WeakMesonDecayKl2Par: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(WeakMesonDecayKl2Par,
                                    std::string, q1,
                                    std::string, q2,
                                    std::string, lepton,
				    std::string, output);
};

template <typename FImpl>
class TWeakMesonDecayKl2: public Module<WeakMesonDecayKl2Par>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    class Metadata: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Metadata,
                                        int, spinidx1,
                                        int, spinidx2);
    };
    typedef Correlator<Metadata> Result;
public:
    // constructor
    TWeakMesonDecayKl2(const std::string name);
    // destructor
    virtual ~TWeakMesonDecayKl2(void) {};
    // dependencies/products
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // execution
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(WeakMesonDecayKl2, TWeakMesonDecayKl2<FIMPL>, MContraction);

/******************************************************************************
 *                           TWeakMesonDecayKl2 implementation                   *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TWeakMesonDecayKl2<FImpl>::TWeakMesonDecayKl2(const std::string name)
: Module<WeakMesonDecayKl2Par>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TWeakMesonDecayKl2<FImpl>::getInput(void)
{
    std::vector<std::string> input = {par().q1, par().q2, par().lepton};
    
    return input;
}

template <typename FImpl>
std::vector<std::string> TWeakMesonDecayKl2<FImpl>::getOutput(void)
{
    std::vector<std::string> output = {};
    
    return output;
}

// setup ////////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TWeakMesonDecayKl2<FImpl>::setup(void)
{
    envTmpLat(LatticeComplex, "c");
    envTmpLat(PropagatorField, "prop_buf");
    envCreateLat(PropagatorField, getName());
    envTmpLat(LatticeComplex, "buf");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TWeakMesonDecayKl2<FImpl>::execute(void)
{
    LOG(Message) << "Computing QED Kl2 contractions '" << getName() << "' using"
                 << " quarks '" << par().q1 << "' and '" << par().q2 << "' and"
		 << "lepton '"  << par().lepton << "'" << std::endl;


    auto &res = envGet(PropagatorField, getName()); res = zero;
    Gamma                  g5(Gamma::Algebra::Gamma5);
    int                    nt = env().getDim(Tp);

    auto &q1 = envGet(PropagatorField, par().q1);
    auto &q2 = envGet(PropagatorField, par().q2);
    auto &lepton = envGet(PropagatorField, par().lepton);
    envGetTmp(LatticeComplex, buf);
    std::vector<TComplex>  res_summed;
    envGetTmp(LatticeComplex, c);
    envGetTmp(PropagatorField, prop_buf);  

    std::vector<Result>    result;
    Result r;

    for (unsigned int mu = 0; mu < 4; ++mu)
    {
	c = zero;
	//hadronic part: trace(q1*adj(q2)*g5*gL[mu]) 
        c   = trace(q1*adj(q2)*g5*GammaL(Gamma::gmu[mu]));
    	prop_buf = 1.;
	//multiply lepton part
	res += c * prop_buf * GammaL(Gamma::gmu[mu]) * lepton;
    }

    //loop over spinor index of lepton part
    unsigned int i = 0;
    for (unsigned int s1 = 0; s1 < Ns ; ++s1)
    for (unsigned int s2 = 0; s2 < Ns ; ++s2)
    {
	buf = peekColour(peekSpin(res,s1,s2),0,0);

	sliceSum(buf, res_summed, Tp);

	r.corr.clear();
	for (unsigned int t = 0; t < nt; ++t)
	{
              r.corr.push_back(TensorRemove(res_summed[t]));
	}

	r.info.spinidx1 = s1;
	r.info.spinidx2 = s2;
	result.push_back(r);

	i+=1;
    }

    saveResult(par().output, "weakdecay", result);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_WeakMesonDecayKl2_hpp_
