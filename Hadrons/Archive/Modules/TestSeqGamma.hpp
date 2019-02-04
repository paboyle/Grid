/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MUtilities/TestSeqGamma.hpp

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

#ifndef Hadrons_MUtilities_TestSeqGamma_hpp_
#define Hadrons_MUtilities_TestSeqGamma_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                              TestSeqGamma                                  *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)

class TestSeqGammaPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(TestSeqGammaPar,
                                    std::string,    q,
                                    std::string,    qSeq,
                                    std::string,    origin,
                                    Gamma::Algebra, gamma,
                                    unsigned int,   t_g);
};

template <typename FImpl>
class TTestSeqGamma: public Module<TestSeqGammaPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TTestSeqGamma(const std::string name);
    // destructor
    virtual ~TTestSeqGamma(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(TestSeqGamma, TTestSeqGamma<FIMPL>, MUtilities);

/******************************************************************************
 *                      TTestSeqGamma implementation                          *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TTestSeqGamma<FImpl>::TTestSeqGamma(const std::string name)
: Module<TestSeqGammaPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TTestSeqGamma<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().q, par().qSeq};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TTestSeqGamma<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TTestSeqGamma<FImpl>::setup(void)
{
    envTmpLat(LatticeComplex, "c");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TTestSeqGamma<FImpl>::execute(void)
{
    auto                  &q    = envGet(PropagatorField, par().q);
    auto                  &qSeq = envGet(PropagatorField, par().qSeq);
    Gamma                 g5(Gamma::Algebra::Gamma5);
    Gamma                 g(par().gamma);
    SitePropagator        qSite;
    Complex               test, check;
    std::vector<TComplex> check_buf;
    std::vector<int>      siteCoord;

    // Check sequential insertion of gamma matrix gives same result as 
    // insertion of gamma at sink upon contraction. Assume q uses a point 
    // source.
    
    envGetTmp(LatticeComplex, c);
    siteCoord = strToVec<int>(par().origin);
    peekSite(qSite, qSeq, siteCoord);
    test = trace(g*qSite);

    c = trace(adj(g)*g5*adj(q)*g5*g*q);
    sliceSum(c, check_buf, Tp);
    check = TensorRemove(check_buf[par().t_g]);

    LOG(Message) << "Seq Result = " << abs(test)  << std::endl;
    LOG(Message) << "Reference  = " << abs(check) << std::endl;

    // Check difference = 0
    check -= test;

    LOG(Message) << "Consistency check for sequential " << par().gamma  
                 << " insertion = " << abs(check) << std::endl;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_TestSeqGamma_hpp_
