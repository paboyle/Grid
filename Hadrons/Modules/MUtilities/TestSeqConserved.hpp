/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MUtilities/TestSeqConserved.hpp

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

#ifndef Hadrons_MUtilities_TestSeqConserved_hpp_
#define Hadrons_MUtilities_TestSeqConserved_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/*
  Ward Identity contractions using sequential propagators.
 -----------------------------
 
 * options:
 - q:      point source propagator, 5D if available (string)
 - qSeq:   result of sequential insertion of conserved current using q (string)
 - action: action used for computation of q (string)
 - origin: string giving point source origin of q (string)
 - t_J:    time at which sequential current is inserted (int)
 - mu:     Lorentz index of current inserted (int)
 - curr:   current type, e.g. vector/axial (Current)
*/

/******************************************************************************
 *                            TestSeqConserved                                *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)

class TestSeqConservedPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(TestSeqConservedPar,
                                    std::string,  q,
                                    std::string,  qSeq,
                                    std::string,  action,
                                    std::string,  origin,
                                    unsigned int, t_J,
                                    unsigned int, mu,
                                    Current,      curr);
};

template <typename FImpl>
class TTestSeqConserved: public Module<TestSeqConservedPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TTestSeqConserved(const std::string name);
    // destructor
    virtual ~TTestSeqConserved(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(TestSeqConserved, TTestSeqConserved<FIMPL>, MUtilities);

/******************************************************************************
 *                     TTestSeqConserved implementation                       *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TTestSeqConserved<FImpl>::TTestSeqConserved(const std::string name)
: Module<TestSeqConservedPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TTestSeqConserved<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().q, par().qSeq, par().action};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TTestSeqConserved<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TTestSeqConserved<FImpl>::setup(void)
{
    auto Ls = env().getObjectLs(par().q);
    if (Ls != env().getObjectLs(par().action))
    {
        HADRONS_ERROR(Size, "Ls mismatch between quark action and propagator");
    }
    envTmpLat(PropagatorField, "tmp");
    envTmpLat(LatticeComplex, "c");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TTestSeqConserved<FImpl>::execute(void)
{
    // Check sequential insertion of current gives same result as conserved 
    // current sink upon contraction. Assume q uses a point source.

    auto                  &q    = envGet(PropagatorField, par().q);
    auto                  &qSeq = envGet(PropagatorField, par().qSeq);
    auto                  &act  = envGet(FMat, par().action);
    Gamma                 g5(Gamma::Algebra::Gamma5);
    Gamma::Algebra        gA = (par().curr == Current::Axial) ?
                                  Gamma::Algebra::Gamma5 :
                                  Gamma::Algebra::Identity;
    Gamma                 g(gA);
    SitePropagator        qSite;
    Complex               test_S, test_V, check_S, check_V;
    std::vector<TComplex> check_buf;
    std::vector<int>      siteCoord;

    envGetTmp(PropagatorField, tmp);
    envGetTmp(LatticeComplex, c);
    siteCoord = strToVec<int>(par().origin);
    peekSite(qSite, qSeq, siteCoord);
    test_S = trace(qSite*g);
    test_V = trace(qSite*g*Gamma::gmu[par().mu]);
    act.ContractConservedCurrent(q, q, tmp, par().curr, par().mu);
    c = trace(tmp*g);
    sliceSum(c, check_buf, Tp);
    check_S = TensorRemove(check_buf[par().t_J]);

    c = trace(tmp*g*Gamma::gmu[par().mu]);
    sliceSum(c, check_buf, Tp);
    check_V = TensorRemove(check_buf[par().t_J]);

    LOG(Message) << "Test S  = " << abs(test_S)   << std::endl;
    LOG(Message) << "Test V  = " << abs(test_V) << std::endl;
    LOG(Message) << "Check S = " << abs(check_S) << std::endl;
    LOG(Message) << "Check V = " << abs(check_V) << std::endl;

    // Check difference = 0
    check_S -= test_S;
    check_V -= test_V;

    LOG(Message) << "Consistency check for sequential conserved " 
                 << par().curr << " current insertion: " << std::endl; 
    LOG(Message) << "Diff S  = " << abs(check_S) << std::endl;
    LOG(Message) << "Diff V  = " << abs(check_V) << std::endl;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_TestSeqConserved_hpp_
