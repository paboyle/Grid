/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: extras/Hadrons/Modules/MUtilities/TestSeqConserved.hpp

Copyright (C) 2017

Author: Andrew Lawson    <andrew.lawson1991@gmail.com>

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

#ifndef Hadrons_TestSeqConserved_hpp_
#define Hadrons_TestSeqConserved_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/*
  Ward Identity contractions using sequential propagators.
 -----------------------------
 
 * options:
 - q:      point source propagator, 5D if available (string)
 - q4d:    4D point source propagator, duplicate of q if q is 4D (string)
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
                                    std::string,  q4d,
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
    TYPE_ALIASES(FImpl,);
public:
    // constructor
    TTestSeqConserved(const std::string name);
    // destructor
    virtual ~TTestSeqConserved(void) = default;
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_NS(TestSeqConserved, TTestSeqConserved<FIMPL>, MUtilities);

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
    std::vector<std::string> in = {par().q, par().q4d, 
                                   par().qSeq, par().action};
    
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
    
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TTestSeqConserved<FImpl>::execute(void)
{
    PropagatorField tmp(env().getGrid());
    PropagatorField &q    = *env().template getObject<PropagatorField>(par().q);
    PropagatorField &q4d  = *env().template getObject<PropagatorField>(par().q4d);
    PropagatorField &qSeq = *env().template getObject<PropagatorField>(par().qSeq);
    FMat            &act  = *(env().template getObject<FMat>(par().action));
    Gamma           g5(Gamma::Algebra::Gamma5);
    SitePropagator  qSite;
    LatticeComplex  c(env().getGrid());
    Complex         seq_res, check_res;
    std::vector<TComplex> check_buf;

    // Check sequential insertion of current gives same result as conserved 
    // current sink upon contraction. Assume q uses a point source.
    std::vector<int> siteCoord;
    siteCoord = strToVec<int>(par().origin);
    peekSite(qSite, q, siteCoord);
    seq_res = trace(g5*qSite);

    act.ContractConservedCurrent(q, q, tmp, par().curr, par().mu);
    c = trace(tmp);
    sliceSum(c, check_buf, Tp);
    check_res = TensorRemove(check_buf[par().t_J]);

    // Check difference = 0
    check_res -= seq_res;

    LOG(Message) << "Consistency check for sequential conserved " 
                 << par().curr << " current insertion = " << abs(check_res) 
                 << std::endl;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_TestSeqConserved_hpp_
