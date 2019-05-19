/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MContraction/WeakEye3pt.hpp

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
#ifndef Hadrons_MContraction_WeakEye3pt_hpp_
#define Hadrons_MContraction_WeakEye3pt_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/*
 * Weak Hamiltonian meson 3-pt diagrams, eye topologies.
 * 
 * Schematics:       loop                 |                  
 *                  /-<-¬                 |                             
 *                 /     \                |            qbl     G     qbr
 *                 \     /                |        /----<------*------<----¬         
 *            qbl   \   /    qbr          |       /          /-*-¬          \
 *       /-----<-----* *-----<----¬       |      /          /  G  \          \
 *  gIn *            G G           * gOut | gIn *           \     /  loop    * gOut
 *       \                        /       |      \           \->-/          /   
 *        \                      /        |       \                        /       
 *         \---------->---------/         |        \----------->----------/        
 *                   qs                   |                   qs                  
 *                                        |
 *                one trace               |                two traces
 * 
 * one trace : tr(qbr*gOut*qs*adj(gIn)*g5*adj(qbl)*g5*G*loop*G)
 * two traces: tr(qbr*gOut*qs*adj(gIn)*g5*adj(qbl)*g5*G)*tr(loop*G)
 * 
 */

BEGIN_MODULE_NAMESPACE(MContraction)

class WeakEye3ptPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(WeakEye3ptPar,
                                    std::string,    qBarLeft,
                                    std::string,    qBarRight,
                                    std::string,    qSpectator,
                                    std::string,    loop,
                                    unsigned int,   tOut,
                                    Gamma::Algebra, gammaIn,
                                    Gamma::Algebra, gammaOut,
                                    std::string,    output);
};

template <typename FImpl>
class TWeakEye3pt: public Module<WeakEye3ptPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    class Metadata: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Metadata,
                                        Gamma::Algebra, in,
                                        Gamma::Algebra, out,
                                        Gamma::Algebra, op,
                                        unsigned int,   trace);
    };
    typedef Correlator<Metadata> Result;
public:
    // constructor
    TWeakEye3pt(const std::string name);
    // destructor
    virtual ~TWeakEye3pt(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(WeakEye3pt, TWeakEye3pt<FIMPL>, MContraction);

/******************************************************************************
 *                        TWeakEye3pt implementation                          *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TWeakEye3pt<FImpl>::TWeakEye3pt(const std::string name)
: Module<WeakEye3ptPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TWeakEye3pt<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().qBarLeft, par().qBarRight, 
                                   par().qSpectator, par().loop};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TWeakEye3pt<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TWeakEye3pt<FImpl>::setup(void)
{
    envTmpLat(ComplexField, "corr");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TWeakEye3pt<FImpl>::execute(void)
{
    LOG(Message) << "Computing mesonic weak 3pt contractions, eye topologies" << std::endl;
    LOG(Message) << "gIn : " << par().gammaIn << std::endl;
    LOG(Message) << "gOut: " << par().gammaIn << std::endl;
    LOG(Message) << "tOut: " << par().tOut << std::endl;
    LOG(Message) << "qbl : " << par().qBarLeft << std::endl;
    LOG(Message) << "qbr : " << par().qBarRight << std::endl;
    LOG(Message) << "qs  : " << par().qSpectator << std::endl;
    LOG(Message) << "loop: " << par().loop << std::endl;

    std::vector<Result> result;
    Result              r;
    auto                &qbl  = envGet(PropagatorField, par().qBarLeft);
    auto                &qbr  = envGet(PropagatorField, par().qBarRight);
    auto                &loop = envGet(PropagatorField, par().loop);
    auto                &qs   = envGet(SlicedPropagator, par().qSpectator);
    auto                qst   = qs[par().tOut];
    Gamma               gIn(par().gammaIn), gOut(par().gammaOut);
    Gamma               g5(Gamma::Algebra::Gamma5);

    envGetTmp(ComplexField, corr);
    r.info.in  = par().gammaIn;
    r.info.out = par().gammaOut;
    for (auto &G: Gamma::gall)
    {
        SlicedComplex buf;

        r.info.op = G.g;
        // one trace
        corr = trace(qbr*gOut*qst*adj(gIn)*g5*adj(qbl)*g5*G*loop*G);
        sliceSum(corr, buf, Tp);
        r.corr.clear();
        for (unsigned int t = 0; t < buf.size(); ++t)
        {
            r.corr.push_back(TensorRemove(buf[t]));
        }
        r.info.trace = 1;
        result.push_back(r);
        // two traces
        corr = trace(qbr*gOut*qst*adj(gIn)*g5*adj(qbl)*g5*G)*trace(loop*G);
        sliceSum(corr, buf, Tp);
        r.corr.clear();
        for (unsigned int t = 0; t < buf.size(); ++t)
        {
            r.corr.push_back(TensorRemove(buf[t]));
        }
        r.info.trace = 2;
        result.push_back(r);
    }
    saveResult(par().output, "weakEye3pt", result);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_WeakEye3pt_hpp_
