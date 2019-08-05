/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MContraction/Baryon.hpp

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

#ifndef Hadrons_MContraction_SelfContract_hpp_
#define Hadrons_MContraction_SelfContract_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                                SelfContract                                    *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class SelfContractPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(SelfContractPar,
                                    std::string,    q_self,
                                    Gamma::Algebra, gamma,
                                    std::string,    sink,
                                    std::string,    output);
};

template <typename FImpl>
class TSelfContract: public Module<SelfContractPar>
{
    FERM_TYPE_ALIASES(FImpl,);
    BASIC_TYPE_ALIASES(ScalarImplCR, Scalar);
    SINK_TYPE_ALIASES(Scalar);
    class Result: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
                                        Gamma::Algebra, gamma,
                                        std::vector<Complex>, corr);
    };
public:
    // constructor
    TSelfContract(const std::string name);
    // destructor
    virtual ~TSelfContract(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(SelfContract, TSelfContract<FIMPL>, MContraction);

/******************************************************************************
 *                       TSelfContract implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TSelfContract<FImpl>::TSelfContract(const std::string name)
: Module<SelfContractPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TSelfContract<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().q_self, par().sink};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TSelfContract<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TSelfContract<FImpl>::setup(void)
{
    envTmpLat(LatticeComplex, "c");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TSelfContract<FImpl>::execute(void)
{
    LOG(Message) << "Computing self-contracted loop '" << getName() 
                 << "' using '" << par().q_self << "' with " << par().gamma 
                 << " insertion and sink " << par().sink << std::endl;

    auto                  &q_self = envGet(PropagatorField, par().q_self);
    Gamma                 gamma(par().gamma);
    std::vector<TComplex> buf;
    Result                result;
    envGetTmp(LatticeComplex, c);
    SinkFnScalar &sink = envGet(SinkFnScalar, par().sink);

    c = trace(gamma*q_self);
    buf = sink(c);
    result.gamma = par().gamma;
    result.corr.resize(buf.size());
    for (unsigned int t = 0; t < buf.size(); ++t)
    {
        result.corr[t] = TensorRemove(buf[t]);
    }
    saveResult(par().output, "selfContr", result);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_SelfContract_hpp_
