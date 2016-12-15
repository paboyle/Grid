/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: extras/Hadrons/Modules/MContraction/Meson.hpp

Copyright (C) 2015
Copyright (C) 2016

Author: Antonin Portelli <antonin.portelli@me.com>

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

#ifndef Hadrons_Meson_hpp_
#define Hadrons_Meson_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                                TMeson                                       *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class MesonPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(MesonPar,
                                    std::string, q1,
                                    std::string, q2,
                                    std::string, output);
};

template <typename FImpl1, typename FImpl2>
class TMeson: public Module<MesonPar>
{
public:
    TYPE_ALIASES(FImpl1, 1);
    TYPE_ALIASES(FImpl2, 2);
    class Result: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
                                        std::vector<std::vector<std::vector<Complex>>>, corr);
    };
public:
    // constructor
    TMeson(const std::string name);
    // destructor
    virtual ~TMeson(void) = default;
    // dependencies/products
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_NS(Meson, ARG(TMeson<FIMPL, FIMPL>), MContraction);

/******************************************************************************
 *                           TMeson implementation                            *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
TMeson<FImpl1, FImpl2>::TMeson(const std::string name)
: Module<MesonPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
std::vector<std::string> TMeson<FImpl1, FImpl2>::getInput(void)
{
    std::vector<std::string> input = {par().q1, par().q2};
    
    return input;
}

template <typename FImpl1, typename FImpl2>
std::vector<std::string> TMeson<FImpl1, FImpl2>::getOutput(void)
{
    std::vector<std::string> output = {getName()};
    
    return output;
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
void TMeson<FImpl1, FImpl2>::execute(void)
{
    LOG(Message) << "Computing meson contractions '" << getName() << "' using"
                 << " quarks '" << par().q1 << "' and '" << par().q2 << "'"
                 << std::endl;
    
    XmlWriter             writer(par().output);
    PropagatorField1      &q1 = *env().template getObject<PropagatorField1>(par().q1);
    PropagatorField2      &q2 = *env().template getObject<PropagatorField2>(par().q2);
    LatticeComplex        c(env().getGrid());
    SpinMatrix            g[Ns*Ns], g5;
    std::vector<TComplex> buf;
    Result                result;
    
    g5 = makeGammaProd(Ns*Ns - 1);
    result.corr.resize(Ns*Ns);
    for (unsigned int i = 0; i < Ns*Ns; ++i)
    {
        g[i] = makeGammaProd(i);
    }
    for (unsigned int iSink = 0; iSink < Ns*Ns; ++iSink)
    {
        result.corr[iSink].resize(Ns*Ns);
        for (unsigned int iSrc = 0; iSrc < Ns*Ns; ++iSrc)
        {
            c = trace(g[iSink]*q1*g[iSrc]*g5*adj(q2)*g5);
            sliceSum(c, buf, Tp);
            result.corr[iSink][iSrc].resize(buf.size());
            for (unsigned int t = 0; t < buf.size(); ++t)
            {
                result.corr[iSink][iSrc][t] = TensorRemove(buf[t]);
            }
        }
    }
    write(writer, "meson", result);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_Meson_hpp_
