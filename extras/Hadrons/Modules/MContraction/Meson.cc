/*******************************************************************************
Grid physics library, www.github.com/paboyle/Grid 

Source file: programs/Hadrons/Meson.cc

Copyright (C) 2015

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

See the full license in the file "LICENSE" in the top level distribution 
directory.
*******************************************************************************/

#include <Grid/Hadrons/Modules/MContraction/Meson.hpp>

using namespace Grid;
using namespace QCD;
using namespace Hadrons;
using namespace MContraction;

/******************************************************************************
 *                           Meson implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
Meson::Meson(const std::string name)
: Module<MesonPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
std::vector<std::string> Meson::getInput(void)
{
    std::vector<std::string> input = {par().q1, par().q2};
    
    return input;
}

std::vector<std::string> Meson::getOutput(void)
{
    std::vector<std::string> output = {getName()};
    
    return output;
}

// execution ///////////////////////////////////////////////////////////////////
void Meson::execute(void)
{
    LOG(Message) << "Computing meson contraction '" << getName() << "' using"
                 << " quarks '" << par().q1 << "' and '" << par().q2 << "'"
                 << std::endl;
    
    XmlWriter             writer(par().output);
    PropagatorField       &q1 = *env().getObject<PropagatorField>(par().q1);
    PropagatorField       &q2 = *env().getObject<PropagatorField>(par().q2);
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
