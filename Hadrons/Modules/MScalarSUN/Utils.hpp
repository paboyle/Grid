/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MScalarSUN/Utils.hpp

Copyright (C) 2015-2018

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
#ifndef Hadrons_MScalarSUN_Utils_hpp_
#define Hadrons_MScalarSUN_Utils_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>

BEGIN_HADRONS_NAMESPACE

BEGIN_MODULE_NAMESPACE(MScalarSUN)

GRID_SERIALIZABLE_ENUM(DiffType, undef, forward, 1, backward, 2, central, 3);

template <typename Field>
inline void dmu(Field &out, const Field &in, const unsigned int mu, const DiffType type)
{
    auto & env = Environment::getInstance();

    if (mu >= env.getNd())
    {
        HADRONS_ERROR(Range, "Derivative direction out of range");
    }
    switch(type)
    {
        case DiffType::backward:
            out = in - Cshift(in, mu, -1);
            break;
        case DiffType::forward:
            out = Cshift(in, mu, 1) - in;
            break;
        case DiffType::central:
            out = 0.5*(Cshift(in, mu, 1) - Cshift(in, mu, -1));
            break;
        default:
            HADRONS_ERROR(Argument, "Derivative type invalid");
            break;
    }
}

template <typename Field>
inline void dmuAcc(Field &out, const Field &in, const unsigned int mu, const DiffType type)
{
    auto & env = Environment::getInstance();

    if (mu >= env.getNd())
    {
        HADRONS_ERROR(Range, "Derivative direction out of range");
    }
    switch(type)
    {
        case DiffType::backward:
            out += in - Cshift(in, mu, -1);
            break;
        case DiffType::forward:
            out += Cshift(in, mu, 1) - in;
            break;
        case DiffType::central:
            out += 0.5*(Cshift(in, mu, 1) - Cshift(in, mu, -1));
            break;
        default:
            HADRONS_ERROR(Argument, "Derivative type invalid");
            break;
    }
}

template <class SinkSite, class SourceSite>
std::vector<Complex> makeTwoPoint(const std::vector<SinkSite>   &sink,
                                  const std::vector<SourceSite> &source,
                                  const double factor = 1.)
{
    assert(sink.size() == source.size());
    
    unsigned int         nt = sink.size();
    std::vector<Complex> res(nt, 0.);
    
    for (unsigned int dt = 0; dt < nt; ++dt)
    {
        for (unsigned int t  = 0; t < nt; ++t)
        {
            res[dt] += trace(sink[(t+dt)%nt]*adj(source[t]));
        }
        res[dt] *= factor/static_cast<double>(nt);
    }
    
    return res;
}

inline std::string varName(const std::string name, const std::string suf)
{
    return name + "_" + suf;
}

inline std::string varName(const std::string name, const unsigned int mu)
{
    return varName(name, std::to_string(mu));
}

inline std::string varName(const std::string name, const unsigned int mu, 
                           const unsigned int nu)
{
    return varName(name, std::to_string(mu) + "_" + std::to_string(nu));
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MScalarSUN_Utils_hpp_
