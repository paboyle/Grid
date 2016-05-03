/*******************************************************************************
Grid physics library, www.github.com/paboyle/Grid 

Source file: programs/Hadrons/MSource.cc

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

See the full license in the file "LICENSE" in the top level distribution 
directory.
*******************************************************************************/

#include <Hadrons/MSource.hpp>

#define ERROR_SUF " (source '" << getName() << "')"

using namespace Grid;
using namespace QCD;
using namespace Hadrons;

/******************************************************************************
*                         MSource implementation                              *
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
MSource::MSource(const std::string name)
: Module(name)
{}

// parse parameters ////////////////////////////////////////////////////////////
void MSource::parseParameters(XmlReader &reader, const std::string name)
{
    read(reader, name, par_);
}

// dependencies/products ///////////////////////////////////////////////////////
std::vector<std::string> MSource::getInput(void)
{
    return std::vector<std::string>();
}

std::vector<std::string> MSource::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// allocation //////////////////////////////////////////////////////////////////
void MSource::allocate(Environment &env)
{
    switch (par_.sourceType)
    {
        // 4D sources
        case Grid::SourceType::point:
        case Grid::SourceType::z2Band:
            env.createProp(getName());
            src_ = env.getProp(getName());
            break;
        // error
        default:
            HADRON_ERROR("no allocation implemented for source type '"
                         << par_.sourceType << "'" << ERROR_SUF);
            break;
    }
}

// execution
#define ARG_CHECK(n)\
if (par_.arguments.size() != (n))\
{\
    HADRON_ERROR("source type '" << par_.sourceType << "' expect "\
                 << (n) << " arguments (got "\
                 << par_.arguments.size() << ")" << ERROR_SUF);\
}

void MSource::execute(Environment &env)
{
    LOG(Message) << "generating source '" << getName() << "' of type '"
                 << par_.sourceType << "'" << std::endl;
    switch (par_.sourceType)
    {
        // point source
        case Grid::SourceType::point:
        {
            ARG_CHECK(1);
            
            std::vector<int> origin = strToVec<int>(par_.arguments[0]);
            SpinColourMatrix id(1.);
            
            if (origin.size() != Nd)
            {
                HADRON_ERROR("point source origin dimension different from "
                             << Nd << ERROR_SUF);
            }
            *src_ = zero;
            pokeSite(id, *src_, origin);
            
            break;
        }
        // z2Band source
        case Grid::SourceType::z2Band:
        {
            ARG_CHECK(2);
            
            int                        ta = std::stoi(par_.arguments[0]);
            int                        tb = std::stoi(par_.arguments[1]);
            Lattice<iScalar<vInteger>> t(env.getGrid());
            LatticeComplex             eta(env.getGrid());
            LatticeFermion             phi(env.getGrid());
            Complex                    shift(1., 1.);
            
            LatticeCoordinate(t, Tp);
            bernoulli(*env.get4dRng(), eta);
            eta   = (2.*eta - shift)*(1./::sqrt(2.));
            eta   = where((t >= ta) and (t <= tb), eta, 0.*eta);
            *src_ = 1.;
            *src_ = (*src_)*eta;
            
            break;
        }
        // error
        default:
        {
            HADRON_ERROR("no definition implemented for source type '"
                         << par_.sourceType << "'" << ERROR_SUF);
            break;
        }
    }
}
