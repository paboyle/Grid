/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: extras/Hadrons/Global.cc

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

#include <Grid/Hadrons/Global.hpp>

using namespace Grid;
using namespace QCD;
using namespace Hadrons;

HadronsLogger Hadrons::HadronsLogError(1,"Error");
HadronsLogger Hadrons::HadronsLogWarning(1,"Warning");
HadronsLogger Hadrons::HadronsLogMessage(1,"Message");
HadronsLogger Hadrons::HadronsLogIterative(1,"Iterative");
HadronsLogger Hadrons::HadronsLogDebug(1,"Debug");
HadronsLogger Hadrons::HadronsLogIRL(1,"IRL");

void Hadrons::initLogger(void)
{
    auto w = std::string("Hadrons").length();
    GridLogError.setTopWidth(w);
    GridLogWarning.setTopWidth(w);
    GridLogMessage.setTopWidth(w);
    GridLogIterative.setTopWidth(w);
    GridLogDebug.setTopWidth(w);
    GridLogIRL.setTopWidth(w);
    HadronsLogError.Active(GridLogError.isActive());
    HadronsLogWarning.Active(GridLogWarning.isActive());
    HadronsLogMessage.Active(GridLogMessage.isActive());
    HadronsLogIterative.Active(GridLogIterative.isActive());
    HadronsLogDebug.Active(GridLogDebug.isActive());
    HadronsLogIRL.Active(GridLogIRL.isActive());
}

// type utilities //////////////////////////////////////////////////////////////
constexpr unsigned int maxNameSize = 1024u;

std::string Hadrons::typeName(const std::type_info *info)
{
    char        *buf;
    std::string name;
    
    buf  = abi::__cxa_demangle(info->name(), nullptr, nullptr, nullptr);
    name = buf;
    free(buf);
    
    return name;
}

// default writers/readers /////////////////////////////////////////////////////
#ifdef HAVE_HDF5
const std::string Hadrons::resultFileExt = "h5";
#else
const std::string Hadrons::resultFileExt = "xml";
#endif
