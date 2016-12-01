/*******************************************************************************
Grid physics library, www.github.com/paboyle/Grid 

Source file: programs/Hadrons/Global.cc

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

#include <Grid/Hadrons/Global.hpp>

using namespace Grid;
using namespace QCD;
using namespace Hadrons;

HadronsLogger Hadrons::HadronsLogError(1,"Error");
HadronsLogger Hadrons::HadronsLogWarning(1,"Warning");
HadronsLogger Hadrons::HadronsLogMessage(1,"Message");
HadronsLogger Hadrons::HadronsLogIterative(1,"Iterative");
HadronsLogger Hadrons::HadronsLogDebug(1,"Debug");

// pretty size formatting //////////////////////////////////////////////////////
std::string Hadrons::sizeString(long unsigned int bytes)

{
    constexpr unsigned int bufSize = 256;
    const char             *suffixes[7] = {"", "K", "M", "G", "T", "P", "E"};
    char                   buf[256];
    long unsigned int      s     = 0;
    double                 count = bytes;
    
    while (count >= 1024 && s < 7)
    {
        s++;
        count /= 1024;
    }
    if (count - floor(count) == 0.0)
    {
        snprintf(buf, bufSize, "%d %sB", (int)count, suffixes[s]);
    }
    else
    {
        snprintf(buf, bufSize, "%.1f %sB", count, suffixes[s]);
    }
    
    return std::string(buf);
}
