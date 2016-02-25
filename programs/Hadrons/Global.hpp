/*******************************************************************************
Grid physics library, www.github.com/paboyle/Grid 

Source file: programs/Hadrons/Global.hpp

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

#ifndef Hadrons_Global_hpp_
#define Hadrons_Global_hpp_

#include <Grid.h>
#include <set>
#include <stack>

#define BEGIN_HADRONS_NAMESPACE \
namespace Hadrons {\
using namespace Grid;\
using namespace QCD;
#define END_HADRONS_NAMESPACE }

BEGIN_HADRONS_NAMESPACE

class HadronsLogger: public Logger
{
public:
    HadronsLogger(int on, std::string nm): Logger("Hadrons", on, nm,
                                                  Logger::BLACK){};
};

#define LOG(channel) std::cout << HadronsLog##channel
#define HADRON_ERROR(msg)\
LOG(Error) << msg << std::endl;\
abort();

#define DEBUG_VAR(var) LOG(Debug) << #var << "= " << (var) << std::endl;

extern HadronsLogger HadronsLogError;
extern HadronsLogger HadronsLogWarning;
extern HadronsLogger HadronsLogMessage;
extern HadronsLogger HadronsLogDebug;

// singleton pattern
#define SINGLETON(name)\
public:\
    name(const name &e) = delete;\
    void operator=(const name &e) = delete;\
    static name & getInstance(void)\
    {\
        static name e;\
        return e;\
    }\
private:\
    name(void);

END_HADRONS_NAMESPACE

#endif // Hadrons_Global_hpp_
