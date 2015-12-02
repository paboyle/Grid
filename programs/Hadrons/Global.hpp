/*
 * Globals.hpp, part of Grid
 *
 * Copyright (C) 2015 Antonin Portelli
 *
 * LatCore is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * LatCore is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with LatCore.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef Hadrons_Global_hpp_
#define Hadrons_Global_hpp_

#include <Grid.h>
#include <set>
#include <stack>

#define BEGIN_HADRONS_NAMESPACE \
namespace Hadrons {\
using namespace Grid;
#define END_HADRONS_NAMESPACE }

BEGIN_HADRONS_NAMESPACE

class HadronsLogger: public Logger
{
public:
    HadronsLogger(int on, std::string nm): Logger("Hadrons", on, nm){};
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

END_HADRONS_NAMESPACE

#endif // Hadrons_Global_hpp_
