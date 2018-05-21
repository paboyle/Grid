/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: extras/Hadrons/Exceptions.cc

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

#include <Grid/Hadrons/Exceptions.hpp>
#include <Grid/Hadrons/VirtualMachine.hpp>
#include <Grid/Hadrons/Module.hpp>

#ifndef ERR_SUFF
#define ERR_SUFF " (" + loc + ")"
#endif

#define CONST_EXC(name, init) \
name::name(std::string msg, std::string loc)\
:init\
{}

using namespace Grid;
using namespace Hadrons;
using namespace Exceptions;

// logic errors
CONST_EXC(Logic, logic_error(msg + ERR_SUFF))
CONST_EXC(Definition, Logic("definition error: " + msg, loc))
CONST_EXC(Implementation, Logic("implementation error: " + msg, loc))
CONST_EXC(Range, Logic("range error: " + msg, loc))
CONST_EXC(Size, Logic("size error: " + msg, loc))

// runtime errors
CONST_EXC(Runtime, runtime_error(msg + ERR_SUFF))
CONST_EXC(Argument, Runtime("argument error: " + msg, loc))
CONST_EXC(Io, Runtime("IO error: " + msg, loc))
CONST_EXC(Memory, Runtime("memory error: " + msg, loc))
CONST_EXC(Parsing, Runtime("parsing error: " + msg, loc))
CONST_EXC(Program, Runtime("program error: " + msg, loc))
CONST_EXC(System, Runtime("system error: " + msg, loc))

// abort functions
void Grid::Hadrons::Exceptions::abort(const std::exception& e)
{
    auto &vm = VirtualMachine::getInstance();
    int  mod = vm.getCurrentModule();

    LOG(Error) << "FATAL ERROR -- Exception " << typeName(&typeid(e)) 
               << std::endl;
    if (mod >= 0)
    {
        LOG(Error) << "During execution of module '"
                    << vm.getModuleName(mod) << "' (address " << mod << ")"
                    << std::endl;
    }
    LOG(Error) << e.what() << std::endl;
    LOG(Error) << "Aborting program" << std::endl;
    Grid_finalize();

    exit(EXIT_FAILURE);
}
