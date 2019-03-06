/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Exceptions.hpp

Copyright (C) 2015-2019

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

#ifndef Hadrons_Exceptions_hpp_
#define Hadrons_Exceptions_hpp_

#include <stdexcept>
#include <execinfo.h>
#ifndef Hadrons_Global_hpp_
#include <Hadrons/Global.hpp>
#endif

#define HADRONS_SRC_LOC std::string(__FUNCTION__) + " at " \
                        + std::string(__FILE__) + ":" + std::to_string(__LINE__)
#define HADRONS_BACKTRACE_MAX 128
#ifdef HAVE_EXECINFO_H
#define HADRONS_CACHE_BACKTRACE \
{\
    void* _callstack[HADRONS_BACKTRACE_MAX];\
    int _i, _frames = backtrace(_callstack, HADRONS_BACKTRACE_MAX);\
    char** _strs = backtrace_symbols(_callstack, _frames);\
    Grid::Hadrons::Exceptions::backtraceStr.clear();\
    for (_i = 0; _i < _frames; ++_i)\
    {\
        Hadrons::Exceptions::backtraceStr.push_back(std::string(_strs[_i]));\
    }\
    free(_strs);\
}
#else
#define HADRONS_CACHE_BACKTRACE \
Grid::Hadrons::Exceptions::backtraceStr.clear();\
Grid::Hadrons::Exceptions::backtraceStr.push_back("<backtrace not supported>");
#endif

#define HADRONS_ERROR(exc, msg)\
HADRONS_CACHE_BACKTRACE \
throw(Exceptions::exc(msg, HADRONS_SRC_LOC));

#define HADRONS_ERROR_REF(exc, msg, address)\
HADRONS_CACHE_BACKTRACE \
throw(Exceptions::exc(msg, HADRONS_SRC_LOC, address));

#define DECL_EXC(name, base) \
class name: public base\
{\
public:\
    name(std::string msg, std::string loc);\
}

#define DECL_EXC_REF(name, base) \
class name: public base\
{\
public:\
    name(std::string msg, std::string loc, const unsigned int address);\
}

BEGIN_HADRONS_NAMESPACE

namespace Exceptions
{
    // backtrace cache
    extern std::vector<std::string> backtraceStr;

    // logic errors
    DECL_EXC(Logic, std::logic_error);
    DECL_EXC(Definition, Logic);
    DECL_EXC(Implementation, Logic);
    DECL_EXC(Range, Logic);
    DECL_EXC(Size, Logic);

    // runtime errors
    DECL_EXC(Runtime, std::runtime_error);
    DECL_EXC(Argument, Runtime);
    DECL_EXC(Io, Runtime);
    DECL_EXC(Memory, Runtime);
    DECL_EXC(Parsing, Runtime);
    DECL_EXC(Program, Runtime);
    DECL_EXC(System, Runtime);

    // virtual machine errors
    class RuntimeRef: public Runtime
    {
    public:
        RuntimeRef(std::string msg, std::string loc, const unsigned int address)
        : Runtime(msg, loc), address_(address)
        {}
        unsigned int getAddress(void) const
        {
            return address_;
        }
    private:
        unsigned int address_;
    };

    DECL_EXC_REF(ObjectDefinition, RuntimeRef);
    DECL_EXC_REF(ObjectType, RuntimeRef);

    // abort functions
    void abort(const std::exception& e);
}

END_HADRONS_NAMESPACE

#endif // Hadrons_Exceptions_hpp_
