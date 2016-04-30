/*******************************************************************************
Grid physics library, www.github.com/paboyle/Grid 

Source file: programs/Hadrons/Application.hpp

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

#ifndef Hadrons_Application_hpp_
#define Hadrons_Application_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/FermionActionFactory.hpp>
#include <Hadrons/Environment.hpp>
#include <Hadrons/ModuleFactory.hpp>

namespace Grid {
    GRID_SERIALIZABLE_ENUM(ConfigType, undef, load, 1, unit, 2, gen, 3);
}

BEGIN_HADRONS_NAMESPACE

class TrajRange: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(TrajRange,
                                    unsigned int, start,
                                    unsigned int, end,
                                    unsigned int, step);
};

class ConfigPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(ConfigPar,
                                    std::string, ioStem,
                                    TrajRange,   range);
};

class GlobalPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(GlobalPar,
                                    ConfigPar,   configs,
                                    std::string, seed);
};

/******************************************************************************
 *                         Main program manager                               *
 ******************************************************************************/
class Application
{
public:

public:
    // constructor
    Application(const std::string parameterFileName);
    // destructor
    virtual ~Application(void);
    // execute
    void run(void);
private:
    // parse parameter file
    void parseParameterFile(void);
    // schedule computation
    void schedule(void);
    // program execution
    void         configLoop(void);
    unsigned int execute(const std::vector<std::string> &program);
private:
    std::string                                     parameterFileName_;
    GlobalPar                                       par_;
    Environment                                     &env_;
    FermionActionFactory                            &actionFactory_;
    ModuleFactory                                   &modFactory_;
    std::map<std::string, std::unique_ptr<Module>>  module_;
    std::map<std::string, std::string>              associatedModule_;
    std::map<std::string, std::vector<std::string>> input_;
    std::vector<std::string>                        program_;
    std::vector<std::vector<std::string>>           freeProg_;
};

END_HADRONS_NAMESPACE

#endif // Hadrons_Application_hpp_
