/*
 * Application.hpp, part of Grid
 *
 * Copyright (C) 2015 Antonin Portelli
 *
 * Grid is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Grid is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Grid.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef Hadrons_Application_hpp_
#define Hadrons_Application_hpp_

#include <Hadrons/Global.hpp>
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
                                    ConfigPar, configs);
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
    ModuleFactory                                   &modFactory_;
    std::map<std::string, std::unique_ptr<Module>>  module_;
    std::map<std::string, std::string>              associatedModule_;
    std::map<std::string, std::vector<std::string>> input_;
    std::vector<std::string>                        program_;
};

END_HADRONS_NAMESPACE

#endif // Hadrons_Application_hpp_
