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
#include <Hadrons/Environment.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

class TrajRange: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(TrajRange,
                                    unsigned int, start,
                                    unsigned int, end,
                                    unsigned int, step);
};

class GlobalPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(GlobalPar,
                                    TrajRange,   trajCounter,
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
    // loop on configurations
    void configLoop(void);
private:
    long unsigned int                               locVol_;
    std::string                                     parameterFileName_;
    GlobalPar                                       par_;
    Environment                                     &env_;
    std::vector<unsigned int>                       program_;
};

END_HADRONS_NAMESPACE

#endif // Hadrons_Application_hpp_
