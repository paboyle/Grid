/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Utilities/HadronsXmlRun.cc

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

#include <Hadrons/Application.hpp>

using namespace Grid;
using namespace QCD;
using namespace Hadrons;

int main(int argc, char *argv[])
{
    // parse command line
    std::string parameterFileName, scheduleFileName = "";
    
    if (argc < 2)
    {
        std::cerr << "usage: " << argv[0] << " <parameter file> [<precomputed schedule>] [Grid options]";
        std::cerr << std::endl;
        std::exit(EXIT_FAILURE);
    }
    parameterFileName = argv[1];
    if (argc > 2)
    {
        if (argv[2][0] != '-')
        {
            scheduleFileName = argv[2];
        }
    }
    
    // initialization
    Grid_init(&argc, &argv);
    
    // execution
    try
    {
        Application application(parameterFileName);
        
        application.parseParameterFile(parameterFileName);
        if (!scheduleFileName.empty())
        {
            application.loadSchedule(scheduleFileName);
        }
        application.run();
    }
    catch (const std::exception& e)
    {
        Exceptions::abort(e);
    }
    
    // epilogue
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();
    
    return EXIT_SUCCESS;
}
