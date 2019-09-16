/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Utilities/HadronsXmlValidate.cc

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

#include <Hadrons/Application.hpp>
#include <sys/stat.h>
#include <fstream>
#include <string>

using namespace Grid;
using namespace Hadrons;

// Does the specified file exist?
bool FileExists(const std::string& Filename)
{
  struct stat buf;
  return stat(Filename.c_str(), &buf) != -1;
}

void Shorten( Application &app, const std::string &FileList, const std::string OutFileName )
{
  std::vector<std::string> Except;
  std::ifstream list{ FileList };
  for( std::string s; std::getline(list, s); ) {
    //const std::string::size_type l{ s.find_first_of( '.' ) };
    //if( l != std::string::npos )
      //s.resize( l );
    if( s.length() )
      Except.push_back( s );
  }
  std::sort( Except.begin(), Except.end() );
  for( const std::string &s : Except )
    std::cout << s << std::endl;
  app.saveParameterFile( OutFileName, Except );
}

int main(int argc, char *argv[])
{
    // parse command line
    std::string parameterFileName;
    
    if (argc != 2 && argc != 4)
    {
        std::cerr << "usage: " << argv[0] << " <parameter file> [filelist.txt output.xml]";
        std::cerr << std::endl;
        std::exit(EXIT_FAILURE);
    }
    parameterFileName = argv[1];
    if( argc == 4 )
        Grid_init(&argc, &argv);
    try
    {
        Application application(parameterFileName);
        
        application.parseParameterFile(parameterFileName);
        auto &vm = VirtualMachine::getInstance();
        vm.getModuleGraph();
        LOG(Message) << "Application valid (check XML warnings though)" 
                     << std::endl;
      if( argc == 4 ) {
        const std::string FileList{ argv[3] };
        const std::string OutFileName{ argv[2] };
        if( !FileExists( FileList ) )
          std::cout << "File list \"" << FileList << "\" does not exist" << std::endl;
        else {
          Shorten( application, FileList, OutFileName );
        }
      }
    }
    catch (const std::exception& e)
    {
        Exceptions::abort(e);
    }
    if( argc == 4 )
        Grid_finalize();

    return EXIT_SUCCESS;
}
