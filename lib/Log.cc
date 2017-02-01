/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/Log.cc

Copyright (C) 2015

Author: Antonin Portelli <antonin.portelli@me.com>
Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: paboyle <paboyle@ph.ed.ac.uk>

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
directory
*************************************************************************************/
/*  END LEGAL */
#include <Grid/Grid.h>

#include <cxxabi.h>

namespace Grid {

  std::string demangle(const char* name) {
    
    int status = -4; // some arbitrary value to eliminate the compiler warning
    
    // enable c++11 by passing the flag -std=c++11 to g++
    std::unique_ptr<char, void(*)(void*)> res {
      abi::__cxa_demangle(name, NULL, NULL, &status),
	std::free
	};
    
    return (status==0) ? res.get() : name ;
  }
  
GridStopWatch Logger::StopWatch;
int Logger::timestamp;
std::ostream Logger::devnull(0);

void GridLogTimestamp(int on){
  Logger::Timestamp(on);
}

Colours GridLogColours(0);
GridLogger GridLogError(1, "Error", GridLogColours, "RED");
GridLogger GridLogWarning(1, "Warning", GridLogColours, "YELLOW");
GridLogger GridLogMessage(1, "Message", GridLogColours, "NORMAL");
GridLogger GridLogDebug(1, "Debug", GridLogColours, "PURPLE");
GridLogger GridLogPerformance(1, "Performance", GridLogColours, "GREEN");
GridLogger GridLogIterative(1, "Iterative", GridLogColours, "BLUE");
GridLogger GridLogIntegrator(1, "Integrator", GridLogColours, "BLUE");

void GridLogConfigure(std::vector<std::string> &logstreams) {
  GridLogError.Active(0);
  GridLogWarning.Active(0);
  GridLogMessage.Active(1); // at least the messages should be always on
  GridLogIterative.Active(0);
  GridLogDebug.Active(0);
  GridLogPerformance.Active(0);
  GridLogIntegrator.Active(0);
  GridLogColours.Active(0);

  for (int i = 0; i < logstreams.size(); i++) {
    if (logstreams[i] == std::string("Error")) GridLogError.Active(1);
    if (logstreams[i] == std::string("Warning")) GridLogWarning.Active(1);
    if (logstreams[i] == std::string("NoMessage")) GridLogMessage.Active(0);
    if (logstreams[i] == std::string("Iterative")) GridLogIterative.Active(1);
    if (logstreams[i] == std::string("Debug")) GridLogDebug.Active(1);
    if (logstreams[i] == std::string("Performance"))
      GridLogPerformance.Active(1);
    if (logstreams[i] == std::string("Integrator")) GridLogIntegrator.Active(1);
    if (logstreams[i] == std::string("Colours")) GridLogColours.Active(1);
  }
}

////////////////////////////////////////////////////////////
// Verbose limiter on MPI tasks
////////////////////////////////////////////////////////////
void Grid_quiesce_nodes(void) {
  int me = 0;
#if defined(GRID_COMMS_MPI) || defined(GRID_COMMS_MPI3) || defined(GRID_COMMS_MPI3L)
  MPI_Comm_rank(MPI_COMM_WORLD, &me);
#endif
#ifdef GRID_COMMS_SHMEM
  me = shmem_my_pe();
#endif
  if (me) {
    std::cout.setstate(std::ios::badbit);
  }
}

void Grid_unquiesce_nodes(void) {
#ifdef GRID_COMMS_MPI
  std::cout.clear();
#endif
}
}
