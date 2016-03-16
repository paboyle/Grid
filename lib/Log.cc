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

    See the full license in the file "LICENSE" in the top level distribution directory
    *************************************************************************************/
    /*  END LEGAL */
#include <Grid.h>

namespace Grid {

GridStopWatch Logger::StopWatch;
std::ostream  Logger::devnull(0);
std::string Logger::BLACK("\033[30m");
std::string Logger::RED("\033[31m");
std::string Logger::GREEN("\033[32m");
std::string Logger::YELLOW("\033[33m");
std::string Logger::BLUE("\033[34m");
std::string Logger::PURPLE("\033[35m");
std::string Logger::CYAN("\033[36m");
std::string Logger::WHITE("\033[37m");
std::string Logger::NORMAL("\033[0;39m");
std::string EMPTY("");

#if 0  
  GridLogger GridLogError      (1,"Error",Logger::RED);
  GridLogger GridLogWarning    (1,"Warning",Logger::YELLOW);
  GridLogger GridLogMessage    (1,"Message",Logger::BLACK);
  GridLogger GridLogDebug      (1,"Debug",Logger::PURPLE);
  GridLogger GridLogPerformance(1,"Performance",Logger::GREEN);
  GridLogger GridLogIterative  (1,"Iterative",Logger::BLUE);
  GridLogger GridLogIntegrator (1,"Integrator",Logger::BLUE);
#else
  GridLogger GridLogError      (1,"Error",EMPTY);
  GridLogger GridLogWarning    (1,"Warning",EMPTY);
  GridLogger GridLogMessage    (1,"Message",EMPTY);
  GridLogger GridLogDebug      (1,"Debug",EMPTY);
  GridLogger GridLogPerformance(1,"Performance",EMPTY);
  GridLogger GridLogIterative  (1,"Iterative",EMPTY);
  GridLogger GridLogIntegrator (1,"Integrator",EMPTY);
#endif

void GridLogConfigure(std::vector<std::string> &logstreams)
{
  GridLogError.Active(0);
  GridLogWarning.Active(0);
  GridLogMessage.Active(0);
  GridLogIterative.Active(0);
  GridLogDebug.Active(0);
  GridLogPerformance.Active(0);
  GridLogIntegrator.Active(0);

  int blackAndWhite = 1;
  if(blackAndWhite){
    Logger::BLACK = std::string("");
    Logger::RED    =Logger::BLACK;
    Logger::GREEN  =Logger::BLACK;
    Logger::YELLOW =Logger::BLACK;
    Logger::BLUE   =Logger::BLACK;
    Logger::PURPLE =Logger::BLACK;
    Logger::CYAN   =Logger::BLACK;
    Logger::WHITE  =Logger::BLACK;
    Logger::NORMAL =Logger::BLACK;
  }

  for(int i=0;i<logstreams.size();i++){
    if ( logstreams[i]== std::string("Error")       ) GridLogError.Active(1);
    if ( logstreams[i]== std::string("Warning")     ) GridLogWarning.Active(1);
    if ( logstreams[i]== std::string("Message")     ) GridLogMessage.Active(1);
    if ( logstreams[i]== std::string("Iterative")   ) GridLogIterative.Active(1);
    if ( logstreams[i]== std::string("Debug")       ) GridLogDebug.Active(1);
    if ( logstreams[i]== std::string("Performance") ) GridLogPerformance.Active(1);
    if ( logstreams[i]== std::string("Integrator" ) ) GridLogIntegrator.Active(1);
  }
}

////////////////////////////////////////////////////////////
// Verbose limiter on MPI tasks
////////////////////////////////////////////////////////////
void Grid_quiesce_nodes(void)
{
  int me=0;
#ifdef GRID_COMMS_MPI
  MPI_Comm_rank(MPI_COMM_WORLD,&me);
#endif
#ifdef GRID_COMMS_SHMEM
  me = shmem_my_pe();
#endif
  if ( me ) { 
    std::cout.setstate(std::ios::badbit);
  }
}

void Grid_unquiesce_nodes(void)
{
#ifdef GRID_COMMS_MPI
    std::cout.clear();
#endif
}


}

