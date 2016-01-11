    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/Init.h

    Copyright (C) 2015

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
#ifndef GRID_INIT_H
#define GRID_INIT_H

namespace Grid {

  void Grid_init(int *argc,char ***argv);
  void Grid_finalize(void);
  // internal, controled with --handle
  void Grid_sa_signal_handler(int sig,siginfo_t *si,void * ptr);
  void Grid_debug_handler_init(void);
  void Grid_quiesce_nodes(void);
  void Grid_unquiesce_nodes(void);

  const std::vector<int> GridDefaultSimd(int dims,int nsimd);
  const std::vector<int> &GridDefaultLatt(void);
  const std::vector<int> &GridDefaultMpi(void);
  const int              &GridThreads(void)  ;
  void                    GridSetThreads(int t) ;

  // Common parsing chores
  std::string GridCmdOptionPayload(char ** begin, char ** end, const std::string & option);
  bool        GridCmdOptionExists(char** begin, char** end, const std::string& option);
  std::string GridCmdVectorIntToString(const std::vector<int> & vec);
  void GridCmdOptionCSL(std::string str,std::vector<std::string> & vec);
  void GridCmdOptionIntVector(std::string &str,std::vector<int> & vec);

  void GridParseLayout(char **argv,int argc,
		       std::vector<int> &latt,
		       std::vector<int> &simd,
		       std::vector<int> &mpi);


};
#endif
