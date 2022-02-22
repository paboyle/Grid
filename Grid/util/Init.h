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
#pragma once


NAMESPACE_BEGIN(Grid);

void Grid_init(int *argc,char ***argv);
void Grid_finalize(void);

// internal, controled with --handle
void Grid_sa_signal_handler(int sig,siginfo_t *si,void * ptr);
void Grid_debug_handler_init(void);
void Grid_quiesce_nodes(void);
void Grid_unquiesce_nodes(void);

const Coordinate  GridDefaultSimd(int dims,int nsimd);
const Coordinate &GridDefaultLatt(void);
const Coordinate &GridDefaultMpi(void);
const int        &GridThreads(void)  ;
void              GridSetThreads(int t) ;
void GridLogTimestamp(int);
void GridLogLayout();

// Common parsing chores
std::string GridCmdOptionPayload(char ** begin, char ** end, const std::string & option);
bool        GridCmdOptionExists(char** begin, char** end, const std::string& option);
template<class VectorInt>
std::string GridCmdVectorIntToString(const VectorInt & vec);
void GridCmdOptionCSL(std::string str,std::vector<std::string> & vec);
template<class VectorInt>
void GridCmdOptionIntVector(const std::string &str,VectorInt & vec);
void GridCmdOptionInt(std::string &str,int & val);
void GridCmdOptionFloat(std::string &str,float & val);


void GridParseLayout(char **argv,int argc,
		     std::vector<int> &latt,
		     std::vector<int> &simd,
		     std::vector<int> &mpi);

void printHash(void);


NAMESPACE_END(Grid);

