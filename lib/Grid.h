//
//  Grid.h
//  simd
//
//  Created by Peter Boyle on 09/05/2014.
//  Copyright (c) 2014 University of Edinburgh. All rights reserved.
//


#ifndef GRID_H
#define GRID_H

#include <cassert>

#include <complex>
#include <vector>

#include <iostream>
#include <iomanip>
#include <random>
#include <functional>

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <stdio.h>
#include <signal.h>

#ifndef MAX
#define MAX(x,y) ((x)>(y)?(x):(y))
#define MIN(x,y) ((x)>(y)?(y):(x))
#endif

#define strong_inline __attribute__((always_inline)) inline

#include <GridConfig.h>

#include <GridTime.h>
#include <GridLog.h>

////////////////////////////////////////////////////////////
// Tunable header includes
////////////////////////////////////////////////////////////

#ifdef HAVE_MALLOC_MALLOC_H
#include <malloc/malloc.h>
#endif
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif

#include <AlignedAllocator.h>

#include <Simd.h>
#include <Threads.h>

#include <Communicator.h> // subdir aggregate
#include <Cartesian.h> // subdir aggregate
#include <Tensors.h>   // subdir aggregate
#include <Lattice.h>   // subdir aggregate
#include <Cshift.h>    // subdir aggregate
#include <Stencil.h>   // subdir aggregate
#include <Algorithms.h>// subdir aggregate

#include <qcd/QCD.h>
#include <parallelIO/NerscIO.h>

namespace Grid {

  void Grid_init(int *argc,char ***argv);
  void Grid_finalize(void);
  // internal, controled with --handle
  void Grid_sa_signal_handler(int sig,siginfo_t *si,void * ptr);
  void Grid_debug_handler_init(void);
  void Grid_quiesce_nodes(void);
  void Grid_unquiesce_nodes(void);

  // C++11 time facilities better?
  double usecond(void);

  const std::vector<int> GridDefaultSimd(int dims,int nsimd);
  const std::vector<int> &GridDefaultLatt(void);
  const std::vector<int> &GridDefaultMpi(void);
  const int              &GridThreads(void)  ;
  void                 GridSetThreads(int t) ;

  // Common parsing chores
  std::string GridCmdOptionPayload(char ** begin, char ** end, const std::string & option);
  bool        GridCmdOptionExists(char** begin, char** end, const std::string& option);
  std::string GridCmdVectorIntToString(const std::vector<int> & vec);

  void GridParseLayout(char **argv,int argc,
		       std::vector<int> &latt,
		       std::vector<int> &simd,
		       std::vector<int> &mpi);


};

#endif
