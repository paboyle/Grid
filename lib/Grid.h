//
//  Grid.h
//  simd
//
//  Created by Peter Boyle on 09/05/2014.
//  Copyright (c) 2014 University of Edinburgh. All rights reserved.
//


#ifndef GRID_H
#define GRID_H

#include <stdio.h>

#include <complex>
#include <vector>
#include <valarray>
#include <iostream>
#include <cassert>
#include <random>
#include <functional>
#include <stdlib.h>
#include <sys/time.h>
#include <stdio.h>
#include <signal.h>

#ifndef MAX
#define MAX(x,y) ((x)>(y)?(x):(y))
#define MIN(x,y) ((x)>(y)?(y):(x))
#endif

#define strong_inline __attribute__((always_inline)) inline

#include <Grid_config.h>

////////////////////////////////////////////////////////////
// Tunable header includes
////////////////////////////////////////////////////////////

#ifdef HAVE_MALLOC_MALLOC_H
#include <malloc/malloc.h>
#endif
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif

#include <Grid_aligned_allocator.h>
#include <Grid_simd.h>
#include <Grid_threads.h>

#include <Grid_cartesian.h>

#include <Grid_math.h>
#include <Grid_lattice.h>
#include <Grid_comparison.h>
#include <Grid_cshift.h>
#include <Grid_stencil.h>
#include <qcd/Grid_qcd.h>
#include <parallelIO/GridNerscIO.h>

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
  void        GridParseIntVector(std::string &str,std::vector<int> & vec);

  void GridParseLayout(char **argv,int argc,
		       std::vector<int> &latt,
		       std::vector<int> &simd,
		       std::vector<int> &mpi);


};

#endif
