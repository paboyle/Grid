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
#include <iostream>
#include <cassert>
#include <random>
#include <functional>
#include <stdlib.h>
#include <sys/time.h>
#include <stdio.h>
#include <signal.h>

#include <Grid_config.h>

////////////////////////////////////////////////////////////
// Tunable header includes
////////////////////////////////////////////////////////////

#ifdef HAVE_OPENMP
#define OMP
#include <omp.h>
#endif

#ifdef HAVE_MALLOC_MALLOC_H
#include <malloc/malloc.h>
#endif

#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif

#include <Grid_aligned_allocator.h>
#include <Grid_simd.h>
#include <Grid_math.h>
#include <Grid_cartesian.h>
#include <Grid_lattice.h>
#include <Grid_comparison.h>
#include <Grid_cshift.h>
#include <Grid_where.h>
#include <Grid_stencil.h>
#include <qcd/Grid_QCD.h>
#include <parallelIO/GridNerscIO.h>

namespace Grid {

  void Grid_init(int *argc,char ***argv);
  void Grid_finalize(void);
  double usecond(void);
  void Grid_sa_signal_handler(int sig,siginfo_t *si,void * ptr);
  void Grid_debug_handler_init(void);

};

#endif
