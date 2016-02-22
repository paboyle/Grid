    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/Threads.h

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
#ifndef GRID_THREADS_H
#define GRID_THREADS_H

#ifdef _OPENMP
#define GRID_OMP
#endif

#define UNROLL  _Pragma("unroll")

#ifdef GRID_OMP
#include <omp.h>
#define PARALLEL_FOR_LOOP _Pragma("omp parallel for ")
#define PARALLEL_NESTED_LOOP2 _Pragma("omp parallel for collapse(2)")
#else
#define PARALLEL_FOR_LOOP 
#define PARALLEL_NESTED_LOOP2
#endif

namespace Grid {

  // Introduce a class to gain deterministic bit reproducible reduction.
  // make static; perhaps just a namespace is required.

class GridThread {
 public:
  static int _threads;
  static int _hyperthreads;
  static int _cores;

  static void SetCores(int cr) { 
#ifdef GRID_OMP
    _cores = cr;
#else 
    _cores = 1;
#endif
  }
  static void SetThreads(int thr) { 
#ifdef GRID_OMP
    _threads = MIN(thr,omp_get_max_threads()) ;
    omp_set_num_threads(_threads);
#else 
    _threads = 1;
#endif
  };
  static void SetMaxThreads(void) { 
#ifdef GRID_OMP
    //    setenv("KMP_AFFINITY","balanced",1);
    _threads = omp_get_max_threads();
    omp_set_num_threads(_threads);
#else 
    _threads = 1;
#endif
  };
  static int GetHyperThreads(void) { assert(_threads%_cores ==0); return _threads/_cores; };
  static int GetCores(void)   { return _cores; };
  static int GetThreads(void) { return _threads; };
  static int SumArraySize(void) {return _threads;};

  static void GetWork(int nwork, int me, int & mywork, int & myoff){
    GetWork(nwork,me,mywork,myoff,_threads);
  }
  static void GetWork(int nwork, int me, int & mywork, int & myoff,int units){
    int basework = nwork/units;
    int backfill = units-(nwork%units);
    if ( me >= units ) { 
      mywork = myoff = 0;
    } else { 
      mywork = (nwork+me)/units;
      myoff  = basework * me;
      if ( me > backfill ) 
	myoff+= (me-backfill);
    }
    return;
  };

  static void GetWorkBarrier(int nwork, int &me, int & mywork, int & myoff){
    me     = ThreadBarrier();
    GetWork(nwork,me,mywork,myoff);
  };

  static int  ThreadBarrier(void) {
#ifdef GRID_OMP
#pragma omp barrier
    return omp_get_thread_num();
#else
    return 0;
#endif
  };
  
  template<class obj> static void ThreadSum( std::vector<obj> &sum_array,obj &val,int me){
    sum_array[me] = val;
    val=zero;
    ThreadBarrier();
    for(int i=0;i<_threads;i++) val+= sum_array[i];
    ThreadBarrier();
  };

};

}
#endif
