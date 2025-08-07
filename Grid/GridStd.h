#ifndef GRID_STD_H
#define GRID_STD_H

///////////////////
// Grid config
///////////////////
#include "Config.h"

///////////////////
// Std C++ dependencies
///////////////////
#define _NBACKTRACE (256)
extern void * Grid_backtrace_buffer[_NBACKTRACE];

#include <cassert>
#include <complex>
#include <memory>
#include <vector>
#include <array>
#include <string>
#include <iostream>
#include <iomanip>
#include <random>
#include <functional>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <strings.h>
#include <stdio.h>
#include <signal.h>
#include <ctime>
#include <sys/time.h>
#include <chrono>
#include <zlib.h>
#ifdef HAVE_EXECINFO_H
#include <execinfo.h>
#endif

void GridAbort(void);

#define ASSLOG(A) ::write(STDERR_FILENO,A,strlen(A));
#ifdef HAVE_EXECINFO_H
#define GRID_ASSERT(b) if(!(b)) {					\
    ASSLOG(" GRID_ASSERT failure: ");					\
    ASSLOG(__FILE__);							\
    ASSLOG(" : ");							\
    ASSLOG(#b);								\
    ASSLOG(" : ");							\
    int symbols = backtrace(Grid_backtrace_buffer,_NBACKTRACE);		\
    backtrace_symbols_fd(Grid_backtrace_buffer,symbols,STDERR_FILENO);	\
    GridAbort();							\
  };
#else
#define GRID_ASSERT(b) if(!(b)) {					\
    ASSLOG(" GRID_ASSERT failure: ");					\
    ASSLOG(__FILE__);							\
    ASSLOG(" : ");							\
    ASSLOG(#b);								\
    ASSLOG(" : ");							\
    GridAbort();							\
  };
#endif


#ifdef TOFU
#undef GRID_COMMS_THREADS
#endif
#endif /* GRID_STD_H */
