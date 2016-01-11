    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/PerfCount.h

    Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
Author: Peter Boyle <peterboyle@MacBook-Pro.local>
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
#ifndef GRID_PERFCOUNT_H
#define GRID_PERFCOUNT_H

#include <sys/time.h>
#include <ctime>
#include <chrono>
#include <string.h>

#include <sys/ioctl.h>

#ifdef __linux__
#include <syscall.h>
#include <linux/perf_event.h>
#else
#include <sys/syscall.h>
#endif
namespace Grid {


#ifdef __linux__
static long perf_event_open(struct perf_event_attr *hw_event, pid_t pid,
			    int cpu, int group_fd, unsigned long flags)
{
  int ret=0;

  ret = syscall(__NR_perf_event_open, hw_event, pid, cpu,
		group_fd, flags);
  return ret;
}
#endif


class PerformanceCounter {
private:

  typedef struct { 
  public:
    uint32_t type;
    uint64_t config;
    const char *name;
  } PerformanceCounterConfig; 
  
  static const PerformanceCounterConfig PerformanceCounterConfigs [];

public:

  enum PerformanceCounterType {
    CPUCYCLES=0,
    INSTRUCTIONS,
    //    STALL_CYCLES,
    CACHE_REFERENCES,
    CACHE_MISSES,
    L1D_READ_MISS,
    L1D_READ_ACCESS,
    L1D_WRITE_MISS,
    L1D_WRITE_ACCESS,
    L1D_PREFETCH_MISS,
    L1D_PREFETCH_ACCESS,
    LL_READ_MISS,
    //    LL_READ_ACCESS,
    LL_WRITE_MISS,
    LL_WRITE_ACCESS,
    LL_PREFETCH_MISS,
    LL_PREFETCH_ACCESS,
    L1I_READ_MISS,
    L1I_READ_ACCESS,
    PERFORMANCE_COUNTER_NUM_TYPES
  };

public:
    
  int PCT;

  long long count;
  int fd;
  uint64_t elapsed;
  uint64_t begin;

  static int NumTypes(void){ 
    return PERFORMANCE_COUNTER_NUM_TYPES;
  }

  PerformanceCounter(int _pct) {
#ifdef __linux__
    assert(_pct>=0);
    assert(_pct<PERFORMANCE_COUNTER_NUM_TYPES);
    fd=-1;
    count=0;
    PCT =_pct;
    Open();
#endif
  }
  void Open(void) 
  {
#ifdef __linux__
    struct perf_event_attr pe;
    memset(&pe, 0, sizeof(struct perf_event_attr));
    pe.size = sizeof(struct perf_event_attr);

    pe.disabled = 1;
    pe.exclude_kernel = 1;
    pe.exclude_hv = 1;
    pe.inherit    = 1;

    pe.type  = PerformanceCounterConfigs[PCT].type;
    pe.config= PerformanceCounterConfigs[PCT].config;
    const char * name = PerformanceCounterConfigs[PCT].name;
    fd = perf_event_open(&pe, 0, -1, -1, 0); // pid 0, cpu -1 current process any cpu. group -1
    if (fd == -1) {
      fprintf(stderr, "Error opening leader %llx for event %s\n", pe.config,name);
      perror("Error is");
    }
#endif
  }

  void Start(void)
  {
#ifdef __linux__
    if ( fd!= -1) {
      ioctl(fd, PERF_EVENT_IOC_RESET, 0);
      ioctl(fd, PERF_EVENT_IOC_ENABLE, 0);
    }
    begin  =__rdtsc();
#else
    begin = 0;
#endif
  }

  void Stop(void) {
    count=0;
#ifdef __linux__
    if ( fd!= -1) {
      ioctl(fd, PERF_EVENT_IOC_DISABLE, 0);
      ::read(fd, &count, sizeof(long long));
    }
    elapsed = __rdtsc() - begin;
#else
    elapsed = 0;
#endif

  }
  void Report(void) {
#ifdef __linux__
    printf("%llu cycles %s = %20llu\n", elapsed , PerformanceCounterConfigs[PCT].name, count);
#else
    printf("%llu cycles \n", elapsed );
#endif
  }

  ~PerformanceCounter()
  {
#ifdef __linux__
    close(fd);
#endif
  }

};

}
#endif
