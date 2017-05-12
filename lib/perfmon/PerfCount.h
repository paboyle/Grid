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
#include <unistd.h>
#include <sys/ioctl.h>

#ifdef __linux__
#include <syscall.h>
#include <linux/perf_event.h>
#else
#include <sys/syscall.h>
#endif
#ifdef __x86_64__
#include <x86intrin.h>
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

#ifdef TIMERS_OFF


inline uint64_t cyclecount(void){ 
  return 0;
}
#define __SSC_MARK(mark) __asm__ __volatile__ ("movl %0, %%ebx; .byte 0x64, 0x67, 0x90 " ::"i"(mark):"%ebx")
#define __SSC_STOP  __SSC_MARK(0x110)
#define __SSC_START __SSC_MARK(0x111)


#else

#define __SSC_MARK(mark) 
#define __SSC_STOP  
#define __SSC_START 

/*
 * cycle counters arch dependent
 */

#ifdef __bgq__
inline uint64_t cyclecount(void){ 
   uint64_t tmp;
   asm volatile ("mfspr %0,0x10C" : "=&r" (tmp)  );
   return tmp;
}
#elif defined __x86_64__
inline uint64_t cyclecount(void){ 
  return __rdtsc();
  //  unsigned int dummy;
  // return __rdtscp(&dummy);
}
#else

inline uint64_t cyclecount(void){ 
   return 0;
}

#endif

#endif

class PerformanceCounter {
private:

  typedef struct { 
  public:
    uint32_t type;
    uint64_t config;
    const char *name;
    int normalisation;
  } PerformanceCounterConfig; 
  
  static const PerformanceCounterConfig PerformanceCounterConfigs [];

public:

  enum PerformanceCounterType {
    CACHE_REFERENCES=0,
    CACHE_MISSES=1,
    CPUCYCLES=2,
    INSTRUCTIONS=3,
    L1D_READ_ACCESS=4,
    PERFORMANCE_COUNTER_NUM_TYPES=19
  };

public:
    
  int PCT;

  long long count;
  long long cycles;
  int fd;
  int cyclefd;
  unsigned long long elapsed;
  uint64_t begin;

  static int NumTypes(void){ 
    return PERFORMANCE_COUNTER_NUM_TYPES;
  }

  PerformanceCounter(int _pct) {
#ifdef __linux__
    assert(_pct>=0);
    assert(_pct<PERFORMANCE_COUNTER_NUM_TYPES);
    fd=-1;
    cyclefd=-1;
    count=0;
    cycles=0;
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
      fprintf(stderr, "Error opening leader %llx for event %s\n",(long long) pe.config,name);
      perror("Error is");
    }
    int norm = PerformanceCounterConfigs[PCT].normalisation;
    pe.type  = PerformanceCounterConfigs[norm].type;
    pe.config= PerformanceCounterConfigs[norm].config;
    name = PerformanceCounterConfigs[norm].name;
    cyclefd = perf_event_open(&pe, 0, -1, -1, 0); // pid 0, cpu -1 current process any cpu. group -1
    if (cyclefd == -1) {
      fprintf(stderr, "Error opening leader %llx for event %s\n",(long long) pe.config,name);
      perror("Error is");
    }
#endif
  }

  void Start(void)
  {
#ifdef __linux__
    if ( fd!= -1) {
      ::ioctl(fd, PERF_EVENT_IOC_RESET, 0);
      ::ioctl(fd, PERF_EVENT_IOC_ENABLE, 0);
      ::ioctl(cyclefd, PERF_EVENT_IOC_RESET, 0);
      ::ioctl(cyclefd, PERF_EVENT_IOC_ENABLE, 0);
    }
    begin  =cyclecount();
#else
    begin = 0;
#endif
  }

  void Stop(void) {
    count=0;
    cycles=0;
#ifdef __linux__
    ssize_t ign;
    if ( fd!= -1) {
      ::ioctl(fd, PERF_EVENT_IOC_DISABLE, 0);
      ::ioctl(cyclefd, PERF_EVENT_IOC_DISABLE, 0);
      ign=::read(fd, &count, sizeof(long long));
      ign+=::read(cyclefd, &cycles, sizeof(long long));
      assert(ign=2*sizeof(long long));
    }
    elapsed = cyclecount() - begin;
#else
    elapsed = 0;
#endif

  }
  void Report(void) {
#ifdef __linux__
    int N = PerformanceCounterConfigs[PCT].normalisation;
    const char * sn = PerformanceCounterConfigs[N].name ;
    const char * sc = PerformanceCounterConfigs[PCT].name;
      std::printf("tsc = %llu %s = %llu  %s = %20llu\n (%s/%s) rate = %lf\n", elapsed,sn ,cycles, 
		  sc, count, sc,sn, (double)count/(double)cycles);
#else
    std::printf("%llu cycles \n", elapsed );
#endif
  }

  ~PerformanceCounter()
  {
#ifdef __linux__
    ::close(fd);    ::close(cyclefd);
#endif
  }

};

}
#endif
