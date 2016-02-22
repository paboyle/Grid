    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/PerfCount.cc

    Copyright (C) 2015

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
#include <PerfCount.h>

namespace Grid {

#define CacheControl(L,O,R) ((PERF_COUNT_HW_CACHE_##L)|(PERF_COUNT_HW_CACHE_OP_##O<<8)| (PERF_COUNT_HW_CACHE_RESULT_##R<<16))

const PerformanceCounter::PerformanceCounterConfig PerformanceCounter::PerformanceCounterConfigs [] = {
#ifdef __linux__
  { PERF_TYPE_HARDWARE, PERF_COUNT_HW_CPU_CYCLES          ,  "CPUCYCLES.........." },
  { PERF_TYPE_HARDWARE, PERF_COUNT_HW_INSTRUCTIONS        ,  "INSTRUCTIONS......." },
  { PERF_TYPE_HARDWARE, PERF_COUNT_HW_CACHE_REFERENCES    ,  "CACHE_REFERENCES..." },
  { PERF_TYPE_HARDWARE, PERF_COUNT_HW_CACHE_MISSES        ,  "CACHE_MISSES......." },
  { PERF_TYPE_HW_CACHE, CacheControl(L1D,READ,MISS)       ,  "L1D_READ_MISS......"},
  { PERF_TYPE_HW_CACHE, CacheControl(L1D,READ,ACCESS)     ,  "L1D_READ_ACCESS...."},
  { PERF_TYPE_HW_CACHE, CacheControl(L1D,WRITE,MISS)      ,  "L1D_WRITE_MISS....."},
  { PERF_TYPE_HW_CACHE, CacheControl(L1D,WRITE,ACCESS)    ,  "L1D_WRITE_ACCESS..."},
  { PERF_TYPE_HW_CACHE, CacheControl(L1D,PREFETCH,MISS)   ,  "L1D_PREFETCH_MISS.."},
  { PERF_TYPE_HW_CACHE, CacheControl(L1D,PREFETCH,ACCESS) ,  "L1D_PREFETCH_ACCESS"},
  { PERF_TYPE_HW_CACHE, CacheControl(LL,READ,MISS)        ,  "LL_READ_MISS......."},
  //  { PERF_TYPE_HW_CACHE, CacheControl(LL,READ,ACCESS)      ,  "LL_READ_ACCESS....."},
  { PERF_TYPE_HW_CACHE, CacheControl(LL,WRITE,MISS)       ,  "LL_WRITE_MISS......"},
  { PERF_TYPE_HW_CACHE, CacheControl(LL,WRITE,ACCESS)     ,  "LL_WRITE_ACCESS...."},
  { PERF_TYPE_HW_CACHE, CacheControl(LL,PREFETCH,MISS)    ,  "LL_PREFETCH_MISS..."},
  { PERF_TYPE_HW_CACHE, CacheControl(LL,PREFETCH,ACCESS)  ,  "LL_PREFETCH_ACCESS."},
  { PERF_TYPE_HW_CACHE, CacheControl(L1I,READ,MISS)       ,  "L1I_READ_MISS......"},
  { PERF_TYPE_HW_CACHE, CacheControl(L1I,READ,ACCESS)     ,  "L1I_READ_ACCESS...."}
#endif
  //  { PERF_TYPE_HARDWARE, PERF_COUNT_HW_STALLED_CYCLES_FRONTEND, "STALL_CYCLES" },
};
}
