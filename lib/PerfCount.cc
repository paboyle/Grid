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

#include <Grid/Grid.h>
#include <Grid/PerfCount.h>

namespace Grid {

#define CacheControl(L,O,R) ((PERF_COUNT_HW_CACHE_##L)|(PERF_COUNT_HW_CACHE_OP_##O<<8)| (PERF_COUNT_HW_CACHE_RESULT_##R<<16))
#define RawConfig(A,B) (A<<8|B)
const PerformanceCounter::PerformanceCounterConfig PerformanceCounter::PerformanceCounterConfigs [] = {
#ifdef __linux__
  { PERF_TYPE_HARDWARE, PERF_COUNT_HW_CACHE_REFERENCES    ,  "CACHE_REFERENCES..." , INSTRUCTIONS},
  { PERF_TYPE_HARDWARE, PERF_COUNT_HW_CACHE_MISSES        ,  "CACHE_MISSES......." , CACHE_REFERENCES},
  { PERF_TYPE_HARDWARE, PERF_COUNT_HW_CPU_CYCLES          ,  "CPUCYCLES.........." , INSTRUCTIONS},
  { PERF_TYPE_HARDWARE, PERF_COUNT_HW_INSTRUCTIONS        ,  "INSTRUCTIONS......." , CPUCYCLES   },
    // 4
#ifdef AVX512
    { PERF_TYPE_RAW, RawConfig(0x40,0x04), "ALL_LOADS..........", CPUCYCLES    },
    { PERF_TYPE_RAW, RawConfig(0x01,0x04), "L1_MISS_LOADS......", L1D_READ_ACCESS  },
    { PERF_TYPE_RAW, RawConfig(0x40,0x04), "ALL_LOADS..........", L1D_READ_ACCESS    },
    { PERF_TYPE_RAW, RawConfig(0x02,0x04), "L2_HIT_LOADS.......", L1D_READ_ACCESS  },
    { PERF_TYPE_RAW, RawConfig(0x04,0x04), "L2_MISS_LOADS......", L1D_READ_ACCESS  },
    { PERF_TYPE_RAW, RawConfig(0x10,0x04), "UTLB_MISS_LOADS....", L1D_READ_ACCESS },
    { PERF_TYPE_RAW, RawConfig(0x08,0x04), "DTLB_MISS_LOADS....", L1D_READ_ACCESS },
    // 11
#else
  { PERF_TYPE_HW_CACHE, CacheControl(L1D,READ,ACCESS)     ,  "L1D_READ_ACCESS....",INSTRUCTIONS},
  { PERF_TYPE_HW_CACHE, CacheControl(L1D,READ,MISS)       ,  "L1D_READ_MISS......",L1D_READ_ACCESS},
  { PERF_TYPE_HW_CACHE, CacheControl(L1D,WRITE,MISS)      ,  "L1D_WRITE_MISS.....",L1D_READ_ACCESS},
  { PERF_TYPE_HW_CACHE, CacheControl(L1D,WRITE,ACCESS)    ,  "L1D_WRITE_ACCESS...",L1D_READ_ACCESS},
  { PERF_TYPE_HW_CACHE, CacheControl(L1D,PREFETCH,MISS)   ,  "L1D_PREFETCH_MISS..",L1D_READ_ACCESS},
  { PERF_TYPE_HW_CACHE, CacheControl(L1D,PREFETCH,ACCESS) ,  "L1D_PREFETCH_ACCESS",L1D_READ_ACCESS},
  { PERF_TYPE_HW_CACHE, CacheControl(L1D,PREFETCH,ACCESS) ,  "L1D_PREFETCH_ACCESS",L1D_READ_ACCESS},
    // 11
#endif
  { PERF_TYPE_HW_CACHE, CacheControl(LL,READ,MISS)        ,  "LL_READ_MISS.......",L1D_READ_ACCESS},
  { PERF_TYPE_HW_CACHE, CacheControl(LL,READ,ACCESS)      ,  "LL_READ_ACCESS.....",L1D_READ_ACCESS},
  { PERF_TYPE_HW_CACHE, CacheControl(LL,WRITE,MISS)       ,  "LL_WRITE_MISS......",L1D_READ_ACCESS},
  { PERF_TYPE_HW_CACHE, CacheControl(LL,WRITE,ACCESS)     ,  "LL_WRITE_ACCESS....",L1D_READ_ACCESS},
    //15
  { PERF_TYPE_HW_CACHE, CacheControl(LL,PREFETCH,MISS)    ,  "LL_PREFETCH_MISS...",L1D_READ_ACCESS},
  { PERF_TYPE_HW_CACHE, CacheControl(LL,PREFETCH,ACCESS)  ,  "LL_PREFETCH_ACCESS.",L1D_READ_ACCESS},
  { PERF_TYPE_HW_CACHE, CacheControl(L1I,READ,MISS)       ,  "L1I_READ_MISS......",INSTRUCTIONS},
  { PERF_TYPE_HW_CACHE, CacheControl(L1I,READ,ACCESS)     ,  "L1I_READ_ACCESS....",INSTRUCTIONS}
    //19
  //  { PERF_TYPE_HARDWARE, PERF_COUNT_HW_STALLED_CYCLES_FRONTEND, "STALL_CYCLES" },
#endif
};
}
