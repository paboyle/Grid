/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/AllocationCache.h

    Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>

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

// Move control to configure.ac and Config.h?

#define ALLOCATION_CACHE
#define GRID_ALLOC_ALIGN (2*1024*1024)
#define GRID_ALLOC_SMALL_LIMIT (4096)

/*Pinning pages is costly*/

class AllocationCache {
private:

  ////////////////////////////////////////////////////////////
  // For caching recently freed allocations
  ////////////////////////////////////////////////////////////
  typedef struct { 
    void *address;
    size_t bytes;
    int valid;
  } AllocationCacheEntry;

  static const int NallocCacheMax=128; 
  static const int NallocType=4;
  static AllocationCacheEntry Entries[NallocType][NallocCacheMax];
  static int Victim[NallocType];
  static int Ncache[NallocType];

  /////////////////////////////////////////////////
  // Free pool
  /////////////////////////////////////////////////
  static void *Insert(void *ptr,size_t bytes,int type) ;
  static void *Insert(void *ptr,size_t bytes,AllocationCacheEntry *entries,int ncache,int &victim) ;
  static void *Lookup(size_t bytes,int type) ;
  static void *Lookup(size_t bytes,AllocationCacheEntry *entries,int ncache) ;

  /////////////////////////////////////////////////
  // Internal device view
  /////////////////////////////////////////////////
  static void *AcceleratorAllocate(size_t bytes);
  static void  AcceleratorFree    (void *ptr,size_t bytes);
  static int   ViewVictim(void);
  static void  Evict(int e);
  static void  Flush(int e);
  static void  Clone(int e);
  static int   CpuViewLookup(void *CpuPtr);
  static int   AccViewLookup(void *AccPtr);

public:
  static void Init(void);

  static void  AccViewClose(void* AccPtr);
  static void  CpuViewClose(void* CpuPtr);
  static void *AccViewOpen(void* CpuPtr,size_t bytes,int mode,int transient);
  static void *CpuViewOpen(void* CpuPtr,size_t bytes,int mode,int transient);

  static void *CpuAllocate(size_t bytes);
  static void  CpuFree    (void *ptr,size_t bytes);
};

NAMESPACE_END(Grid);


