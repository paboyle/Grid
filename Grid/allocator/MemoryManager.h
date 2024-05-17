/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/MemoryManager.h

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
#include <list> 
#include <unordered_map>  

NAMESPACE_BEGIN(Grid);

// Move control to configure.ac and Config.h?

#define GRID_ALLOC_SMALL_LIMIT (4096)
#define GRID_ALLOC_HUGE_LIMIT  (2147483648)

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)
#define FILE_LINE __FILE__ ":" TOSTRING(__LINE__)
#define AUDIT(a) MemoryManager::Audit(FILE_LINE)

/*Pinning pages is costly*/
////////////////////////////////////////////////////////////////////////////
// Advise the LatticeAccelerator class
////////////////////////////////////////////////////////////////////////////
enum ViewAdvise {
 AdviseDefault       = 0x0,    // Regular data
 AdviseInfrequentUse = 0x1     // Advise that the data is used infrequently.  This can
                               // significantly influence performance of bulk storage.
 
 // AdviseTransient      = 0x2,   // Data will mostly be read.  On some architectures
                               // enables read-only copies of memory to be kept on
                               // host and device.

 // AdviseAcceleratorWriteDiscard = 0x4  // Field will be written in entirety on device

};

////////////////////////////////////////////////////////////////////////////
// View Access Mode
////////////////////////////////////////////////////////////////////////////
enum ViewMode {
  AcceleratorRead  = 0x01,
  AcceleratorWrite = 0x02,
  AcceleratorWriteDiscard = 0x04,
  CpuRead  = 0x08,
  CpuWrite = 0x10,
  CpuWriteDiscard = 0x10 // same for now
};

struct MemoryStatus {
  uint64_t     DeviceBytes;
  uint64_t     DeviceLRUBytes;
  uint64_t     DeviceMaxBytes;
  uint64_t     HostToDeviceBytes;
  uint64_t     DeviceToHostBytes;
  uint64_t     HostToDeviceXfer;
  uint64_t     DeviceToHostXfer;
  uint64_t     DeviceEvictions;
  uint64_t     DeviceDestroy;
  uint64_t     DeviceAllocCacheBytes;
  uint64_t     HostAllocCacheBytes;
};


class MemoryManager {
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
  static const int NallocType=9;
  static AllocationCacheEntry Entries[NallocType][NallocCacheMax];
  static int Victim[NallocType];
  static int Ncache[NallocType];
  static uint64_t CacheBytes[NallocType];

  /////////////////////////////////////////////////
  // Free pool
  /////////////////////////////////////////////////
  static void *Insert(void *ptr,size_t bytes,int type) ;
  static void *Lookup(size_t bytes,int type) ;
  static void *Insert(void *ptr,size_t bytes,AllocationCacheEntry *entries,int ncache,int &victim,uint64_t &cbytes) ;
  static void *Lookup(size_t bytes,AllocationCacheEntry *entries,int ncache,uint64_t &cbytes) ;

 public:
  static void PrintBytes(void);
  static void Audit(std::string s);
  static void Init(void);
  static void InitMessage(void);
  static void *AcceleratorAllocate(size_t bytes);
  static void  AcceleratorFree    (void *ptr,size_t bytes);
  static void *SharedAllocate(size_t bytes);
  static void  SharedFree    (void *ptr,size_t bytes);
  static void *CpuAllocate(size_t bytes);
  static void  CpuFree    (void *ptr,size_t bytes);

  ////////////////////////////////////////////////////////
  // Footprint tracking
  ////////////////////////////////////////////////////////
  static uint64_t     DeviceBytes;
  static uint64_t     DeviceLRUBytes;
  static uint64_t     DeviceMaxBytes;
  static uint64_t     HostToDeviceBytes;
  static uint64_t     DeviceToHostBytes;
  static uint64_t     HostToDeviceXfer;
  static uint64_t     DeviceToHostXfer;
  static uint64_t     DeviceEvictions;
  static uint64_t     DeviceDestroy;
  
  static uint64_t     DeviceCacheBytes();
  static uint64_t     HostCacheBytes();

  static MemoryStatus GetFootprint(void) {
    MemoryStatus stat;
    stat.DeviceBytes       = DeviceBytes;
    stat.DeviceLRUBytes    = DeviceLRUBytes;
    stat.DeviceMaxBytes    = DeviceMaxBytes;
    stat.HostToDeviceBytes = HostToDeviceBytes;
    stat.DeviceToHostBytes = DeviceToHostBytes;
    stat.HostToDeviceXfer  = HostToDeviceXfer;
    stat.DeviceToHostXfer  = DeviceToHostXfer;
    stat.DeviceEvictions   = DeviceEvictions;
    stat.DeviceDestroy     = DeviceDestroy;
    stat.DeviceAllocCacheBytes = DeviceCacheBytes();
    stat.HostAllocCacheBytes   = HostCacheBytes();
    return stat;
  };
  
 private:
#ifndef GRID_UVM
  //////////////////////////////////////////////////////////////////////
  // Data tables for ViewCache
  //////////////////////////////////////////////////////////////////////
  typedef std::list<uint64_t> LRU_t;
  typedef typename LRU_t::iterator LRUiterator;
  typedef struct { 
    int        LRU_valid;
    LRUiterator LRU_entry;
    uint64_t CpuPtr;
    uint64_t AccPtr;
    size_t   bytes;
    uint32_t transient;
    uint32_t state;
    uint32_t accLock;
    uint32_t cpuLock;
  } AcceleratorViewEntry;
  
  typedef std::unordered_map<uint64_t,AcceleratorViewEntry> AccViewTable_t;
  typedef typename AccViewTable_t::iterator AccViewTableIterator ;

  static AccViewTable_t AccViewTable;
  static LRU_t LRU;

  /////////////////////////////////////////////////
  // Device motion
  /////////////////////////////////////////////////
  static void  Create(uint64_t CpuPtr,size_t bytes,ViewMode mode,ViewAdvise hint);
  static void  EvictVictims(uint64_t bytes); // Frees up <bytes>
  static void  Evict(AcceleratorViewEntry &AccCache);
  static void  Flush(AcceleratorViewEntry &AccCache);
  static void  Clone(AcceleratorViewEntry &AccCache);
  static void  AccDiscard(AcceleratorViewEntry &AccCache);
  static void  CpuDiscard(AcceleratorViewEntry &AccCache);

  //  static void  LRUupdate(AcceleratorViewEntry &AccCache);
  static void  LRUinsert(AcceleratorViewEntry &AccCache);
  static void  LRUremove(AcceleratorViewEntry &AccCache);
  
  // manage entries in the table
  static int                  EntryPresent(uint64_t CpuPtr);
  static void                 EntryCreate(uint64_t CpuPtr,size_t bytes,ViewMode mode,ViewAdvise hint);
  static void                 EntryErase (uint64_t CpuPtr);
  static AccViewTableIterator EntryLookup(uint64_t CpuPtr);
  static void                 EntrySet   (uint64_t CpuPtr,AcceleratorViewEntry &entry);

  static void     AcceleratorViewClose(uint64_t AccPtr);
  static uint64_t AcceleratorViewOpen(uint64_t  CpuPtr,size_t bytes,ViewMode mode,ViewAdvise hint);
  static void     CpuViewClose(uint64_t Ptr);
  static uint64_t CpuViewOpen(uint64_t  CpuPtr,size_t bytes,ViewMode mode,ViewAdvise hint);
#endif

 public:
  static void NotifyDeletion(void * CpuPtr);
  static void Print(void);
  static void PrintAll(void);
  static void PrintState( void* CpuPtr);
  static int   isOpen   (void* CpuPtr);
  static void  ViewClose(void* CpuPtr,ViewMode mode);
  static void *ViewOpen (void* CpuPtr,size_t bytes,ViewMode mode,ViewAdvise hint);

};

NAMESPACE_END(Grid);


