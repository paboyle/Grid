/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/DeviceMemoryAllocator.h

    Copyright (C) 2025

Author: Christoph Lehner <christoph@lhnr.de>

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

#include <Grid/GridCore.h>

NAMESPACE_BEGIN(Grid);

#define DEVICE_MEMORY_ALLOCATOR_PAGE_SIZE (64*1024)
#define OVERALLOCATION_FACTOR 1.2

#ifdef GRID_DEVICE_MEMORY_ALLOCATOR
struct DeviceMemoryAllocator {

  bool initialized;
  char* base;
  size_t size;
  size_t offset;
  bool verbose;

  DeviceMemoryAllocator() {
    initialized = false;
    base = 0;
    size = 0;
    offset = 0;
    verbose = false;
  }

  ~DeviceMemoryAllocator() {
    if (initialized) {
      acceleratorFreeDevice(base);
      initialized = false;
    }
  }

  std::vector<size_t> pages;
  std::map<size_t, std::vector<size_t> > size_map;

  void Init(size_t _size) {
    assert(!initialized);

    char* str;
    if ((str = getenv("GRID_OVERALLOCATION_FACTOR"))) {
      _size = (size_t)(_size * atof(str));
    } else {
      _size = (size_t)(_size * OVERALLOCATION_FACTOR);
    }

    verbose = (getenv("GRID_DEBUG_DEVICE_ALLOCATOR") != 0);
    
    size_t n_pages = (_size + DEVICE_MEMORY_ALLOCATOR_PAGE_SIZE - 1) / DEVICE_MEMORY_ALLOCATOR_PAGE_SIZE;
    size = n_pages * DEVICE_MEMORY_ALLOCATOR_PAGE_SIZE;
    std::cout << GridLogMessage << "Init device allocator with " << size << " bytes" << std::endl;

    base = (char*)acceleratorAllocDeviceInternal(size);
    assert(base);

    if (verbose)
      std::cout << GridLogMessage << "Initialize memory to zero" << std::endl;

    {
      uint64_t* ba = (uint64_t*)base;
      size_t n = size / sizeof(uint64_t);
      size_t MAX_BLOCK_INIT = 128*1024*1024;
      while (n > 0) {
	size_t n0 = n;
	if (n0 > MAX_BLOCK_INIT)
	  n0 = MAX_BLOCK_INIT;
	accelerator_for(i, n0, 1, {
	    ba[i] = (uint64_t)-1;
	  });
	ba += n0;
	n -= n0;
      }
    }  

    if (verbose)
      std::cout << GridLogMessage << "Done" << std::endl;
    
    offset = 0;
    
    pages.resize(n_pages, 0);

    if (verbose)
      std::cout << GridLogMessage << "Pages initialized" << std::endl;
  
    initialized = true;
  }

  void* attemptReuseExactSize(size_t n_pages) {
    auto sm = size_map.find(n_pages);
    if (sm != size_map.end() && sm->second.size() > 0) {
      size_t index = sm->second.back();
      sm->second.pop_back();
      
      if (sm->second.size() == 0)
	size_map.erase(sm);
      
      assert(pages[index] == 0);
      pages[index] = n_pages;

      return base + index * DEVICE_MEMORY_ALLOCATOR_PAGE_SIZE;
    }
    return 0;
  }

  void* attemptAllocUnused(size_t n_pages) {
    size_t end = (offset + n_pages) * DEVICE_MEMORY_ALLOCATOR_PAGE_SIZE;
    void* ptr = 0;
    
    if (end <= size) {
      pages[offset] = n_pages;
    
      ptr = base + offset * DEVICE_MEMORY_ALLOCATOR_PAGE_SIZE;
      offset += n_pages;

      if (verbose) {
	size_t reusable_pages = 0;
	for (auto & sm : size_map)
	  reusable_pages += sm.first * sm.second.size();
	
	std::cout << GridLogMessage << (size - end) / DEVICE_MEMORY_ALLOCATOR_PAGE_SIZE << " pages left to allocate ("
		  << (size - end) * 100 / size << "% unallocated, " << reusable_pages << " reusable pages)" << std::endl;
      }
    }

    return ptr;
  }

  void* alloc(size_t bytes) {
    if (!initialized)
      Init(MemoryManager::DeviceMaxBytes);
    
    if (!bytes)
      bytes++;
    
    size_t n_pages = (bytes + DEVICE_MEMORY_ALLOCATOR_PAGE_SIZE - 1) / DEVICE_MEMORY_ALLOCATOR_PAGE_SIZE;
    
    // first check if block of perfect size is available
    void* ptr;
    if ((ptr = attemptReuseExactSize(n_pages))) {
      
      if (verbose)
	std::cout << GridLogMessage << "Can re-use perfect pointer for " << n_pages << " pages" << std::endl;
      
      return ptr;
    }

    // if not, attempt to allocate in the unused area
    if ((ptr = attemptAllocUnused(n_pages)))
      return ptr;

    // last attempt, find a re-usable region that barely fits and return it
    // for loop of std::map iterates in ascending order
    size_t reusable_pages = 0;
    size_t n_pages_usable = 0;
    for (auto & sm : size_map) {
      assert(sm.second.size() > 0); // should never be empty
      reusable_pages += sm.first * sm.second.size();
      if (n_pages_usable == 0 && sm.first > n_pages)
	n_pages_usable = sm.first;
    }    

    if (n_pages_usable == 0) {
      std::cout << GridLogMessage << "Out of memory for " << n_pages << " pages!  Re-usable pages at time of death:" << std::endl;

      for (auto & sm : size_map) {
	std::cout << GridLogMessage << sm.second.size() << " x " << sm.first << " pages" << std::endl;
      }
	
      exit(1);
    }

    if ((ptr = attemptReuseExactSize(n_pages_usable))) {
      
      if (verbose)
	std::cout << GridLogMessage << "Can re-use pointer for " << n_pages_usable << " pages when " << n_pages << " were needed; " << reusable_pages << " reusable pages" << std::endl;
      
      return ptr;
    }

    // this should never be reached
    assert(0);
    return ptr;
  }

  void free(void* ptr) {
    if (!initialized)
      return;
    
    size_t index = ((size_t)((char*)ptr - base)) / DEVICE_MEMORY_ALLOCATOR_PAGE_SIZE;
    size_t n_pages = pages[index];
    //std::cout << GridLogMessage << "Freeing ptr " << ptr << " has " << n_pages << " pages" << std::endl;
    pages[index] = 0;
    auto & sm = size_map[n_pages];
    sm.push_back(index);
  }
};

static DeviceMemoryAllocator dma;

void *acceleratorAllocDevice(size_t bytes) {
  return dma.alloc(bytes);
}

void acceleratorFreeDevice(void *ptr) {
  dma.free(ptr);
}
#endif

NAMESPACE_END(Grid);
