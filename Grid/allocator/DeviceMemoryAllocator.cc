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

  DeviceMemoryAllocator() {
    initialized = false;
    base = 0;
    size = 0;
    offset = 0;
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
    
    size_t n_pages = (_size + DEVICE_MEMORY_ALLOCATOR_PAGE_SIZE - 1) / DEVICE_MEMORY_ALLOCATOR_PAGE_SIZE;
    size = n_pages * DEVICE_MEMORY_ALLOCATOR_PAGE_SIZE;
    std::cout << GridLogMessage << "Init with " <<
      size << " bytes" << std::endl;

    base = (char*)acceleratorAllocDeviceInternal(size);
    assert(base);
    
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
    
    std::cout << GridLogMessage << "Done" << std::endl;
    
    offset = 0;
    
    pages.resize(n_pages, 0);
    
    std::cout << GridLogMessage << "Pages initialized" << std::endl;
  
    initialized = true;
  }
};

static DeviceMemoryAllocator dma;

void *acceleratorAllocDevice(size_t bytes) {
  if (!dma.initialized)
    dma.Init(MemoryManager::DeviceMaxBytes);
  
  if (!bytes)
    bytes++;
  
  size_t n_pages = (bytes + DEVICE_MEMORY_ALLOCATOR_PAGE_SIZE - 1) / DEVICE_MEMORY_ALLOCATOR_PAGE_SIZE;
  
  // first check if
  auto & sm = dma.size_map[n_pages];
  if (sm.size()) {
    size_t index = sm.back();
    sm.pop_back();
    assert(dma.pages[index] == 0);
    dma.pages[index] = n_pages;

    std::cout << GridLogMessage << "Can re-use pointer for " << n_pages << " pages" << std::endl;
    
    return dma.base + index * DEVICE_MEMORY_ALLOCATOR_PAGE_SIZE;
  }

  size_t end = (dma.offset + n_pages) * DEVICE_MEMORY_ALLOCATOR_PAGE_SIZE;
  assert(end <= dma.size);
  dma.pages[dma.offset] = n_pages;

  void* ptr = dma.base + dma.offset * DEVICE_MEMORY_ALLOCATOR_PAGE_SIZE;
  dma.offset += n_pages;

  std::cout << GridLogMessage << (dma.size - end) / DEVICE_MEMORY_ALLOCATOR_PAGE_SIZE << " pages left to allocate ("
	    << (dma.size - end) * 100 / dma.size << "% free)" << std::endl;

  return ptr;

}

void acceleratorFreeDevice(void *ptr) {
  if (!dma.initialized)
    return;
  
  size_t index = ((size_t)((char*)ptr - dma.base)) / DEVICE_MEMORY_ALLOCATOR_PAGE_SIZE;
  size_t n_pages = dma.pages[index];
  //std::cout << GridLogMessage << "Freeing ptr " << ptr << " has " << n_pages << " pages" << std::endl;
  dma.pages[index] = 0;
  auto & sm = dma.size_map[n_pages];
  sm.push_back(index);
  
}
#endif

NAMESPACE_END(Grid);
