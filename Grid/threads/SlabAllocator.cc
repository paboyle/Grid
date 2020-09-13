/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./Grid/threads/SlabAllocator.cc

    Copyright (C) 2020

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

#include <unordered_set>

NAMESPACE_BEGIN(Grid);

#ifdef GRID_CUDA

#define GRID_DEVICE_HEAP_SLAB_THRESHOLD (1024*1024)
#define GRID_DEVICE_HEAP_SLAB_SIZE (2*1024*1024)

void *acceleratorAllocDeviceCUDA(size_t bytes) {
  void *ptr=NULL;
  auto err = cudaMalloc((void **)&ptr,bytes);
  if( err != cudaSuccess ) {
    ptr = (void *) NULL;
    printf(" cudaMalloc failed for %d %s \n",bytes,cudaGetErrorString(err));
  }
  return ptr;
}

void acceleratorFreeDeviceCUDA(void *ptr) {
  cudaFree(ptr);
}

struct grid_device_heap_slab_t {
  void* Ptr;
  size_t ElementSize;
  size_t Elements;
  std::unordered_set<uint32_t> Allocated;
  std::unordered_set<uint32_t> Available;
};

std::unordered_map<void*, grid_device_heap_slab_t*> DeviceHeapPtrTable;
std::unordered_map<size_t, std::unordered_set<grid_device_heap_slab_t*> > DeviceHeapSlabTable;

void* SlabAllocateElement(grid_device_heap_slab_t* slab) {
  assert(!slab->Available.empty());
  auto available = slab->Available.begin();
  auto slot = *available;
  slab->Allocated.insert(slot);
  slab->Available.erase(available);
  
  void* Ptr = (void*)((char*)slab->Ptr + slot * slab->ElementSize);
  DeviceHeapPtrTable[Ptr] = slab;
  
  //std::cout << "Allocate element " << slot << " of slab " << slab << " of size " << slab->ElementSize << " with elements " << slab->Elements << 
  //  " (allocated = " << slab->Allocated.size() << ", available = " << slab->Available.size() << ")" << std::endl;
  
  return Ptr;
}

void SlabRemove(grid_device_heap_slab_t* slab) {
  auto & t = DeviceHeapSlabTable[slab->ElementSize];
  assert(slab->Ptr);
  DeviceHeapPtrTable.erase(slab->Ptr);
  acceleratorFreeDeviceCUDA(slab->Ptr);
  assert(t.count(slab) == 1);
  t.erase(slab);
  delete slab;
  //std::cout << "Remove slab " << slab << std::endl;
}

void SlabFreeElement(grid_device_heap_slab_t* slab, void* ElementPtr) {
  size_t Offset = (size_t)ElementPtr - (size_t)slab->Ptr;
  //std::cout << "SlabFreeElement offset " << Offset << std::endl;
  assert(Offset < GRID_DEVICE_HEAP_SLAB_SIZE);
  assert(Offset % slab->ElementSize == 0);
  size_t slot = Offset / slab->ElementSize;
  assert(slot >= 0);
  assert(slab->Allocated.count(slot) == 1 && slab->Available.count(slot) == 0);
  slab->Allocated.erase(slot);
  slab->Available.insert(slot);

  //std::cout << "Free element " << slot << " of slab" << slab << std::endl;
  
  if (slab->Allocated.empty()) {
    SlabRemove(slab);
  }
}

grid_device_heap_slab_t* SlabFind(size_t bytes) {

  grid_device_heap_slab_t* slab = 0;
  std::unordered_set<grid_device_heap_slab_t*>* slab_set = 0;

  decltype(DeviceHeapSlabTable.begin()) slabs = DeviceHeapSlabTable.find(bytes);
  if (slabs == DeviceHeapSlabTable.end()) {
    slab_set = &DeviceHeapSlabTable[bytes];
  } else {
    slab_set = &slabs->second;
  }

  for (auto& s : *slab_set) {
    if (!s->Available.empty()) {
      slab = &(*s);
      break;
    }
  }  

  if (!slab) {
    slab = new grid_device_heap_slab_t;
    slab_set->insert(slab);
    slab->Ptr = acceleratorAllocDeviceCUDA(GRID_DEVICE_HEAP_SLAB_SIZE);
    slab->ElementSize = bytes;
    slab->Elements = GRID_DEVICE_HEAP_SLAB_SIZE / bytes;
    for (size_t i=0;i<slab->Elements;i++)
      slab->Available.insert(i);
    //std::cout << "New slab" << slab << std::endl;
  }

  return slab;
}

void *acceleratorAllocDevice(size_t bytes) {
  if (bytes >= GRID_DEVICE_HEAP_SLAB_THRESHOLD) {
    return acceleratorAllocDeviceCUDA(bytes);
  }

  return SlabAllocateElement(SlabFind(bytes));
}

void acceleratorFreeDevice(void *ptr) {
  auto p = DeviceHeapPtrTable.find(ptr);
  if (p == DeviceHeapPtrTable.end()) {
    acceleratorFreeDeviceCUDA(ptr);
  } else {
    SlabFreeElement(p->second,ptr);
  }
}

#endif

NAMESPACE_END(Grid);
