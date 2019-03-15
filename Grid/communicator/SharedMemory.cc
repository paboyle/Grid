/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/communicator/SharedMemory.cc

    Copyright (C) 2015

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

#include <Grid/GridCore.h>

namespace Grid { 

// static data

uint64_t            GlobalSharedMemory::MAX_MPI_SHM_BYTES   = 1024LL*1024LL*1024LL; 
int                 GlobalSharedMemory::Hugepages = 0;
int                 GlobalSharedMemory::_ShmSetup;
int                 GlobalSharedMemory::_ShmAlloc;
uint64_t            GlobalSharedMemory::_ShmAllocBytes;

std::vector<void *> GlobalSharedMemory::WorldShmCommBufs;

Grid_MPI_Comm       GlobalSharedMemory::WorldShmComm;
int                 GlobalSharedMemory::WorldShmRank;
int                 GlobalSharedMemory::WorldShmSize;
std::vector<int>    GlobalSharedMemory::WorldShmRanks;

Grid_MPI_Comm       GlobalSharedMemory::WorldComm;
int                 GlobalSharedMemory::WorldSize;
int                 GlobalSharedMemory::WorldRank;

int                 GlobalSharedMemory::WorldNodes;
int                 GlobalSharedMemory::WorldNode;

void GlobalSharedMemory::SharedMemoryFree(void)
{
  assert(_ShmAlloc);
  assert(_ShmAllocBytes>0);
  for(int r=0;r<WorldShmSize;r++){
    munmap(WorldShmCommBufs[r],_ShmAllocBytes);
  }
  _ShmAlloc = 0;
  _ShmAllocBytes = 0;
}
/////////////////////////////////
// Alloc, free shmem region
/////////////////////////////////
void *SharedMemory::ShmBufferMalloc(size_t bytes){
  //  bytes = (bytes+sizeof(vRealD))&(~(sizeof(vRealD)-1));// align up bytes
  void *ptr = (void *)heap_top;
  heap_top  += bytes;
  heap_bytes+= bytes;
  if (heap_bytes >= heap_size) {
    std::cout<< " ShmBufferMalloc exceeded shared heap size -- try increasing with --shm <MB> flag" <<std::endl;
    std::cout<< " Parameter specified in units of MB (megabytes) " <<std::endl;
    std::cout<< " Current value is " << (heap_size/(1024*1024)) <<std::endl;
    assert(heap_bytes<heap_size);
  }
  return ptr;
}
void SharedMemory::ShmBufferFreeAll(void) { 
  heap_top  =(size_t)ShmBufferSelf();
  heap_bytes=0;
}
void *SharedMemory::ShmBufferSelf(void)
{
  return ShmCommBufs[ShmRank];
}



}
