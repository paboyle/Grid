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

/*Construct from an MPI communicator*/
void GlobalSharedMemory::Init(Grid_MPI_Comm comm)
{
  WorldComm = 0;
  WorldRank = 0;
  WorldSize = 1;
  WorldShmComm = 0 ;
  WorldShmRank = 0 ;
  WorldShmSize = 1 ;
  WorldNodes   = 1 ;
  WorldNode    = 0 ;
  WorldShmRanks.resize(WorldSize); WorldShmRanks[0] = 0;
  WorldShmCommBufs.resize(1);
}

void GlobalSharedMemory::OptimalCommunicator(const std::vector<int> &processors,Grid_MPI_Comm & optimal_comm)
{
  optimal_comm = WorldComm;
}

////////////////////////////////////////////////////////////////////////////////////////////
// Hugetlbfs mapping intended, use anonymous mmap
////////////////////////////////////////////////////////////////////////////////////////////
void GlobalSharedMemory::SharedMemoryAllocate(uint64_t bytes, int flags)
{
  void * ShmCommBuf ; 
  MAX_MPI_SHM_BYTES=bytes;
  int mmap_flag =0;
#ifdef MAP_ANONYMOUS
  mmap_flag = mmap_flag| MAP_SHARED | MAP_ANONYMOUS;
#endif
#ifdef MAP_ANON
  mmap_flag = mmap_flag| MAP_SHARED | MAP_ANON;
#endif
#ifdef MAP_HUGETLB
  if ( flags ) mmap_flag |= MAP_HUGETLB;
#endif
  ShmCommBuf =(void *) mmap(NULL, bytes, PROT_READ | PROT_WRITE, mmap_flag, -1, 0); 
  if (ShmCommBuf == (void *)MAP_FAILED) {
    perror("mmap failed ");
    exit(EXIT_FAILURE);  
  }
#ifdef MADV_HUGEPAGE
  if (!Hugepages ) madvise(ShmCommBuf,bytes,MADV_HUGEPAGE);
#endif
  bzero(ShmCommBuf,bytes);
  WorldShmCommBufs[0] = ShmCommBuf;
};

void GlobalSharedMemory::SharedMemoryFree(void)
{
  assert(ShmSetup);
  assert(0); // unimplemented
}

  ////////////////////////////////////////////////////////
  // Global shared functionality finished
  // Now move to per communicator functionality
  ////////////////////////////////////////////////////////
void SharedMemory::SetCommunicator(Grid_MPI_Comm comm)
{
  ShmRanks.resize(1);
  ShmCommBufs.resize(1);
  ShmRanks[0] = 0;
  ShmRank     = 0;
  ShmSize     = 1;
  //////////////////////////////////////////////////////////////////////
  // Map ShmRank to WorldShmRank and use the right buffer
  //////////////////////////////////////////////////////////////////////
  ShmCommBufs[0] = GlobalSharedMemory::WorldShmCommBufs[0];
  heap_size      = GlobalSharedMemory::MAX_MPI_SHM_BYTES;
  ShmBufferFreeAll();
  return;
}
//////////////////////////////////////////////////////////////////
// On node barrier
//////////////////////////////////////////////////////////////////
void SharedMemory::ShmBarrier(void){ return ; }

//////////////////////////////////////////////////////////////////////////////////////////////////////////
// Test the shared memory is working
//////////////////////////////////////////////////////////////////////////////////////////////////////////
void SharedMemory::SharedMemoryTest(void) { return; }

void *SharedMemory::ShmBufferSelf(void)
{
  return ShmCommBufs[ShmRank];
}
void *SharedMemory::ShmBuffer(int rank)
{
  return NULL;
}
void *SharedMemory::ShmBufferTranslate(int rank,void * local_p)
{
  return NULL;
}

/////////////////////////////////
// Alloc, free shmem region ; common to MPI and none?
/////////////////////////////////
void *SharedMemory::ShmBufferMalloc(size_t bytes){
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

}
