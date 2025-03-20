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
#pragma once 

#include <Grid/GridCore.h>

#if defined (GRID_COMMS_MPI3) 
#include <mpi.h>
#endif 
#include <semaphore.h>
#include <fcntl.h>
#include <unistd.h>
#include <limits.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <sys/mman.h>
#include <zlib.h>

NAMESPACE_BEGIN(Grid);

#if defined (GRID_COMMS_MPI3) 
typedef MPI_Comm    Grid_MPI_Comm;
typedef MPI_Request MpiCommsRequest_t;
#ifdef ACCELERATOR_AWARE_MPI
typedef MPI_Request CommsRequest_t;
#else
/*
 * Enable state transitions as each packet flows.
 */
enum PacketType_t {
  FaceGather,
  InterNodeXmit,
  InterNodeRecv,
  IntraNodeXmit,
  IntraNodeRecv,
  InterNodeXmitISend,
  InterNodeReceiveHtoD
};
/*
 *Package arguments needed for various actions along packet flow
 */
typedef struct {
  PacketType_t PacketType;
  void *host_buf;
  void *device_buf;
  int dest;
  int tag;
  int commdir;
  unsigned long bytes;
  acceleratorEvent_t ev;
  MpiCommsRequest_t req;
} CommsRequest_t;
#endif

#else 
typedef int MpiCommsRequest_t;
typedef int CommsRequest_t;
typedef int Grid_MPI_Comm;
#endif

class GlobalSharedMemory {
private:
  static const int     MAXLOG2RANKSPERNODE = 16;            


  // Init once lock on the buffer allocation
  static int      _ShmSetup;
  static int      _ShmAlloc;
  static uint64_t _ShmAllocBytes;

public:
  ///////////////////////////////////////
  // HPE 8600 hypercube optimisation
  ///////////////////////////////////////
  static int HPEhypercube;

  static int      ShmSetup(void)      { return _ShmSetup; }
  static int      ShmAlloc(void)      { return _ShmAlloc; }
  static uint64_t ShmAllocBytes(void) { return _ShmAllocBytes; }
  static uint64_t      MAX_MPI_SHM_BYTES;
  static int           Hugepages;

  static std::vector<void *> WorldShmCommBufs;
#ifndef ACCELERATOR_AWARE_MPI
  static void *HostCommBuf;
#endif
  static Grid_MPI_Comm WorldComm;
  static int           WorldRank;
  static int           WorldSize;

  static Grid_MPI_Comm WorldShmComm;
  static int           WorldShmRank;
  static int           WorldShmSize;

  static int           WorldNodes;
  static int           WorldNode;

  static std::vector<int>  WorldShmRanks;

  //////////////////////////////////////////////////////////////////////////////////////
  // Create an optimal reordered communicator that makes MPI_Cart_create get it right
  //////////////////////////////////////////////////////////////////////////////////////
  static void Init(Grid_MPI_Comm comm); // Typically MPI_COMM_WORLD
  // Turns MPI_COMM_WORLD into right layout for Cartesian
  static void OptimalCommunicator            (const Coordinate &processors,Grid_MPI_Comm & optimal_comm,Coordinate &ShmDims); 
  static void OptimalCommunicatorHypercube   (const Coordinate &processors,Grid_MPI_Comm & optimal_comm,Coordinate &ShmDims); 
  static void OptimalCommunicatorSharedMemory(const Coordinate &processors,Grid_MPI_Comm & optimal_comm,Coordinate &ShmDims); 
  static void GetShmDims(const Coordinate &WorldDims,Coordinate &ShmDims);
  ///////////////////////////////////////////////////
  // Provide shared memory facilities off comm world
  ///////////////////////////////////////////////////
  static void SharedMemoryAllocate(uint64_t bytes, int flags);
  static void SharedMemoryFree(void);
  static void SharedMemoryCopy(void *dest,void *src,size_t bytes);
  static void SharedMemoryZero(void *dest,size_t bytes);

};

//////////////////////////////
// one per communicator
//////////////////////////////
class SharedMemory 
{
private:
  static const int     MAXLOG2RANKSPERNODE = 16;            

  size_t heap_top;
  size_t heap_bytes;
  size_t heap_size;

#ifndef ACCELERATOR_AWARE_MPI
  size_t host_heap_top;  // set in free all
  size_t host_heap_bytes;// set in free all
  void *HostCommBuf;     // set in SetCommunicator
  size_t host_heap_size; // set in SetCommunicator
#endif
  
protected:

  Grid_MPI_Comm    ShmComm; // for barriers
  int    ShmRank; 
  int    ShmSize;
  std::vector<void *> ShmCommBufs;
  std::vector<int>    ShmRanks;// Mapping comm ranks to Shm ranks

public:
  SharedMemory() {};
  ~SharedMemory();
  ///////////////////////////////////////////////////////////////////////////////////////
  // set the buffers & sizes
  ///////////////////////////////////////////////////////////////////////////////////////
  void SetCommunicator(Grid_MPI_Comm comm);

  ////////////////////////////////////////////////////////////////////////
  // For this instance ; disjoint buffer sets between splits if split grid
  ////////////////////////////////////////////////////////////////////////
  void ShmBarrier(void); 

  ///////////////////////////////////////////////////
  // Call on any instance
  ///////////////////////////////////////////////////
  void SharedMemoryTest(void);
  
  void *ShmBufferSelf(void);
  void *ShmBuffer    (int rank);
  void *ShmBufferTranslate(int rank,void * local_p);
  void *ShmBufferMalloc(size_t bytes);
  void  ShmBufferFreeAll(void) ;
#ifndef ACCELERATOR_AWARE_MPI
  void *HostBufferMalloc(size_t bytes);
  void HostBufferFreeAll(void);
#endif  
  //////////////////////////////////////////////////////////////////////////
  // Make info on Nodes & ranks and Shared memory available
  //////////////////////////////////////////////////////////////////////////
  int NodeCount(void) { return GlobalSharedMemory::WorldNodes;};
  int RankCount(void) { return GlobalSharedMemory::WorldSize;};

};

NAMESPACE_END(Grid);

