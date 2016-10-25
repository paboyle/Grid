    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/communicator/Communicator_none.cc

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
#include "Grid.h"
namespace Grid {

///////////////////////////////////////////////////////////////
// Info that is setup once and indept of cartesian layout
///////////////////////////////////////////////////////////////
int CartesianCommunicator::ShmRank;
int CartesianCommunicator::ShmSize;
int CartesianCommunicator::GroupRank;
int CartesianCommunicator::GroupSize;
int CartesianCommunicator::WorldRank;
int CartesianCommunicator::WorldSize;
int CartesianCommunicator::Slave;
void *              CartesianCommunicator::ShmCommBuf;

/////////////////////////////////
// Alloc, free shmem region
/////////////////////////////////
void *CartesianCommunicator::ShmBufferMalloc(size_t bytes){
  //  bytes = (bytes+sizeof(vRealD))&(~(sizeof(vRealD)-1));// align up bytes
  void *ptr = (void *)heap_top;
  heap_top  += bytes;
  heap_bytes+= bytes;
  assert(heap_bytes < MAX_MPI_SHM_BYTES);
  return ptr;
}
void *CartesianCommunicator::ShmBufferFreeAll(void) { 
  heap_top  =(size_t)ShmBufferSelf();
  heap_bytes=0;
}

/////////////////////////////////
// Grid information queries
/////////////////////////////////
int                      CartesianCommunicator::IsBoss(void)            { return _processor==0; };
int                      CartesianCommunicator::BossRank(void)          { return 0; };
int                      CartesianCommunicator::ThisRank(void)          { return _processor; };
const std::vector<int> & CartesianCommunicator::ThisProcessorCoor(void) { return _processor_coor; };
const std::vector<int> & CartesianCommunicator::ProcessorGrid(void)     { return _processors; };
int                      CartesianCommunicator::ProcessorCount(void)    { return _Nprocessors; };

////////////////////////////////////////////////////////////////////////////////
// very VERY rarely (Log, serial RNG) we need world without a grid
////////////////////////////////////////////////////////////////////////////////
int  CartesianCommunicator::RankWorld(void){ return WorldRank; };
int CartesianCommunicator::Ranks    (void) { return WorldSize; };
int CartesianCommunicator::Nodes    (void) { return GroupSize; };
int CartesianCommunicator::Cores    (void) { return ShmSize;   };
int CartesianCommunicator::NodeRank (void) { return GroupRank; };
int CartesianCommunicator::CoreRank (void) { return ShmRank;   };

void CartesianCommunicator::GlobalSum(ComplexF &c)
{
  GlobalSumVector((float *)&c,2);
}
void CartesianCommunicator::GlobalSumVector(ComplexF *c,int N)
{
  GlobalSumVector((float *)c,2*N);
}
void CartesianCommunicator::GlobalSum(ComplexD &c)
{
  GlobalSumVector((double *)&c,2);
}
void CartesianCommunicator::GlobalSumVector(ComplexD *c,int N)
{
  GlobalSumVector((double *)c,2*N);
}

#ifndef GRID_COMMS_MPI3

void CartesianCommunicator::StencilSendToRecvFromBegin(std::vector<CommsRequest_t> &list,
						       void *xmit,
						       int xmit_to_rank,
						       void *recv,
						       int recv_from_rank,
						       int bytes)
{
  SendToRecvFromBegin(list,xmit,xmit_to_rank,recv,recv_from_rank,bytes);
}
void CartesianCommunicator::StencilSendToRecvFromComplete(std::vector<CommsRequest_t> &waitall)
{
  SendToRecvFromComplete(waitall);
}
void CartesianCommunicator::StencilBarrier(void){};

commVector<uint8_t> CartesianCommunicator::ShmBufStorageVector;

void *CartesianCommunicator::ShmBufferSelf(void) { return ShmCommBuf; }

void *CartesianCommunicator::ShmBuffer(int rank) {
  return NULL;
}
void *CartesianCommunicator::ShmBufferTranslate(int rank,void * local_p) { 
  return NULL;
}
void CartesianCommunicator::ShmInitGeneric(void){
  ShmBufStorageVector.resize(MAX_MPI_SHM_BYTES);
  ShmCommBuf=(void *)&ShmBufStorageVector[0];
}

#endif
  
}

