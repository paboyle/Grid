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
#include <Grid/GridCore.h>

void GridAbort(void) { abort(); }

NAMESPACE_BEGIN(Grid);

///////////////////////////////////////////////////////////////////////////////////////////////////
// Info that is setup once and indept of cartesian layout
///////////////////////////////////////////////////////////////////////////////////////////////////
Grid_MPI_Comm       CartesianCommunicator::communicator_world;


void CartesianCommunicator::Init(int *argc, char *** arv)
{
  GlobalSharedMemory::Init(communicator_world);
  GlobalSharedMemory::SharedMemoryAllocate(
					   GlobalSharedMemory::MAX_MPI_SHM_BYTES,
					   GlobalSharedMemory::Hugepages);
}

CartesianCommunicator::CartesianCommunicator(const Coordinate &processors,const CartesianCommunicator &parent,int &srank) 
  : CartesianCommunicator(processors) 
{
  _shm_processors = Coordinate(processors.size(),1);
  srank=0;
  SetCommunicator(communicator_world);
}

CartesianCommunicator::CartesianCommunicator(const Coordinate &processors)
{
  _shm_processors = Coordinate(processors.size(),1);
  _processors = processors;
  _ndimension = processors.size();  GRID_ASSERT(_ndimension>=1);
  _processor_coor.resize(_ndimension);
  
  // Require 1^N processor grid for fake
  _Nprocessors=1;
  _processor = 0;
  for(int d=0;d<_ndimension;d++) {
    GRID_ASSERT(_processors[d]==1);
    _processor_coor[d] = 0;
  }
  SetCommunicator(communicator_world);
}

CartesianCommunicator::~CartesianCommunicator(){}

void CartesianCommunicator::GlobalMax(float &){}
void CartesianCommunicator::GlobalMax(double &){}
void CartesianCommunicator::GlobalSum(float &){}
void CartesianCommunicator::GlobalSumVector(float *,int N){}
void CartesianCommunicator::GlobalSum(double &){}
void CartesianCommunicator::GlobalSumVector(double *,int N){}
void CartesianCommunicator::GlobalSum(uint32_t &){}
void CartesianCommunicator::GlobalSum(uint64_t &){}
void CartesianCommunicator::GlobalSumVector(uint64_t *,int N){}
void CartesianCommunicator::GlobalXOR(uint32_t &){}
void CartesianCommunicator::GlobalXOR(uint64_t &){}


// Basic Halo comms primitive -- should never call in single node
void CartesianCommunicator::SendToRecvFrom(void *xmit,
					   int dest,
					   void *recv,
					   int from,
					   uint64_t bytes)
{
  GRID_ASSERT(0);
}
void CartesianCommunicator::CommsComplete(std::vector<CommsRequest_t> &list){ GRID_ASSERT(list.size()==0);}
void CartesianCommunicator::SendToRecvFromBegin(std::vector<CommsRequest_t> &list,
						void *xmit,
						int dest,
						void *recv,
						int from,
						uint64_t bytes,int dir)
{
  GRID_ASSERT(0);
}

void CartesianCommunicator::AllToAll(int dim,void  *in,void *out,uint64_t words,uint64_t bytes)
{
  bcopy(in,out,bytes*words);
}
void CartesianCommunicator::AllToAll(void  *in,void *out,uint64_t words,uint64_t bytes)
{
  bcopy(in,out,bytes*words);
}
void CartesianCommunicator::AllToAllV(void *in ,const std::vector<int> &sendcounts,const std::vector<int> &senddispls,
                                      void *out,const std::vector<int> &recvcounts,const std::vector<int> &recvdispls,
                                      uint64_t bytes)
{
  // Single rank: the exchange degenerates to a copy of our own segment
  GRID_ASSERT(sendcounts.size()==1);
  GRID_ASSERT(recvcounts.size()==1);
  GRID_ASSERT(sendcounts[0]==recvcounts[0]);
  bcopy((char *)in +(uint64_t)senddispls[0]*bytes,
        (char *)out+(uint64_t)recvdispls[0]*bytes,bytes*(uint64_t)sendcounts[0]);
}

void CartesianCommunicator::AllGatherV(void *in ,int sendcount,
                                       void *out,const std::vector<int> &recvcounts,const std::vector<int> &recvdispls,
                                       uint64_t bytes)
{
  // Single rank: the gather degenerates to a copy of our own contribution
  GRID_ASSERT(recvcounts.size()==1);
  GRID_ASSERT(recvdispls.size()==1);
  GRID_ASSERT(recvcounts[0]==sendcount);
  bcopy((char *)in,
        (char *)out+(uint64_t)recvdispls[0]*bytes,bytes*(uint64_t)sendcount);
}

void CartesianCommunicator::AllGather(void *in,void *out,uint64_t words,uint64_t bytes)
{
  bcopy((char *)in,(char *)out,bytes*words);
}

int  CartesianCommunicator::RankWorld(void){return 0;}
void CartesianCommunicator::Barrier(void){}
void CartesianCommunicator::Broadcast(int root,void* data, uint64_t bytes) {}
void CartesianCommunicator::BroadcastWorld(int root,void* data, uint64_t bytes) { }
void CartesianCommunicator::BarrierWorld(void) { }
int  CartesianCommunicator::RankFromProcessorCoor(Coordinate &coor) {  return 0;}
void CartesianCommunicator::ProcessorCoorFromRank(int rank, Coordinate &coor){  coor = _processor_coor; }
void CartesianCommunicator::ShiftedRanks(int dim,int shift,int &source,int &dest)
{
  source =0;
  dest=0;
}

int CartesianCommunicator::IsOffNode(int rank) { return false; }

double CartesianCommunicator::StencilSendToRecvFrom( void *xmit,
						     int xmit_to_rank,int dox,
						     void *recv,
						     int recv_from_rank,int dor,
						     uint64_t bytes, int dir)
{
  return 2.0*bytes;
}
void CartesianCommunicator::StencilSendToRecvFromPollIRecv(std::vector<CommsRequest_t> &list) {};
void CartesianCommunicator::StencilSendToRecvFromPollDtoH(std::vector<CommsRequest_t> &list) {};
double CartesianCommunicator::StencilSendToRecvFromPrepare(std::vector<CommsRequest_t> &list,
							   void *xmit,
							   int xmit_to_rank,int dox,
							   void *recv,
							   int recv_from_rank,int dor,
							   uint64_t xbytes,uint64_t rbytes, int dir)
{
  return 0.0;
}
double CartesianCommunicator::StencilSendToRecvFromBegin(std::vector<CommsRequest_t> &list,
							 void *xmit, void *xmit_comp,
							 int xmit_to_rank,int dox,
							 void *recv, void *recv_comp,
							 int recv_from_rank,int dor,
							 uint64_t xbytes,uint64_t rbytes, int dir)
{
  return xbytes+rbytes;
}
void CartesianCommunicator::StencilSendToRecvFromComplete(std::vector<CommsRequest_t> &waitall,int dir)
{
}

void CartesianCommunicator::StencilBarrier(void){};

NAMESPACE_END(Grid);


