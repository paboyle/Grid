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
  srank=0;
  SetCommunicator(communicator_world);
}

CartesianCommunicator::CartesianCommunicator(const Coordinate &processors)
{
  _processors = processors;
  _ndimension = processors.size();  assert(_ndimension>=1);
  _processor_coor.resize(_ndimension);
  
  // Require 1^N processor grid for fake
  _Nprocessors=1;
  _processor = 0;
  for(int d=0;d<_ndimension;d++) {
    assert(_processors[d]==1);
    _processor_coor[d] = 0;
  }
  SetCommunicator(communicator_world);
}

CartesianCommunicator::~CartesianCommunicator(){}

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
					   int bytes)
{
  assert(0);
}
void CartesianCommunicator::AllToAll(int dim,void  *in,void *out,uint64_t words,uint64_t bytes)
{
  bcopy(in,out,bytes*words);
}
void CartesianCommunicator::AllToAll(void  *in,void *out,uint64_t words,uint64_t bytes)
{
  bcopy(in,out,bytes*words);
}

int  CartesianCommunicator::RankWorld(void){return 0;}
void CartesianCommunicator::Barrier(void){}
void CartesianCommunicator::Broadcast(int root,void* data, int bytes) {}
void CartesianCommunicator::BroadcastWorld(int root,void* data, int bytes) { }
int  CartesianCommunicator::RankFromProcessorCoor(Coordinate &coor) {  return 0;}
void CartesianCommunicator::ProcessorCoorFromRank(int rank, Coordinate &coor){  coor = _processor_coor; }
void CartesianCommunicator::ShiftedRanks(int dim,int shift,int &source,int &dest)
{
  source =0;
  dest=0;
}

double CartesianCommunicator::StencilSendToRecvFrom( void *xmit,
						     int xmit_to_rank,
						     void *recv,
						     int recv_from_rank,
						     int bytes, int dir)
{
  return 2.0*bytes;
}
double CartesianCommunicator::StencilSendToRecvFromBegin(std::vector<CommsRequest_t> &list,
							 void *xmit,
							 int xmit_to_rank,
							 void *recv,
							 int recv_from_rank,
							 int bytes, int dir)
{
  return 2.0*bytes;
}
void CartesianCommunicator::StencilSendToRecvFromComplete(std::vector<CommsRequest_t> &waitall,int dir)
{
}

void CartesianCommunicator::StencilBarrier(void){};

NAMESPACE_END(Grid);


