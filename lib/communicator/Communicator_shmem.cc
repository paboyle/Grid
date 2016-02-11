    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/communicator/Communicator_shmem.cc

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
#include <mpp/shmem.h>

namespace Grid {

  // Should error check all MPI calls.

CartesianCommunicator::CartesianCommunicator(const std::vector<int> &processors)
{
  _ndimension = processors.size();
  std::vector<int> periodic(_ndimension,1);

  _Nprocessors=1;
  _processors = processors;
  _processor_coor.resize(_ndimension);

  //  shmem_init_thread(SHMEM_THREAD_FUNNELED);
  start_pes(0);
  _processor = shmem_my_pe();
  
  Lexicographic::CoorFromIndex(_processor_coor,_processor,_processors);

  for(int i=0;i<_ndimension;i++){
    _Nprocessors*=_processors[i];
  }
  if ( _processor == 0 ) {
    printf("I'm running SHMEM communications %d  \n",_processor);
  }
  int Size = shmem_n_pes(); 
  assert(Size==_Nprocessors);
}

void CartesianCommunicator::GlobalSum(uint32_t &u){
  static long long source = (long long) u;
  static long long dest   = 0 ;
  static long long llwrk[_SHMEM_REDUCE_MIN_WRKDATA_SIZE];
  static long      psync[_SHMEM_REDUCE_SYNC_SIZE];

  //  int nreduce=1;
  //  int pestart=0;
  //  int logStride=0;
  shmem_longlong_sum_to_all(&dest,&source,1,0,0,_Nprocessors,llwrk,psync);
  u = dest;
}
void CartesianCommunicator::GlobalSum(float &f){
  static float source = f;
  static float dest   = 0 ;
  static float llwrk[_SHMEM_REDUCE_MIN_WRKDATA_SIZE];
  static long  psync[_SHMEM_REDUCE_SYNC_SIZE];

  shmem_float_sum_to_all(&dest,&source,1,0,0,_Nprocessors,llwrk,psync);
  f = dest;
}
void CartesianCommunicator::GlobalSumVector(float *f,int N)
{
  static float source ;
  static float dest   = 0 ;
  static float llwrk[_SHMEM_REDUCE_MIN_WRKDATA_SIZE];
  static long  psync[_SHMEM_REDUCE_SYNC_SIZE];

  for(int i=0;i<N;i++){
    source = f[i];
    shmem_float_sum_to_all(&dest,&source,1,0,0,_Nprocessors,llwrk,psync);
    f[i] = dest;
  }
}
void CartesianCommunicator::GlobalSum(double &d)
{
  static double source = d;
  static double dest   = 0 ;
  static double llwrk[_SHMEM_REDUCE_MIN_WRKDATA_SIZE];
  static long  psync[_SHMEM_REDUCE_SYNC_SIZE];

  shmem_double_sum_to_all(&dest,&source,1,0,0,_Nprocessors,llwrk,psync);
  d = dest;
}
void CartesianCommunicator::GlobalSumVector(double *d,int N)
{
  static double source ;
  static double dest   = 0 ;
  static double llwrk[_SHMEM_REDUCE_MIN_WRKDATA_SIZE];
  static long  psync[_SHMEM_REDUCE_SYNC_SIZE];

  for(int i=0;i<N;i++){
    source = d[i];
    shmem_double_sum_to_all(&dest,&source,1,0,0,_Nprocessors,llwrk,psync);
    d[i] = dest;
  }
}
void CartesianCommunicator::ShiftedRanks(int dim,int shift,int &source,int &dest)
{
  std::vector<int> coor = _processor_coor;
  assert(std::abs(shift) <_processors[dim]);

  coor[dim] = (coor[dim] + shift + _processors[dim])%_processors[dim];
  Lexicographic::IndexFromCoor(coor,source,_processors);

  coor[dim] = (coor[dim] - shift + _processors[dim])%_processors[dim];
  Lexicographic::IndexFromCoor(coor,dest,_processors);

}
int CartesianCommunicator::RankFromProcessorCoor(std::vector<int> &coor)
{
  int rank;
  Lexicographic::IndexFromCoor(coor,rank,_processors);
  return rank;
}
void  CartesianCommunicator::ProcessorCoorFromRank(int rank, std::vector<int> &coor)
{
  Lexicographic::CoorFromIndex(coor,rank,_processors);
}

// Basic Halo comms primitive
void CartesianCommunicator::SendToRecvFrom(void *xmit,
					   int dest,
					   void *recv,
					   int from,
					   int bytes)
{
  std::vector<CommsRequest_t> reqs(0);
  SendToRecvFromBegin(reqs,xmit,dest,recv,from,bytes);
  SendToRecvFromComplete(reqs);
}
void CartesianCommunicator::RecvFrom(void *recv,
				     int from,
				     int bytes) 
{
  // Need to change interface to know send buffer; change to a get/put interface.
  assert(0);
}
void CartesianCommunicator::SendTo(void *xmit,
				   int dest,
				   int bytes)
{
  // Need to change interface to know destination buffer... likely needed for I/O
  assert(0);
}

// Basic Halo comms primitive
void CartesianCommunicator::SendToRecvFromBegin(std::vector<CommsRequest_t> &list,
						void *xmit,
						int dest,
						void *recv,
						int from,
						int bytes)
{
  shmem_putmem(recv,xmit,bytes,dest);
}
void CartesianCommunicator::SendToRecvFromComplete(std::vector<CommsRequest_t> &list)
{
  //  shmem_quiet();      // I'm done
  shmem_barrier_all();// He's done too
}
void CartesianCommunicator::Barrier(void)
{
  shmem_barrier_all();
}
void CartesianCommunicator::Broadcast(int root,void* data, int bytes)
{
  static long  psync[_SHMEM_REDUCE_SYNC_SIZE];
  assert( (bytes % 4)==0);
  shmem_broadcast32(data,data,bytes/4,root,0,0,_Nprocessors,psync);
}
void CartesianCommunicator::BroadcastWorld(int root,void* data, int bytes)
{
  static long  psync[_SHMEM_REDUCE_SYNC_SIZE];
  assert( (bytes % 4)==0);
  shmem_broadcast32(data,data,bytes/4,root,0,0,shmem_n_pes(),psync);
}

}

