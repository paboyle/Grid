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
#define SHMEM_VET(addr) 

#define SHMEM_VET_DEBUG(addr) {				\
  if ( ! shmem_addr_accessible(addr,_processor) ) {\
    std::fprintf(stderr,"%d Inaccessible shmem address %lx %s %s\n",_processor,addr,__FUNCTION__,#addr); \
    BACKTRACEFILE();		   \
  }\
}
  int Rank(void) {
    return shmem_my_pe();
  }
void CartesianCommunicator::Init(int *argc, char ***argv) {
  shmem_init();
}
CartesianCommunicator::CartesianCommunicator(const std::vector<int> &processors)
{
  _ndimension = processors.size();
  std::vector<int> periodic(_ndimension,1);

  _Nprocessors=1;
  _processors = processors;
  _processor_coor.resize(_ndimension);

  _processor = shmem_my_pe();
  
  Lexicographic::CoorFromIndex(_processor_coor,_processor,_processors);

  for(int i=0;i<_ndimension;i++){
    _Nprocessors*=_processors[i];
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

  // Inefficient, but don't want to dynamic alloc
  if ( shmem_addr_accessible(f,_processor)  ){
    shmem_float_sum_to_all(f,f,N,0,0,_Nprocessors,llwrk,psync);
    return;
  }

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

  if ( shmem_addr_accessible(d,_processor)  ){
    shmem_double_sum_to_all(d,d,N,0,0,_Nprocessors,llwrk,psync);
    return;
  }

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

  coor[dim] = (_processor_coor[dim] + shift + _processors[dim])%_processors[dim];
  Lexicographic::IndexFromCoor(coor,source,_processors);

  coor[dim] = (_processor_coor[dim] - shift + _processors[dim])%_processors[dim];
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
  SHMEM_VET(xmit);
  SHMEM_VET(recv);
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
  SHMEM_VET(xmit);
  SHMEM_VET(recv);
  //  shmem_putmem_nb(recv,xmit,bytes,dest,NULL);
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
  static uint32_t word;
  uint32_t *array = (uint32_t *) data;
  assert( (bytes % 4)==0);
  int words = bytes/4;
  
  for(int w=0;w<words;w++){
    word = array[w];
    shmem_broadcast32((void *)&word,(void *)&word,1,root,0,0,shmem_n_pes(),psync);
    if ( shmem_my_pe() != root ) {
      array[w] = word;
    }
    shmem_barrier_all();
  }

}
void CartesianCommunicator::BroadcastWorld(int root,void* data, int bytes)
{
  static long  psync[_SHMEM_REDUCE_SYNC_SIZE];
  static uint32_t word;
  uint32_t *array = (uint32_t *) data;
  assert( (bytes % 4)==0);
  int words = bytes/4;

  for(int w=0;w<words;w++){
    word = array[w];
    shmem_broadcast32((void *)&word,(void *)&word,1,root,0,0,shmem_n_pes(),psync);
    if ( shmem_my_pe() != root ) {
      array[w]= word;
    }
    shmem_barrier_all();
  }
}

}

