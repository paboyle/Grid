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
typedef struct HandShake_t { 
  uint64_t seq_local;
  uint64_t seq_remote;
} HandShake;

static Vector< HandShake > XConnections;
static Vector< HandShake > RConnections;

void CartesianCommunicator::Init(int *argc, char ***argv) {
  shmem_init();
  XConnections.resize(shmem_n_pes());
  RConnections.resize(shmem_n_pes());
  for(int pe =0 ; pe<shmem_n_pes();pe++){
    XConnections[pe].seq_local = 0;
    XConnections[pe].seq_remote= 0;
    RConnections[pe].seq_local = 0;
    RConnections[pe].seq_remote= 0;
  }
  shmem_barrier_all();
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
  static long long source ;
  static long long dest   ;
  static long long llwrk[_SHMEM_REDUCE_MIN_WRKDATA_SIZE];
  static long      psync[_SHMEM_REDUCE_SYNC_SIZE];

  //  int nreduce=1;
  //  int pestart=0;
  //  int logStride=0;

  source = u;
  dest   = 0;
  shmem_longlong_sum_to_all(&dest,&source,1,0,0,_Nprocessors,llwrk,psync);
  shmem_barrier_all(); // necessary?
  u = dest;
}
void CartesianCommunicator::GlobalSum(uint64_t &u){
  static long long source ;
  static long long dest   ;
  static long long llwrk[_SHMEM_REDUCE_MIN_WRKDATA_SIZE];
  static long      psync[_SHMEM_REDUCE_SYNC_SIZE];

  //  int nreduce=1;
  //  int pestart=0;
  //  int logStride=0;

  source = u;
  dest   = 0;
  shmem_longlong_sum_to_all(&dest,&source,1,0,0,_Nprocessors,llwrk,psync);
  shmem_barrier_all(); // necessary?
  u = dest;
}
void CartesianCommunicator::GlobalSum(float &f){
  static float source ;
  static float dest   ;
  static float llwrk[_SHMEM_REDUCE_MIN_WRKDATA_SIZE];
  static long  psync[_SHMEM_REDUCE_SYNC_SIZE];

  source = f;
  dest   =0.0;
  shmem_float_sum_to_all(&dest,&source,1,0,0,_Nprocessors,llwrk,psync);
  f = dest;
}
void CartesianCommunicator::GlobalSumVector(float *f,int N)
{
  static float source ;
  static float dest   = 0 ;
  static float llwrk[_SHMEM_REDUCE_MIN_WRKDATA_SIZE];
  static long  psync[_SHMEM_REDUCE_SYNC_SIZE];

  if ( shmem_addr_accessible(f,_processor)  ){
    shmem_float_sum_to_all(f,f,N,0,0,_Nprocessors,llwrk,psync);
    return;
  }

  for(int i=0;i<N;i++){
    dest   =0.0;
    source = f[i];
    shmem_float_sum_to_all(&dest,&source,1,0,0,_Nprocessors,llwrk,psync);
    f[i] = dest;
  }
}
void CartesianCommunicator::GlobalSum(double &d)
{
  static double source;
  static double dest  ;
  static double llwrk[_SHMEM_REDUCE_MIN_WRKDATA_SIZE];
  static long  psync[_SHMEM_REDUCE_SYNC_SIZE];

  source = d;
  dest   = 0;
  shmem_double_sum_to_all(&dest,&source,1,0,0,_Nprocessors,llwrk,psync);
  d = dest;
}
void CartesianCommunicator::GlobalSumVector(double *d,int N)
{
  static double source ;
  static double dest   ;
  static double llwrk[_SHMEM_REDUCE_MIN_WRKDATA_SIZE];
  static long  psync[_SHMEM_REDUCE_SYNC_SIZE];

  if ( shmem_addr_accessible(d,_processor)  ){
    shmem_double_sum_to_all(d,d,N,0,0,_Nprocessors,llwrk,psync);
    return;
  }

  for(int i=0;i<N;i++){
    source = d[i];
    dest   =0.0;
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

void CartesianCommunicator::SendRecvPacket(void *xmit,
					   void *recv,
					   int sender,
					   int receiver,
					   int bytes)
{
  static uint64_t seq;

  assert(recv!=xmit);
  volatile HandShake *RecvSeq = (volatile HandShake *) & RConnections[sender];
  volatile HandShake *SendSeq = (volatile HandShake *) & XConnections[receiver];

  if ( _processor == sender ) {

    printf("Sender SHMEM pt2pt %d -> %d\n",sender,receiver);
    // Check he has posted a receive
    while(SendSeq->seq_remote == SendSeq->seq_local);

    printf("Sender receive %d posted\n",sender,receiver);

    // Advance our send count
    seq = ++(SendSeq->seq_local);
    
    // Send this packet 
    SHMEM_VET(recv);
    shmem_putmem(recv,xmit,bytes,receiver);
    shmem_fence();

    printf("Sender sent payload %d\n",seq);
    //Notify him we're done
    shmem_putmem((void *)&(RecvSeq->seq_remote),&seq,sizeof(seq),receiver);
    shmem_fence();
    printf("Sender ringing door bell  %d\n",seq);
  }
  if ( _processor == receiver ) {

    printf("Receiver SHMEM pt2pt %d->%d\n",sender,receiver);
    // Post a receive
    seq = ++(RecvSeq->seq_local);
    shmem_putmem((void *)&(SendSeq->seq_remote),&seq,sizeof(seq),sender);

    printf("Receiver Opening letter box %d\n",seq);

    
    // Now wait until he has advanced our reception counter
    while(RecvSeq->seq_remote != RecvSeq->seq_local);

    printf("Receiver Got the mail %d\n",seq);
  }
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

  if ( shmem_addr_accessible(data,_processor)  ){
    shmem_broadcast32(data,data,words,root,0,0,shmem_n_pes(),psync);
    return;
  }

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

