/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/communicator/Communicator_base.h

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
#ifndef GRID_COMMUNICATOR_BASE_H
#define GRID_COMMUNICATOR_BASE_H

///////////////////////////////////
// Processor layout information
///////////////////////////////////
#include <Grid/communicator/SharedMemory.h>

NAMESPACE_BEGIN(Grid);

extern bool Stencil_force_mpi ;

class CartesianCommunicator : public SharedMemory {

public:    

  ////////////////////////////////////////////
  // Policies
  ////////////////////////////////////////////
  enum CommunicatorPolicy_t { CommunicatorPolicyConcurrent, CommunicatorPolicySequential };
  static CommunicatorPolicy_t CommunicatorPolicy;
  static void SetCommunicatorPolicy(CommunicatorPolicy_t policy ) { CommunicatorPolicy = policy; }
  static int       nCommThreads;

  ////////////////////////////////////////////
  // Communicator should know nothing of the physics grid, only processor grid.
  ////////////////////////////////////////////
  int              _Nprocessors;     // How many in all
  int              _processor;       // linear processor rank
  unsigned long    _ndimension;
  Coordinate _shm_processors;  // Which dimensions get relayed out over processors lanes.
  Coordinate _processors;      // Which dimensions get relayed out over processors lanes.
  Coordinate _processor_coor;  // linear processor coordinate
  static Grid_MPI_Comm      communicator_world;
  Grid_MPI_Comm             communicator;
  std::vector<Grid_MPI_Comm> communicator_halo;
  
  ////////////////////////////////////////////////
  // Must call in Grid startup
  ////////////////////////////////////////////////
  static void Init(int *argc, char ***argv);

  ////////////////////////////////////////////////
  // Constructors to sub-divide a parent communicator
  // and default to comm world
  ////////////////////////////////////////////////
  CartesianCommunicator(const Coordinate &processors,const CartesianCommunicator &parent,int &srank);
  CartesianCommunicator(const Coordinate &pdimensions_in);
  virtual ~CartesianCommunicator();

private:

  ////////////////////////////////////////////////
  // Private initialise from an MPI communicator
  // Can use after an MPI_Comm_split, but hidden from user so private
  ////////////////////////////////////////////////
  void InitFromMPICommunicator(const Coordinate &processors, Grid_MPI_Comm communicator_base);

public:
  
  
  ////////////////////////////////////////////////////////////////////////////////////////
  // Wraps MPI_Cart routines, or implements equivalent on other impls
  ////////////////////////////////////////////////////////////////////////////////////////
  void ShiftedRanks(int dim,int shift,int & source, int & dest);
  int  RankFromProcessorCoor(Coordinate &coor);
  void ProcessorCoorFromRank(int rank,Coordinate &coor);
  
  int                      Dimensions(void)        ;
  int                      IsBoss(void)            ;
  int                      BossRank(void)          ;
  int                      ThisRank(void)          ;
  const Coordinate & ThisProcessorCoor(void) ;
  const Coordinate & ShmGrid(void)  { return _shm_processors; }  ;
  const Coordinate & ProcessorGrid(void)     ;
  int                ProcessorCount(void)    ;

  ////////////////////////////////////////////////////////////////////////////////
  // very VERY rarely (Log, serial RNG) we need world without a grid
  ////////////////////////////////////////////////////////////////////////////////
  static int  RankWorld(void) ;
  static void BroadcastWorld(int root,void* data, int bytes);
  static void BarrierWorld(void);
  
  ////////////////////////////////////////////////////////////
  // Reduction
  ////////////////////////////////////////////////////////////
  void GlobalMax(RealD &);
  void GlobalMax(RealF &);
  void GlobalSum(RealF &);
  void GlobalSumVector(RealF *,int N);
  void GlobalSum(RealD &);
  void GlobalSumVector(RealD *,int N);
  void GlobalSum(uint32_t &);
  void GlobalSum(uint64_t &);
  void GlobalSumVector(uint64_t*,int N);
  void GlobalSum(ComplexF &c);
  void GlobalSumVector(ComplexF *c,int N);
  void GlobalSum(ComplexD &c);
  void GlobalSumVector(ComplexD *c,int N);
  void GlobalXOR(uint32_t &);
  void GlobalXOR(uint64_t &);

  template<class obj> void GlobalSumP2P(obj &o)
  {
    std::vector<obj> column;
    obj accum = o;
    int source,dest;
    for(int d=0;d<_ndimension;d++){
      column.resize(_processors[d]);
      column[0] = accum;
      std::vector<CommsRequest_t> list;
      for(int p=1;p<_processors[d];p++){
	ShiftedRanks(d,p,source,dest);
	SendToRecvFromBegin(list,
			    &column[0],
			    dest,
			    &column[p],
			    source,
			    sizeof(obj),d*100+p);

      }
      CommsComplete(list);
      for(int p=1;p<_processors[d];p++){
	accum = accum + column[p];
      }
    }
    Broadcast(0,accum);
    o=accum;
  }

  template<class obj> void GlobalSum(obj &o){
    typedef typename obj::scalar_type scalar_type;
    int words = sizeof(obj)/sizeof(scalar_type);
    scalar_type * ptr = (scalar_type *)& o; // Safe alias 
    GlobalSumVector(ptr,words);
  }
  
  ////////////////////////////////////////////////////////////
  // Face exchange, buffer swap in translational invariant way
  ////////////////////////////////////////////////////////////
  void CommsComplete(std::vector<CommsRequest_t> &list);
  void SendToRecvFromBegin(std::vector<CommsRequest_t> &list,
			   void *xmit,
			   int dest,
			   void *recv,
			   int from,
			   int bytes,int dir);
  
  void SendToRecvFrom(void *xmit,
		      int xmit_to_rank,
		      void *recv,
		      int recv_from_rank,
		      int bytes);
  
  double StencilSendToRecvFrom(void *xmit,
			       int xmit_to_rank,int do_xmit,
			       void *recv,
			       int recv_from_rank,int do_recv,
			       int bytes,int dir);

  double StencilSendToRecvFromBegin(std::vector<CommsRequest_t> &list,
				    void *xmit,
				    int xmit_to_rank,int do_xmit,
				    void *recv,
				    int recv_from_rank,int do_recv,
				    int xbytes,int rbytes,int dir);
  
  
  void StencilSendToRecvFromComplete(std::vector<CommsRequest_t> &waitall,int i);
  void StencilBarrier(void);

  ////////////////////////////////////////////////////////////
  // Barrier
  ////////////////////////////////////////////////////////////
  void Barrier(void);
  
  ////////////////////////////////////////////////////////////
  // Broadcast a buffer and composite larger
  ////////////////////////////////////////////////////////////
  void Broadcast(int root,void* data, int bytes);

  ////////////////////////////////////////////////////////////
  // All2All down one dimension
  ////////////////////////////////////////////////////////////
  template<class T> void AllToAll(int dim,std::vector<T> &in, std::vector<T> &out){
    assert(dim>=0);
    assert(dim<_ndimension);
    assert(in.size()==out.size());
    int numnode = _processors[dim];
    uint64_t bytes=sizeof(T);
    uint64_t words=in.size()/numnode;
    assert(numnode * words == in.size());
    assert(words < (1ULL<<31));
    AllToAll(dim,(void *)&in[0],(void *)&out[0],words,bytes);
  }
  void AllToAll(int dim  ,void *in,void *out,uint64_t words,uint64_t bytes);
  void AllToAll(void  *in,void *out,uint64_t words         ,uint64_t bytes);
  
  template<class obj> void Broadcast(int root,obj &data)
  {
    Broadcast(root,(void *)&data,sizeof(data));
  }

}; 

NAMESPACE_END(Grid);

#endif
