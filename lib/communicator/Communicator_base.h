
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

namespace Grid {

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
  std::vector<int> _processors;      // Which dimensions get relayed out over processors lanes.
  int              _processor;       // linear processor rank
  std::vector<int> _processor_coor;  // linear processor coordinate
  unsigned long    _ndimension;
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
  CartesianCommunicator(const std::vector<int> &processors,const CartesianCommunicator &parent,int &srank);
  CartesianCommunicator(const std::vector<int> &pdimensions_in);
  virtual ~CartesianCommunicator();

 private:

  ////////////////////////////////////////////////
  // Private initialise from an MPI communicator
  // Can use after an MPI_Comm_split, but hidden from user so private
  ////////////////////////////////////////////////
  void InitFromMPICommunicator(const std::vector<int> &processors, Grid_MPI_Comm communicator_base);

 public:

  
  ////////////////////////////////////////////////////////////////////////////////////////
  // Wraps MPI_Cart routines, or implements equivalent on other impls
  ////////////////////////////////////////////////////////////////////////////////////////
  void ShiftedRanks(int dim,int shift,int & source, int & dest);
  int  RankFromProcessorCoor(std::vector<int> &coor);
  void ProcessorCoorFromRank(int rank,std::vector<int> &coor);
  
  int                      Dimensions(void)        ;
  int                      IsBoss(void)            ;
  int                      BossRank(void)          ;
  int                      ThisRank(void)          ;
  const std::vector<int> & ThisProcessorCoor(void) ;
  const std::vector<int> & ProcessorGrid(void)     ;
  int                      ProcessorCount(void)    ;

  ////////////////////////////////////////////////////////////////////////////////
  // very VERY rarely (Log, serial RNG) we need world without a grid
  ////////////////////////////////////////////////////////////////////////////////
  static int  RankWorld(void) ;
  static void BroadcastWorld(int root,void* data, int bytes);
  
  ////////////////////////////////////////////////////////////
  // Reduction
  ////////////////////////////////////////////////////////////
  void GlobalSum(RealF &);
  void GlobalSumVector(RealF *,int N);
  void GlobalSum(RealD &);
  void GlobalSumVector(RealD *,int N);
  void GlobalSum(uint32_t &);
  void GlobalSum(uint64_t &);
  void GlobalSum(ComplexF &c);
  void GlobalSumVector(ComplexF *c,int N);
  void GlobalSum(ComplexD &c);
  void GlobalSumVector(ComplexD *c,int N);
  void GlobalXOR(uint32_t &);
  void GlobalXOR(uint64_t &);
  
  template<class obj> void GlobalSum(obj &o){
    typedef typename obj::scalar_type scalar_type;
    int words = sizeof(obj)/sizeof(scalar_type);
    scalar_type * ptr = (scalar_type *)& o;
    GlobalSumVector(ptr,words);
  }
  
  ////////////////////////////////////////////////////////////
  // Face exchange, buffer swap in translational invariant way
  ////////////////////////////////////////////////////////////
  void SendToRecvFrom(void *xmit,
		      int xmit_to_rank,
		      void *recv,
		      int recv_from_rank,
		      int bytes);
  
  void SendRecvPacket(void *xmit,
		      void *recv,
		      int xmit_to_rank,
		      int recv_from_rank,
		      int bytes);
  
  void SendToRecvFromBegin(std::vector<CommsRequest_t> &list,
			   void *xmit,
			   int xmit_to_rank,
			   void *recv,
			   int recv_from_rank,
			   int bytes);
  
  void SendToRecvFromComplete(std::vector<CommsRequest_t> &waitall);

  double StencilSendToRecvFrom(void *xmit,
			       int xmit_to_rank,
			       void *recv,
			       int recv_from_rank,
			       int bytes,int dir);

  double StencilSendToRecvFromBegin(std::vector<CommsRequest_t> &list,
				    void *xmit,
				    int xmit_to_rank,
				    void *recv,
				    int recv_from_rank,
				    int bytes,int dir);
  
  
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
    };

}; 
}

#endif
