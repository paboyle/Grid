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
#ifdef GRID_COMMS_MPI
#include <mpi.h>
#endif
namespace Grid {
class CartesianCommunicator {
  public:    

  // Communicator should know nothing of the physics grid, only processor grid.

    int              _Nprocessors;     // How many in all
    std::vector<int> _processors;      // Which dimensions get relayed out over processors lanes.
    int              _processor;       // linear processor rank
    std::vector<int> _processor_coor;  // linear processor coordinate
    unsigned long _ndimension;

#ifdef GRID_COMMS_MPI
    MPI_Comm communicator;
    typedef MPI_Request CommsRequest_t;
#else 
    typedef int CommsRequest_t;
#endif

    // Constructor
    CartesianCommunicator(const std::vector<int> &pdimensions_in);

    // Wraps MPI_Cart routines
    void ShiftedRanks(int dim,int shift,int & source, int & dest);
    int  RankFromProcessorCoor(std::vector<int> &coor);
    void ProcessorCoorFromRank(int rank,std::vector<int> &coor);

    /////////////////////////////////
    // Grid information queries
    /////////////////////////////////
    int                      IsBoss(void)            { return _processor==0; };
    int                      BossRank(void)          { return 0; };
    int                      ThisRank(void)          { return _processor; };
    const std::vector<int> & ThisProcessorCoor(void) { return _processor_coor; };
    const std::vector<int> & ProcessorGrid(void)     { return _processors; };
    int                      ProcessorCount(void)    { return _Nprocessors; };

    ////////////////////////////////////////////////////////////
    // Reduction
    ////////////////////////////////////////////////////////////
    void GlobalSum(RealF &);
    void GlobalSumVector(RealF *,int N);

    void GlobalSum(RealD &);
    void GlobalSumVector(RealD *,int N);

    void GlobalSum(uint32_t &);

    void GlobalSum(ComplexF &c)
    {
      GlobalSumVector((float *)&c,2);
    }
    void GlobalSumVector(ComplexF *c,int N)
    {
      GlobalSumVector((float *)c,2*N);
    }

    void GlobalSum(ComplexD &c)
    {
      GlobalSumVector((double *)&c,2);
    }
    void GlobalSumVector(ComplexD *c,int N)
    {
      GlobalSumVector((double *)c,2*N);
    }
    
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

    void RecvFrom(void *recv,
		  int recv_from_rank,
		  int bytes);
    void SendTo(void *xmit,
		int xmit_to_rank,
		int bytes);

    void SendToRecvFromBegin(std::vector<CommsRequest_t> &list,
			 void *xmit,
			 int xmit_to_rank,
			 void *recv,
			 int recv_from_rank,
			 int bytes);
    void SendToRecvFromComplete(std::vector<CommsRequest_t> &waitall);

    ////////////////////////////////////////////////////////////
    // Barrier
    ////////////////////////////////////////////////////////////
    void Barrier(void);

    ////////////////////////////////////////////////////////////
    // Broadcast a buffer and composite larger
    ////////////////////////////////////////////////////////////
    void Broadcast(int root,void* data, int bytes);
    template<class obj> void Broadcast(int root,obj &data)
    {
      Broadcast(root,(void *)&data,sizeof(data));
    };

    static void BroadcastWorld(int root,void* data, int bytes);

}; 
}

#endif
