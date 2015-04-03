#ifndef GRID_COMMUNICATOR_H
#define GRID_COMMUNICATOR_H
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
#endif

    CartesianCommunicator(std::vector<int> &pdimensions_in);

    void ShiftedRanks(int dim,int shift,int & source, int & dest);
    int  Rank(std::vector<int> coor);
    int  MyRank(void);
    void GlobalSumF(float &);
    void GlobalSumFVector(float *,int N);

    void GlobalSumF(double &);
    void GlobalSumFVector(double *,int N);

    void SendToRecvFrom(void *xmit,
			int xmit_to_rank,
			void *recv,
			int recv_from_rank,
			int bytes);
    
}; 
}

#endif
