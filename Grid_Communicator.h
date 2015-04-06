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

    // Constructor
    CartesianCommunicator(std::vector<int> &pdimensions_in);

    // Wraps MPI_Cart routines
    void ShiftedRanks(int dim,int shift,int & source, int & dest);
    int  RankFromProcessorCoor(std::vector<int> &coor);
    void ProcessorCoorFromRank(int rank,std::vector<int> &coor);

    /////////////////////////////////
    // Grid information queries
    /////////////////////////////////
    int                      IsBoss(void)            { return _processor==0; };
    int                      ThisRank(void)          { return _processor; };
    const std::vector<int> & ThisProcessorCoor(void) { return _processor_coor; };
    const std::vector<int> & ProcessorGrid(void)     { return _processors; };
    int                      ProcessorCount(void)    { return _Nprocessors; };

    ////////////////////////////////////////////////////////////
    // Reduction
    ////////////////////////////////////////////////////////////
    void GlobalSum(float &);
    void GlobalSumVector(float *,int N);

    void GlobalSum(double &);
    void GlobalSumVector(double *,int N);
    
    template<class obj> void GlobalSumObj(obj &o){

      typedef typename obj::scalar_type scalar_type;
      int words = sizeof(obj)/sizeof(scalar_type);

      scalar_type * ptr = (scalar_type *)& o;
      GlobalSum(ptr,words);
    }
    ////////////////////////////////////////////////////////////
    // Face exchange
    ////////////////////////////////////////////////////////////
    void SendToRecvFrom(void *xmit,
			int xmit_to_rank,
			void *recv,
			int recv_from_rank,
			int bytes);

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

}; 
}

#endif
