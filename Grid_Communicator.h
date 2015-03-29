#ifndef GRID_COMMUNICATOR_H
#define GRID_COMMUNICATOR_H
///////////////////////////////////
// Processor layout information
///////////////////////////////////
namespace dpo {
class CartesianCommunicator {
  public:    

  // Communicator should know nothing of the physics grid, only processor grid.

    std::vector<int> _processors;      // Which dimensions get relayed out over processors lanes.
    int              _processor;       // linear processor rank
    std::vector<int> _processor_coor;  // linear processor coordinate
    unsigned long _ndimension;

#ifdef GRID_COMMS_MPI
    MPI_Comm communicator;
#endif

    CartesianCommunicator(std::vector<int> &pdimensions_in);

    int  Rank(std::vector<int> coor);
    void GlobalSumF(float &);
    void GlobalSumFVector(float *,int N);

    void GlobalSumF(double &);
    void GlobalSumFVector(double *,int N);

    void SendToRecvFrom(void *xmit,
			std::vector<int> to_coordinate,
			void *recv,
			std::vector<int> from_coordinate,
			int bytes);
    
}; 
}

#endif
