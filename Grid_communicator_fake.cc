#include "Grid.h"
namespace Grid {

CartesianCommunicator::CartesianCommunicator(std::vector<int> &processors)
{
  _ndimension = _processors.size();
  _processor_coor.resize(_ndimension);
  _processors = processors;
  
  // Require 1^N processor grid for fake
  for(int d=0;d<_ndimension;d++) if(_processors[d]!=1) exit(-1);

  _processor = 0;// I am the one. The only one..
  for(int d=0;d<_ndimension;d++) _processor_coor[d] = 0;
}

void CartesianCommunicator::GlobalSumF(float &){}
void CartesianCommunicator::GlobalSumFVector(float *,int N){}
void CartesianCommunicator::GlobalSumF(double &){}
void CartesianCommunicator::GlobalSumFVector(double *,int N){}

// Basic Halo comms primitive
void CartesianCommunicator::SendToRecvFrom(void *xmit,
		    std::vector<int> to_coordinate,
		    void *recv,
		    std::vector<int> from_coordinate,
		    int bytes)
{
  exit(-1);
}

}

