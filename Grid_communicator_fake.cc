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

void CartesianCommunicator::GlobalSum(float &){}
void CartesianCommunicator::GlobalSumVector(float *,int N){}
void CartesianCommunicator::GlobalSum(double &){}
void CartesianCommunicator::GlobalSumVector(double *,int N){}

// Basic Halo comms primitive
void CartesianCommunicator::SendToRecvFrom(void *xmit,
					   int dest,
					   void *recv,
					   int from,
					   int bytes)
{
  exit(-1);
}

void CartesianCommunicator::Barrier(void)
{
}

void CartesianCommunicator::Broadcast(int root,void* data, int bytes)
{
}


void CartesianCommunicator::ShiftedRanks(int dim,int shift,int &source,int &dest)
{
  source =1;
  dest=1;
}
int CartesianCommunicator::RankFromProcessorCoor(std::vector<int> &coor)
{
  return 1;
}
void  CartesianCommunicator::ProcessorCoorFromRank(int rank, std::vector<int> &coor)
{
}


}

