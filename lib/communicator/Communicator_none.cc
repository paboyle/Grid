#include "Grid.h"
namespace Grid {

CartesianCommunicator::CartesianCommunicator(const std::vector<int> &processors)
{
  _processors = processors;
  _ndimension = processors.size();
  _processor_coor.resize(_ndimension);
  
  // Require 1^N processor grid for fake
  _Nprocessors=1;
  _processor = 0;
  for(int d=0;d<_ndimension;d++) {
    assert(_processors[d]==1);
    _processor_coor[d] = 0;
  }
}

void CartesianCommunicator::GlobalSum(float &){}
void CartesianCommunicator::GlobalSumVector(float *,int N){}
void CartesianCommunicator::GlobalSum(double &){}
void CartesianCommunicator::GlobalSum(uint32_t &){}
void CartesianCommunicator::GlobalSumVector(double *,int N){}

void CartesianCommunicator::RecvFrom(void *recv,
				     int recv_from_rank,
				     int bytes) 
{
  assert(0);
}
void CartesianCommunicator::SendTo(void *xmit,
				   int xmit_to_rank,
				   int bytes)
{
  assert(0);
}


// Basic Halo comms primitive -- should never call in single node
void CartesianCommunicator::SendToRecvFrom(void *xmit,
					   int dest,
					   void *recv,
					   int from,
					   int bytes)
{
  assert(0);
}
void CartesianCommunicator::SendToRecvFromBegin(std::vector<CommsRequest_t> &list,
						void *xmit,
						int dest,
						void *recv,
						int from,
						int bytes)
{
  assert(0);
}
void CartesianCommunicator::SendToRecvFromComplete(std::vector<CommsRequest_t> &list)
{
  assert(0);
}

void CartesianCommunicator::Barrier(void)
{
}

void CartesianCommunicator::Broadcast(int root,void* data, int bytes)
{
}
void CartesianCommunicator::BroadcastWorld(int root,void* data, int bytes)
{
}


void CartesianCommunicator::ShiftedRanks(int dim,int shift,int &source,int &dest)
{
  source =0;
  dest=0;
}
int CartesianCommunicator::RankFromProcessorCoor(std::vector<int> &coor)
{
  return 0;
}
void  CartesianCommunicator::ProcessorCoorFromRank(int rank, std::vector<int> &coor)
{
}


}

