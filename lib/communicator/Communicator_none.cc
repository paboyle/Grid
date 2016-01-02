    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/communicator/Communicator_none.cc

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

