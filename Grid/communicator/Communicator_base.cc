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
#include <Grid/GridCore.h>
#include <fcntl.h>
#include <unistd.h>
#include <limits.h>
#include <sys/mman.h>

NAMESPACE_BEGIN(Grid);

bool Stencil_force_mpi = true;

///////////////////////////////////////////////////////////////
// Info that is setup once and indept of cartesian layout
///////////////////////////////////////////////////////////////
CartesianCommunicator::CommunicatorPolicy_t  
CartesianCommunicator::CommunicatorPolicy= CartesianCommunicator::CommunicatorPolicyConcurrent;
int CartesianCommunicator::nCommThreads = -1;
CartesianCommunicator::StencilCompressionPolicy_t  
CartesianCommunicator::StencilCompressionPolicy= CartesianCommunicator::StencilCompressionPolicyNone;

/////////////////////////////////
// Grid information queries
/////////////////////////////////
int                      CartesianCommunicator::Dimensions(void)        { return _ndimension; };
int                      CartesianCommunicator::IsBoss(void)            { return _processor==0; };
int                      CartesianCommunicator::BossRank(void)          { return 0; };
int                      CartesianCommunicator::ThisRank(void)          { return _processor; };
const Coordinate & CartesianCommunicator::ThisProcessorCoor(void) { return _processor_coor; };
const Coordinate & CartesianCommunicator::ProcessorGrid(void)     { return _processors; };
int                      CartesianCommunicator::ProcessorCount(void)    { return _Nprocessors; };

////////////////////////////////////////////////////////////////////////////////
// very VERY rarely (Log, serial RNG) we need world without a grid
////////////////////////////////////////////////////////////////////////////////

void CartesianCommunicator::GlobalSum(ComplexF &c)
{
  GlobalSumVector((float *)&c,2);
}
void CartesianCommunicator::GlobalSumVector(ComplexF *c,int N)
{
  GlobalSumVector((float *)c,2*N);
}
void CartesianCommunicator::GlobalSum(ComplexD &c)
{
  GlobalSumVector((double *)&c,2);
}
void CartesianCommunicator::GlobalSumVector(ComplexD *c,int N)
{
  GlobalSumVector((double *)c,2*N);
}
  
NAMESPACE_END(Grid);


