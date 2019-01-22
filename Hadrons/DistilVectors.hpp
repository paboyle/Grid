/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/A2AVectors.hpp

Copyright (C) 2015-2018

Author: Antonin Portelli <antonin.portelli@me.com>
Author: fionnoh <fionnoh@gmail.com>

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
#ifndef Distil_Vectors_hpp_
#define Distil_Vectors_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Environment.hpp>
#include <Hadrons/Solver.hpp>
#include "Grid/lattice/Lattice_peekpoke.h"
#include <Grid/Eigen/unsupported/CXX11/Tensor>

BEGIN_HADRONS_NAMESPACE

template<typename LatticeObj>
class Perambulator : Serializable{
	  // TODO: The next line makes friends across all combinations
	  //     (not much of a problem given all public anyway ...)
	  //      FYI, the bug here was that I forgot that the friend is templated
  template<typename T> friend std::ostream & operator<<(std::ostream &os, const Perambulator<T>& p);
protected:
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS( Perambulator,
                                  std::string,             ID,          // Allows owner to specialise
                                  std::string,             Provenance,  // For info only
                                  std::vector<int>,        dimensions,
                                  std::vector<LatticeObj>, perambulator,
                                  // Following items are redundant, but useful
                                  int,     nd,          // Number of dimensions
                                  size_t,  NumElements); // Number of elements
protected:
  // Constructor common code
  inline void ConstructCommon(const int * Dimensions) {
    assert(nd > 0);
    dimensions.resize(nd);
    NumElements = 1;
    for(int i = 0 ; i < nd ; i++) {
      assert(Dimensions[i] > 0);
      NumElements *= (size_t) Dimensions[i];
      dimensions[i] = Dimensions[i];
    }
    //const LatticeObj perambulatorDefault;
    perambulator.resize(NumElements);//,perambulatorDefault);
  }
public:
  // Constructor with dimensions passed as std::vector<int>
  inline Perambulator(const std::vector<int> & Dimensions)
  : nd {(int) Dimensions.size()} {
    ConstructCommon( &Dimensions[0] ); }
  
  // Constructor with dimensions passed as std::vector<int>
  inline Perambulator(const std::vector<int> & Dimensions, const std::string sID)
  : nd {(int) Dimensions.size()}, ID(sID) {
    ConstructCommon( &Dimensions[0] ); }
  
  // Constructor with dimensions passed as std::vector<int>
  inline Perambulator(const std::vector<int> & Dimensions, const std::string sID, const std::string sProvenance)
  : nd {(int) Dimensions.size()}, ID(sID), Provenance(sProvenance) {
    ConstructCommon( &Dimensions[0] ); }
  
  // Constructor with dimensions passed as individual parameters
  // FYI: The caller is free to ignore the names and use the indices however they see fit
  inline Perambulator(int NumNoise, int NumEvec=1, int NumTime=1, int NumSpin=1, int I_k=1, int I_t=1, int I_s=1) {
    int Dimensions[]={NumNoise,NumEvec,NumTime,NumSpin,I_k,I_t,I_s};
    nd = sizeof(Dimensions)/sizeof(Dimensions[0]);
    while( nd > 1 && Dimensions[nd-1] == 1 )
      nd--;
    ConstructCommon( Dimensions );
  }
  
  inline LatticeObj & operator()(size_t count, const int * Coord) {
    assert( count == nd );
    assert( Coord );
    size_t idx = 0;
    // C memory order (???)
    for( int d = 0 ; d < nd ; d++ ) {
      assert( Coord[d] < dimensions[d] );
      idx *= (size_t) dimensions[d];
      idx += (size_t) Coord[d];
    }
    return perambulator[idx];
  }

  inline LatticeObj & operator()(const std::vector<int> Coord) {
    return operator()(Coord.size(), &Coord[0]);
  }
  
  inline LatticeObj & operator()(int idxNoise, int idxEvec=0, int idxTime=0, int idxSpin=0, int I_k=0, int I_t=0, int I_s=0) {
    int MyIndex[]={idxNoise,idxEvec,idxTime,idxSpin,I_k,I_t,I_s};
    int i = sizeof(MyIndex)/sizeof(MyIndex[0]);
    assert( i >= nd );
    while( i > nd )
      assert(MyIndex[--i] == 0);
    return operator()(i, MyIndex);
  }
};
 /*
#define BEGIN_GRID_NAMESPACE namespace Grid {
BEGIN_GRID_NAMESPACE

void CartesianCommunicatorCandidate::SliceShare( GridBase * gridLowDim, GridBase * gridHighDim, void * Buffer, int BufferSize )
{
  // Work out which dimension is the spread-out dimension
  assert(gridLowDim);
  assert(gridHighDim);
  const int iNumDims{(const int)gridHighDim->_gdimensions.size()};
  assert(iNumDims == gridLowDim->_gdimensions.size());
  int dimSpreadOut = -1;
  std::vector<int> coor(iNumDims);
  for( int i = 0 ; i < iNumDims ; i++ ) {
    coor[i] = gridHighDim->_processor_coor[i];
    if( gridLowDim->_gdimensions[i] != gridHighDim->_gdimensions[i] ) {
      assert( dimSpreadOut == -1 );
      assert( gridLowDim->_processors[i] == 1 ); // easiest assumption to make for now
      dimSpreadOut = i;
    }
  }
  if( dimSpreadOut != -1 && gridHighDim->_processors[dimSpreadOut] != gridLowDim->_processors[dimSpreadOut] ) {
    // Make sure the same number of data elements exist on each slice
    const int NumSlices{gridHighDim->_processors[dimSpreadOut] / gridLowDim->_processors[dimSpreadOut]};
    assert(gridHighDim->_processors[dimSpreadOut] == gridLowDim->_processors[dimSpreadOut] * NumSlices);
    const int SliceSize{BufferSize/NumSlices};
    CCC_DEBUG_DUMP(Buffer, NumSlices, SliceSize);
    assert(BufferSize == SliceSize * NumSlices);
#ifndef USE_LOCAL_SLICES
    assert(0); // Can't do this without MPI (should really test whether MPI is defined)
#else
    const auto MyRank{gridHighDim->ThisRank()};
    std::vector<CommsRequest_t> reqs(0);
    int MySlice{coor[dimSpreadOut]};
    char * const _buffer{(char *)Buffer};
    char * const MyData{_buffer + MySlice * SliceSize};
    for(int i = 1; i < NumSlices ; i++ ){
      int SendSlice = ( MySlice + i ) % NumSlices;
      int RecvSlice = ( MySlice - i + NumSlices ) % NumSlices;
      char * const RecvData{_buffer + RecvSlice * SliceSize};
      coor[dimSpreadOut] = SendSlice;
      const auto SendRank{gridHighDim->RankFromProcessorCoor(coor)};
      coor[dimSpreadOut] = RecvSlice;
      const auto RecvRank{gridHighDim->RankFromProcessorCoor(coor)};
      std::cout << GridLogMessage << "Send slice " << MySlice << " (" << MyRank << ") to " << SendSlice << " (" << SendRank
      << "), receive slice from " << RecvSlice << " (" << RecvRank << ")" << std::endl;
      gridHighDim->SendToRecvFromBegin(reqs,MyData,SendRank,RecvData,RecvRank,SliceSize);
      //memcpy(RecvData,MyData,SliceSize); // Debug
    }
    gridHighDim->SendToRecvFromComplete(reqs);
    std::cout << GridLogMessage << "Slice data shared." << std::endl;
    CCC_DEBUG_DUMP(Buffer, NumSlices, SliceSize);
#endif
  }
}

#define END_GRID_NAMESPACE }

*/

END_HADRONS_NAMESPACE

#endif // Distil_Vectors_hpp_
