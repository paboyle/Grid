/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./lib/algorithms/SplitOperator.h

    Copyright (C) 2026

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
#pragma once

NAMESPACE_BEGIN(Grid);

template<class Field> class LinearOperatorBase;
template<class Field> class CheckerBoardedSparseMatrixBase;

/////////////////////////////////////////////////////////////////////////////////////////////
// A copy of an operator living on grids whose communicator is split into independent
// partitions, together with those grids. Produced by SplitClone().
//
// Owns everything it points to. Destruction is the reverse of creation: linear operator,
// matrix, then grids, so no object outlives a grid it references. For 4d operators the
// fermion grids are the gauge grids and are deleted once.
/////////////////////////////////////////////////////////////////////////////////////////////
template<class Field>
class SplitOperator
{
public:
  GridCartesian         *GaugeGrid     = nullptr;
  GridRedBlackCartesian *GaugeRBGrid   = nullptr;
  GridCartesian         *FermionGrid   = nullptr;
  GridRedBlackCartesian *FermionRBGrid = nullptr;

  int Partition  = 0;  // Grid_split vector index held by this rank's partition
  int Partitions = 1;  // number of partitions

  CheckerBoardedSparseMatrixBase<Field> *Matrix    = nullptr;
  LinearOperatorBase<Field>             *Linop     = nullptr;
  GridBase                              *FieldGrid = nullptr;  // grid of Linop's fields

  SplitOperator(void) {};

  SplitOperator(const SplitOperator &) = delete;

  SplitOperator &operator=(const SplitOperator &) = delete;

  ~SplitOperator(void)
  {
    delete Linop;
    delete Matrix;
    if ( FermionRBGrid != GaugeRBGrid ) {
      delete FermionRBGrid;
    }
    if ( FermionGrid != GaugeGrid ) {
      delete FermionGrid;
    }
    delete GaugeRBGrid;
    delete GaugeGrid;
  }
};

/////////////////////////////////////////////////////////////////////////////////////////////
// Index of the vector that Grid_split(std::vector<Field> full, Field split) delivers to this
// rank's partition. Grid_split orders partitions lexicographically with the first
// dimension fastest; the split communicator's own rank (srank) uses the reversed MPI
// convention, so the two differ whenever more than one dimension is split.
/////////////////////////////////////////////////////////////////////////////////////////////
inline int GridSplitVectorIndex(GridBase *full,GridBase *split)
{
  int nd = full->_ndimension;
  GRID_ASSERT(split->_ndimension == nd);

  Coordinate scoor(nd);
  Coordinate ssize(nd);
  for(int d=0;d<nd;d++){
    scoor[d] = full->ThisProcessorCoor()[d] / split->ProcessorGrid()[d];
    ssize[d] = full->ProcessorGrid()[d]     / split->ProcessorGrid()[d];
  }
  int index;
  Lexicographic::IndexFromCoor(scoor,index,ssize);
  return index;
}

/////////////////////////////////////////////////////////////////////////////////////////////
// Partition MPI layout for a batched solve on grid, from a request as given on the command
// line (--batched-solver-split). Uses the trailing dimensions of grid's processor and shm
// layouts, so it works for 4d and 5d grids. Returns grid's own processor layout (one
// partition) when no split is requested. Asserts divisibility; warns when partitions
// straddle nodes or when node boundaries are not visible.
/////////////////////////////////////////////////////////////////////////////////////////////
inline Coordinate BatchedSolverSplitLayout(GridBase *grid,
                                           const Coordinate &request,
                                           bool node,
                                           int &partitions)
{
  int nd  = GridDefaultMpi().size();
  int pad = grid->_ndimension - nd;
  GRID_ASSERT(pad >= 0);
  GRID_ASSERT(grid->ShmGrid().size() == grid->_ndimension);

  Coordinate processors(nd);
  Coordinate shm(nd);
  for(int d=0;d<nd;d++){
    processors[d] = grid->ProcessorGrid()[pad+d];
    shm[d]        = grid->ShmGrid()[pad+d];
  }

  Coordinate split(nd);
  if ( node ) {
    split = shm;
  } else if ( request.size() == 0 ) {
    split = processors;
  } else {
    GRID_ASSERT(request.size() == nd);
    split = request;
  }

  partitions = 1;
  for(int d=0;d<nd;d++){
    GRID_ASSERT( (processors[d] % split[d]) == 0 );
    partitions *= processors[d] / split[d];
  }
  if ( partitions == 1 ) {
    return split;
  }

  std::cout << GridLogMessage << "BatchedSolverSplit: partition layout " << split
            << " of " << processors << " : " << partitions << " partitions" << std::endl;

  int inside_node  = 1;
  int whole_nodes  = 1;
  int shm_trivial  = 1;
  for(int d=0;d<nd;d++){
    if ( (shm[d] % split[d]) != 0 ) {
      inside_node = 0;
    }
    if ( (split[d] % shm[d]) != 0 ) {
      whole_nodes = 0;
    }
    if ( shm[d] != 1 ) {
      shm_trivial = 0;
    }
  }
  if ( shm_trivial ) {
    std::cout << GridLogWarning << "BatchedSolverSplit: shm layout is trivial (one rank per node,"
              << " or shared memory disabled); node locality of partitions not checked" << std::endl;
  } else if ( !inside_node && !whole_nodes ) {
    std::cout << GridLogWarning << "BatchedSolverSplit: partitions " << split
              << " straddle node boundaries (node layout " << shm << "); inner solves will communicate off node" << std::endl;
  }
  return split;
}

NAMESPACE_END(Grid);
