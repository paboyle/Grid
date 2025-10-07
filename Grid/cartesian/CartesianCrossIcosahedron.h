/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/cartesian/CartesianCrossIcosahedron.h

    Copyright (C) 2025

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
    
/////////////////////////////////////////////////////////////////////////////////////////
// Grid Support.
/////////////////////////////////////////////////////////////////////////////////////////

enum IcosahedralMeshType {
  IcosahedralVertices,
  IcosahedralEdges
} ;
enum NorthSouth {
  North = 1,
  South = 0
};

const int num_icosahedron_tiles = 10;

class GridCartesianCrossIcosahedron: public GridCartesian {

public:

  IcosahedralMeshType meshType;

  IcosahedralMeshType MeshType(void) { return meshType; };
  
  /////////////////////////////////////////////////////////////////////////
  // Constructor takes a parent grid and possibly subdivides communicator.
  /////////////////////////////////////////////////////////////////////////
  /*
  GridCartesian(const Coordinate &dimensions,
		const Coordinate &simd_layout,
		const Coordinate &processor_grid,
		const GridCartesian &parent) : GridBase(processor_grid,parent,dummy)
  {
    assert(0); // No subdivision
  }
  GridCartesian(const Coordinate &dimensions,
		const Coordinate &simd_layout,
		const Coordinate &processor_grid,
		const GridCartesian &parent,int &split_rank) : GridBase(processor_grid,parent,split_rank)
  {
    assert(0); // No subdivision
  }
  */
  /////////////////////////////////////////////////////////////////////////
  // Construct from comm world
  /////////////////////////////////////////////////////////////////////////
  GridCartesianCrossIcosahedron(const Coordinate &dimensions,
				const Coordinate &simd_layout,
				const Coordinate &processor_grid,
				IcosahedralMeshType _meshType) : GridCartesian(dimensions,simd_layout,processor_grid)
  {
    meshType = _meshType;
    Coordinate S2dimensions=dimensions;
    Coordinate S2simd      =simd_layout;
    Coordinate S2procs     =processor_grid;

    assert(simd_layout[0]==1); // Force simd into perpendicular dimensions
    assert(simd_layout[1]==1); // to avoid pole storage complexity interacting with SIMD.
    assert(dimensions[_ndimension-1]==num_icosahedron_tiles);
    assert(processor_grid[_ndimension-1]<=2); // Keeps the patches that need a pole on the same node

    // allocate the pole storage if we are seeking vertex domain data
    if ( meshType == IcosahedralVertices ) {
      InitPoles();
    }
  }

  virtual ~GridCartesianCrossIcosahedron() = default;

  ////////////////////////////////////////////////
  // Use to decide if a given grid is icosahedral
  ////////////////////////////////////////////////
  int hasNorthPole;
  int hasSouthPole;
  int northPoleOsite;
  int southPoleOsite;
  int northPoleOsites;
  int southPoleOsites;

  virtual int Icosahedral(void)     override { return 1;}
  virtual int ownsNorthPole(void)   const override { return hasNorthPole; };
  virtual int NorthPoleOsite(void)  const override { return northPoleOsite; };
  virtual int NorthPoleOsites(void) const override { return northPoleOsites; };
  virtual int ownsSouthPole(void)   const override { return hasSouthPole; };
  virtual int SouthPoleOsite(void)  const override { return southPoleOsite; };
  virtual int SouthPoleOsites(void) const override { return southPoleOsites; };

  void InitPoles(void)
  {
    int Ndm1 = _ndimension-1;
    ///////////////////////
    // Add the extra pole storage
    ///////////////////////
    // Vertices = 1x LxLx D1...Dn + 2.D1...Dn
    // Start after the LxL and don't include the 10 patch dim
    int OrthogSize = 1;
    for (int d = 2; d < Ndm1; d++) {
      OrthogSize *= _gdimensions[d];
    }
    _fsites += OrthogSize*2;
    _gsites += OrthogSize*2;

    // Simd reduced sizes are multiplied up.
    // If the leading LxL are simd-ized, the vector objects will contain "redundant" lanes
    // which should contain identical north (south) pole data
    OrthogSize = 1;
    for (int d = 2; d < Ndm1; d++) {
      OrthogSize *= _rdimensions[d];
    }

    // Grow the local volume to hold pole data
    // on rank (0,0) in the LxL planes
    // since SIMD must be placed in the orthogonal directions
    Coordinate pcoor = this->ThisProcessorCoor();
    Coordinate pgrid = this->ProcessorGrid();

    const int xdim=0;
    const int ydim=1;
    /*
     *
     *  /\/\/\/\/\
     * /\/\/\/\/\/
     * \/\/\/\/\/
     *
     *  y
     * /
     * \x
     *
     * Labelling patches as 5 6 7 8 9
     *                      0 1 2 3 4
     *
     * Will ban distribution of the patch dimension by more than 2.
     *
     * Hence all 5 patches associated with the pole must have the
     * appropriate "corner" of the patch L^2 located on the SAME rank.
     */ 
     
    if( (pcoor[xdim]==pgrid[xdim]-1) && (pcoor[ydim]==0) && (pcoor[Ndm1]==0) ){
      hasSouthPole   =1;
      southPoleOsite=this->_osites;
      southPoleOsites=OrthogSize;
      this->_osites += OrthogSize;
    } else {
      hasSouthPole   =0;
      southPoleOsites=0;
      southPoleOsite=0;
    }
    if( (pcoor[xdim]==0) && (pcoor[ydim]==pgrid[ydim]-1) && (pcoor[Ndm1]==pgrid[Ndm1]-1) ){
      hasNorthPole   =1;
      northPoleOsite=this->_osites;
      northPoleOsites=OrthogSize;
      this->_osites += OrthogSize;
    } else {
      hasNorthPole   =0;
      northPoleOsites=0;
      northPoleOsite=0;
    }
    std::cout << "Icosahedral vertex field volume " << this->_osites<<std::endl;
    std::cout << "Icosahedral south pole offset   " << this->southPoleOsite<<std::endl;
    std::cout << "Icosahedral north pole offset   " << this->northPoleOsite<<std::endl;
    std::cout << "Icosahedral south pole size     " << this->southPoleOsites<<std::endl;
    std::cout << "Icosahedral north pole size     " << this->northPoleOsites<<std::endl;
  };

};

NAMESPACE_END(Grid);
