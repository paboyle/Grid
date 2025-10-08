/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/lattice/Lattice_coordinate.h

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
#pragma once 

NAMESPACE_BEGIN(Grid);

template<class iobj> inline void LatticeCoordinate(Lattice<iobj> &l,int mu)
{
  typedef typename iobj::scalar_type scalar_type;
  typedef typename iobj::vector_type vector_type;

  l=Zero();
  
  GridBase *grid = l.Grid();
  int Nsimd = grid->iSites();

  int cartesian_vol = grid->oSites();
  if ( grid->Icosahedral() ) {
    cartesian_vol = cartesian_vol - grid->NorthPoleOsites()-grid->SouthPoleOsites();
  }
  {
    autoView(l_v, l, CpuWrite);
    thread_for( o, cartesian_vol, {
	vector_type vI;
	Coordinate gcoor;
	ExtractBuffer<scalar_type> mergebuf(Nsimd);
	for(int i=0;i<grid->iSites();i++){
	  grid->RankIndexToGlobalCoor(grid->ThisRank(),o,i,gcoor);
	  mergebuf[i]=(Integer)gcoor[mu];
	}
	merge<vector_type,scalar_type>(vI,mergebuf);
	l_v[o]=vI;
      });
  }

  if (grid->IcosahedralVertices()) {
    uint64_t psites=1;
    Coordinate perpdims;
    typename iobj::scalar_object ss;
    for(int d=2;d<grid->_ndimension-1;d++){
      int pd=grid->_gdimensions[d];
      psites*=pd;
      perpdims.push_back(pd);
    }
    for(uint64_t p=0;p<psites;p++){
      Coordinate orthog;
      Lexicographic::CoorFromIndex(orthog,p,perpdims);

      int icoor;
      if ( mu>=2 && mu < grid->_ndimension-1) {
	icoor = orthog[mu-2];
      } else {
	icoor = -1;
      }

      ss=scalar_type(icoor);

      pokePole(ss,l,orthog,South);
      pokePole(ss,l,orthog,North);
    }
  }
};
template<class iobj> inline void LatticePole(Lattice<iobj> &l,NorthSouth pole)
{
  typedef typename iobj::scalar_object sobj;
  typedef typename iobj::scalar_type scalar_type;
  typedef typename iobj::vector_type vector_type;

  GridBase *grid = l.Grid();

  l=Zero();

  assert(grid->IcosahedralVertices());
  
  if (grid->IcosahedralVertices()) {
    uint64_t psites=1;
    Coordinate perpdims;
    sobj ss;
    scalar_type one(1.0);
    ss=one;
    for(int d=2;d<l.Grid()->_ndimension-1;d++){
      int pd=l.Grid()->_gdimensions[d];
      psites*=pd;
      perpdims.push_back(pd);
    }
    for(uint64_t p=0;p<psites;p++){
      Coordinate orthog;
      Lexicographic::CoorFromIndex(orthog,p,perpdims);
      pokePole(ss,l,orthog,pole);
    }
  }
};

NAMESPACE_END(Grid);
