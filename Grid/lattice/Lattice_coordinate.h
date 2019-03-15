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
#ifndef GRID_LATTICE_COORDINATE_H
#define GRID_LATTICE_COORDINATE_H

namespace Grid {

    template<class iobj> inline void LatticeCoordinate(Lattice<iobj> &l,int mu)
    {
      typedef typename iobj::scalar_type scalar_type;
      typedef typename iobj::vector_type vector_type;

      GridBase *grid = l._grid;
      int Nsimd = grid->iSites();

      std::vector<int> gcoor;
      std::vector<scalar_type> mergebuf(Nsimd);

      vector_type vI;
      for(int o=0;o<grid->oSites();o++){
	for(int i=0;i<grid->iSites();i++){
	  grid->RankIndexToGlobalCoor(grid->ThisRank(),o,i,gcoor);
	  mergebuf[i]=(Integer)gcoor[mu];
	}
	merge<vector_type,scalar_type>(vI,mergebuf);
	l._odata[o]=vI;
      }
    };

}
#endif
