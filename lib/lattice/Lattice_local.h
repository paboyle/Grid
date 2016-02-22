    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/lattice/Lattice_local.h

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
#ifndef GRID_LATTICE_LOCALREDUCTION_H
#define GRID_LATTICE_LOCALREDUCTION_H

///////////////////////////////////////////////
// localInner, localNorm, outerProduct
///////////////////////////////////////////////

namespace Grid {

    /////////////////////////////////////////////////////
    // Non site, reduced locally reduced routines
    /////////////////////////////////////////////////////

    // localNorm2,
    template<class vobj>
    inline auto localNorm2 (const Lattice<vobj> &rhs)-> Lattice<typename vobj::tensor_reduced>
    {
      Lattice<typename vobj::tensor_reduced> ret(rhs._grid);
PARALLEL_FOR_LOOP
        for(int ss=0;ss<rhs._grid->oSites(); ss++){
	  ret._odata[ss]=innerProduct(rhs._odata[ss],rhs._odata[ss]);
        }
        return ret;
    }
    
    // localInnerProduct
    template<class vobj>
    inline auto localInnerProduct (const Lattice<vobj> &lhs,const Lattice<vobj> &rhs) -> Lattice<typename vobj::tensor_reduced>
    {
      Lattice<typename vobj::tensor_reduced> ret(rhs._grid);
PARALLEL_FOR_LOOP
      for(int ss=0;ss<rhs._grid->oSites(); ss++){
	ret._odata[ss]=innerProduct(lhs._odata[ss],rhs._odata[ss]);
      }
      return ret;
    }
    
    // outerProduct Scalar x Scalar -> Scalar
    //              Vector x Vector -> Matrix
    template<class ll,class rr>
    inline auto outerProduct (const Lattice<ll> &lhs,const Lattice<rr> &rhs) -> Lattice<decltype(outerProduct(lhs._odata[0],rhs._odata[0]))>
    {
        Lattice<decltype(outerProduct(lhs._odata[0],rhs._odata[0]))> ret(rhs._grid);
PARALLEL_FOR_LOOP
        for(int ss=0;ss<rhs._grid->oSites(); ss++){
            ret._odata[ss]=outerProduct(lhs._odata[ss],rhs._odata[ss]);
        }
        return ret;
     }

}

#endif
