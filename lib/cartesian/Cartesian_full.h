    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/cartesian/Cartesian_full.h

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
#ifndef GRID_CARTESIAN_FULL_H
#define GRID_CARTESIAN_FULL_H

namespace Grid{
    
/////////////////////////////////////////////////////////////////////////////////////////
// Grid Support.
/////////////////////////////////////////////////////////////////////////////////////////


class GridCartesian: public GridBase {

public:

    virtual int CheckerBoarded(int dim){
      return 0;
    }
    virtual int CheckerBoard(std::vector<int> site){
        return 0;
    }
    virtual int CheckerBoardDestination(int cb,int shift,int dim){
        return 0;
    }
    virtual int CheckerBoardShiftForCB(int source_cb,int dim,int shift, int ocb){
      return shift;
    }
    virtual int CheckerBoardShift(int source_cb,int dim,int shift, int osite){
      return shift;
    }
    GridCartesian(const std::vector<int> &dimensions,
		  const std::vector<int> &simd_layout,
		  const std::vector<int> &processor_grid
		  ) : GridBase(processor_grid)
    {
        ///////////////////////
        // Grid information
        ///////////////////////
        _ndimension = dimensions.size();
            
        _fdimensions.resize(_ndimension);
        _gdimensions.resize(_ndimension);
        _ldimensions.resize(_ndimension);
        _rdimensions.resize(_ndimension);
        _simd_layout.resize(_ndimension);
            
        _ostride.resize(_ndimension);
        _istride.resize(_ndimension);
            
        _fsites = _gsites = _osites = _isites = 1;

        for(int d=0;d<_ndimension;d++){
	  _fdimensions[d] = dimensions[d]; // Global dimensions
	  _gdimensions[d] = _fdimensions[d]; // Global dimensions
	  _simd_layout[d] = simd_layout[d];
	  _fsites = _fsites * _fdimensions[d];
	  _gsites = _gsites * _gdimensions[d];

	  //FIXME check for exact division

	  // Use a reduced simd grid
	  _ldimensions[d]= _gdimensions[d]/_processors[d];  //local dimensions
	  _rdimensions[d]= _ldimensions[d]/_simd_layout[d]; //overdecomposition
	  _osites *= _rdimensions[d];
	  _isites *= _simd_layout[d];
                
	  // Addressing support
	  if ( d==0 ) {
	    _ostride[d] = 1;
	    _istride[d] = 1;
	  } else {
	    _ostride[d] = _ostride[d-1]*_rdimensions[d-1];
	    _istride[d] = _istride[d-1]*_simd_layout[d-1];
	  }
        }
        
        ///////////////////////
        // subplane information
        ///////////////////////
        _slice_block.resize(_ndimension);
        _slice_stride.resize(_ndimension);
        _slice_nblock.resize(_ndimension);
            
        int block =1;
        int nblock=1;
        for(int d=0;d<_ndimension;d++) nblock*=_rdimensions[d];
            
        for(int d=0;d<_ndimension;d++){
            nblock/=_rdimensions[d];
            _slice_block[d] =block;
            _slice_stride[d]=_ostride[d]*_rdimensions[d];
            _slice_nblock[d]=nblock;
            block = block*_rdimensions[d];
        }

    };
};


}
#endif
