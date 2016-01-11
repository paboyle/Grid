    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/cartesian/Cartesian_red_black.h

    Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
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
#ifndef GRID_CARTESIAN_RED_BLACK_H
#define GRID_CARTESIAN_RED_BLACK_H


namespace Grid {

    static const int CbRed  =0;
    static const int CbBlack=1;
    static const int Even   =CbRed;
    static const int Odd    =CbBlack;

    // Perhaps these are misplaced and 
    // should be in sparse matrix.
    // Also should make these a named enum type
    static const int DaggerNo=0;
    static const int DaggerYes=1;

// Specialise this for red black grids storing half the data like a chess board.
class GridRedBlackCartesian : public GridBase
{
public:
    std::vector<int> _checker_dim_mask;
    int              _checker_dim;

    virtual int CheckerBoarded(int dim){
      if( dim==_checker_dim) return 1;
      else return 0;
    }
    virtual int CheckerBoard(std::vector<int> site){
      int linear=0;
      assert(site.size()==_ndimension);
      for(int d=0;d<_ndimension;d++){ 
	if(_checker_dim_mask[d])
	  linear=linear+site[d];
      }
      return (linear&0x1);
    }


    // Depending on the cb of site, we toggle source cb.
    // for block #b, element #e = (b, e)
    // we need 
    virtual int CheckerBoardShiftForCB(int source_cb,int dim,int shift,int ocb){
      if(dim != _checker_dim) return shift;

      int fulldim =_fdimensions[dim];
      shift = (shift+fulldim)%fulldim;

      // Probably faster with table lookup;
      // or by looping over x,y,z and multiply rather than computing checkerboard.
	  
      if ( (source_cb+ocb)&1 ) {

	return (shift)/2;
      } else {
	return (shift+1)/2;
      }
    }
    virtual int CheckerBoardShift(int source_cb,int dim,int shift,int osite){

      if(dim != _checker_dim) return shift;

      int ocb=CheckerBoardFromOindex(osite);
      
      return CheckerBoardShiftForCB(source_cb,dim,shift,ocb);
    }
    
    virtual int CheckerBoardDestination(int source_cb,int shift,int dim){
      if ( _checker_dim_mask[dim]  ) {
	// If _fdimensions[checker_dim] is odd, then shifting by 1 in other dims
	// does NOT cause a parity hop.
	int add=(dim==_checker_dim) ? 0 : _fdimensions[_checker_dim];
        if ( (shift+add) &0x1) {
            return 1-source_cb;
        } else {
            return source_cb;
        }
      } else {
	return source_cb;

      }
    };

    GridRedBlackCartesian(const GridBase *base) : GridRedBlackCartesian(base->_fdimensions,base->_simd_layout,base->_processors)  {};

    GridRedBlackCartesian(const std::vector<int> &dimensions,
			  const std::vector<int> &simd_layout,
			  const std::vector<int> &processor_grid,
			  const std::vector<int> &checker_dim_mask,
			  int checker_dim
			  ) :  GridBase(processor_grid) 
    {
      Init(dimensions,simd_layout,processor_grid,checker_dim_mask,checker_dim);
    }
    GridRedBlackCartesian(const std::vector<int> &dimensions,
			  const std::vector<int> &simd_layout,
			  const std::vector<int> &processor_grid) : GridBase(processor_grid) 
    {
      std::vector<int> checker_dim_mask(dimensions.size(),1);
      Init(dimensions,simd_layout,processor_grid,checker_dim_mask,0);
    }
    void Init(const std::vector<int> &dimensions,
	      const std::vector<int> &simd_layout,
	      const std::vector<int> &processor_grid,
	      const std::vector<int> &checker_dim_mask,
	      int checker_dim)
    {
    ///////////////////////
    // Grid information
    ///////////////////////
      _checker_dim = checker_dim;
      assert(checker_dim_mask[checker_dim]==1);
      _ndimension = dimensions.size();
      assert(checker_dim_mask.size()==_ndimension);
      assert(processor_grid.size()==_ndimension);
      assert(simd_layout.size()==_ndimension);
      
      _fdimensions.resize(_ndimension);
      _gdimensions.resize(_ndimension);
      _ldimensions.resize(_ndimension);
      _rdimensions.resize(_ndimension);
      _simd_layout.resize(_ndimension);
      
      _ostride.resize(_ndimension);
      _istride.resize(_ndimension);
      
      _fsites = _gsites = _osites = _isites = 1;
	
      _checker_dim_mask=checker_dim_mask;

      for(int d=0;d<_ndimension;d++){
	_fdimensions[d] = dimensions[d];
	_gdimensions[d] = _fdimensions[d];
	_fsites = _fsites * _fdimensions[d];
	_gsites = _gsites * _gdimensions[d];
        
	if (d==_checker_dim) {
	  _gdimensions[d] = _gdimensions[d]/2; // Remove a checkerboard
	}
	_ldimensions[d] = _gdimensions[d]/_processors[d];

	// Use a reduced simd grid
	_simd_layout[d] = simd_layout[d];
	_rdimensions[d]= _ldimensions[d]/_simd_layout[d];

	// all elements of a simd vector must have same checkerboard.
	if ( simd_layout[d]>1 ) assert((_rdimensions[d]&0x1)==0); 

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
            
      ////////////////////////////////////////////////////////////////////////////////////////////
      // subplane information
      ////////////////////////////////////////////////////////////////////////////////////////////
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
protected:
    virtual int oIndex(std::vector<int> &coor)
    {
      int idx=0;
      for(int d=0;d<_ndimension;d++) {
	if( d==_checker_dim ) {
	  idx+=_ostride[d]*((coor[d]/2)%_rdimensions[d]);
	} else {
	  idx+=_ostride[d]*(coor[d]%_rdimensions[d]);
	}
      }
        return idx;
    };
        
};

}
#endif
