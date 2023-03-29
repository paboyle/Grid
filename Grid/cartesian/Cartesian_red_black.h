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

NAMESPACE_BEGIN(Grid);

static const int CbRed  =0;
static const int CbBlack=1;
static const int Even   =CbRed;
static const int Odd    =CbBlack;

accelerator_inline int RedBlackCheckerBoardFromOindex (int oindex,const Coordinate &rdim,const Coordinate &chk_dim_msk)
{
  int nd=rdim.size();
  Coordinate coor(nd);

  Lexicographic::CoorFromIndex(coor,oindex,rdim);

  int linear=0;
  for(int d=0;d<nd;d++){
    if(chk_dim_msk[d])
      linear=linear+coor[d];
  }
  return (linear&0x1);
}

    
// Specialise this for red black grids storing half the data like a chess board.
class GridRedBlackCartesian : public GridBase
{
public:
  //  Coordinate _checker_dim_mask;
  int              _checker_dim;
  std::vector<int> _checker_board;

  virtual int CheckerBoarded(int dim){
    if( dim==_checker_dim) return 1;
    else return 0;
  }
  virtual int CheckerBoard(const Coordinate &site){
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
  virtual int  CheckerBoardFromOindexTable (int Oindex) {
    return _checker_board[Oindex];
  }
  virtual int  CheckerBoardFromOindex (int Oindex)
  {
    Coordinate ocoor;
    oCoorFromOindex(ocoor,Oindex);
    return CheckerBoard(ocoor);
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

  ////////////////////////////////////////////////////////////
  // Create Redblack from original grid; require full grid pointer ?
  ////////////////////////////////////////////////////////////
  GridRedBlackCartesian(const GridBase *base) : GridBase(base->_processors,*base)
  {
    int dims = base->_ndimension;
    Coordinate checker_dim_mask(dims,1);
    int checker_dim = 0;
    Init(base->_fdimensions,base->_simd_layout,base->_processors,checker_dim_mask,checker_dim);
  };

  ////////////////////////////////////////////////////////////
  // Create redblack from original grid, with non-trivial checker dim mask
  ////////////////////////////////////////////////////////////
  GridRedBlackCartesian(const GridBase *base,
			const Coordinate &checker_dim_mask,
			int checker_dim
			) :  GridBase(base->_processors,*base) 
  {
    Init(base->_fdimensions,base->_simd_layout,base->_processors,checker_dim_mask,checker_dim)  ;
  }

  virtual ~GridRedBlackCartesian() = default;

  void Init(const Coordinate &dimensions,
	    const Coordinate &simd_layout,
	    const Coordinate &processor_grid,
	    const Coordinate &checker_dim_mask,
	    int checker_dim)
  {

      _isCheckerBoarded = true;
    _checker_dim = checker_dim;
    assert(checker_dim_mask[checker_dim] == 1);
    _ndimension = dimensions.size();
    assert(checker_dim_mask.size() == _ndimension);
    assert(processor_grid.size() == _ndimension);
    assert(simd_layout.size() == _ndimension);

    _fdimensions.resize(_ndimension);
    _gdimensions.resize(_ndimension);
    _ldimensions.resize(_ndimension);
    _rdimensions.resize(_ndimension);
    _simd_layout.resize(_ndimension);
    _lstart.resize(_ndimension);
    _lend.resize(_ndimension);

    _ostride.resize(_ndimension);
    _istride.resize(_ndimension);

    _fsites = _gsites = _osites = _isites = 1;

    _checker_dim_mask = checker_dim_mask;

    for (int d = 0; d < _ndimension; d++)
      {
        _fdimensions[d] = dimensions[d];
        _gdimensions[d] = _fdimensions[d];
        _fsites = _fsites * _fdimensions[d];
        _gsites = _gsites * _gdimensions[d];

        if (d == _checker_dim)
	  {
	    assert((_gdimensions[d] & 0x1) == 0);
	    _gdimensions[d] = _gdimensions[d] / 2; // Remove a checkerboard
	    _gsites /= 2;
	  }
        _ldimensions[d] = _gdimensions[d] / _processors[d];
        assert(_ldimensions[d] * _processors[d] == _gdimensions[d]);
        _lstart[d] = _processor_coor[d] * _ldimensions[d];
        _lend[d] = _processor_coor[d] * _ldimensions[d] + _ldimensions[d] - 1;

        // Use a reduced simd grid
        _simd_layout[d] = simd_layout[d];
        _rdimensions[d] = _ldimensions[d] / _simd_layout[d]; // this is not checking if this is integer
        assert(_rdimensions[d] * _simd_layout[d] == _ldimensions[d]);
        assert(_rdimensions[d] > 0);

        // all elements of a simd vector must have same checkerboard.
        // If Ls vectorised, this must still be the case; e.g. dwf rb5d
        if (_simd_layout[d] > 1)
	  {
	    if (checker_dim_mask[d])
	      {
		assert((_rdimensions[d] & 0x1) == 0);
	      }
	  }

        _osites *= _rdimensions[d];
        _isites *= _simd_layout[d];

        // Addressing support
        if (d == 0)
	  {
	    _ostride[d] = 1;
	    _istride[d] = 1;
	  }
        else
	  {
	    _ostride[d] = _ostride[d - 1] * _rdimensions[d - 1];
	    _istride[d] = _istride[d - 1] * _simd_layout[d - 1];
	  }
      }

    ////////////////////////////////////////////////////////////////////////////////////////////
    // subplane information
    ////////////////////////////////////////////////////////////////////////////////////////////
    _slice_block.resize(_ndimension);
    _slice_stride.resize(_ndimension);
    _slice_nblock.resize(_ndimension);

    int block = 1;
    int nblock = 1;
    for (int d = 0; d < _ndimension; d++)
      nblock *= _rdimensions[d];

    for (int d = 0; d < _ndimension; d++)
      {
        nblock /= _rdimensions[d];
        _slice_block[d] = block;
        _slice_stride[d] = _ostride[d] * _rdimensions[d];
        _slice_nblock[d] = nblock;
        block = block * _rdimensions[d];
      }

    ////////////////////////////////////////////////
    // Create a checkerboard lookup table
    ////////////////////////////////////////////////
    int rvol = 1;
    for (int d = 0; d < _ndimension; d++)
      {
        rvol = rvol * _rdimensions[d];
      }
    _checker_board.resize(rvol);
    for (int osite = 0; osite < _osites; osite++)
      {
        _checker_board[osite] = CheckerBoardFromOindex(osite);
      }
  };

protected:
  virtual int oIndex(Coordinate &coor)
  {
    int idx = 0;
    for (int d = 0; d < _ndimension; d++)
      {
        if (d == _checker_dim)
	  {
	    idx += _ostride[d] * ((coor[d] / 2) % _rdimensions[d]);
	  }
        else
	  {
	    idx += _ostride[d] * (coor[d] % _rdimensions[d]);
	  }
      }
    return idx;
  };

  virtual int iIndex(Coordinate &lcoor)
  {
    int idx = 0;
    for (int d = 0; d < _ndimension; d++)
      {
        if (d == _checker_dim)
	  {
	    idx += _istride[d] * (lcoor[d] / (2 * _rdimensions[d]));
	  }
        else
	  {
	    idx += _istride[d] * (lcoor[d] / _rdimensions[d]);
	  }
      }
    return idx;
  }
};
NAMESPACE_END(Grid);
#endif
