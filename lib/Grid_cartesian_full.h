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
    virtual int CheckerBoardDestination(int cb,int shift){
        return 0;
    }
    virtual int CheckerBoardShift(int source_cb,int dim,int shift, int osite){
        return shift;
    }
    GridCartesian(std::vector<int> &dimensions,
		  std::vector<int> &simd_layout,
		  std::vector<int> &processor_grid
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
            
        _osites = 1;
        _isites = 1;
        for(int d=0;d<_ndimension;d++){
	  _fdimensions[d] = dimensions[d]; // Global dimensions
	  _gdimensions[d] = _fdimensions[d]; // Global dimensions
	  _simd_layout[d] = simd_layout[d];

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
