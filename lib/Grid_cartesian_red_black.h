#ifndef GRID_CARTESIAN_RED_BLACK_H
#define GRID_CARTESIAN_RED_BLACK_H


namespace Grid {

// Specialise this for red black grids storing half the data like a chess board.
class GridRedBlackCartesian : public GridBase
{
public:
    virtual int CheckerBoarded(int dim){
      if( dim==0) return 1;
      else return 0;
    }
    virtual int CheckerBoard(std::vector<int> site){
      return (site[0]+site[1]+site[2]+site[3])&0x1;
    }

    // Depending on the cb of site, we toggle source cb.
    // for block #b, element #e = (b, e)
    // we need 
    virtual int CheckerBoardShift(int source_cb,int dim,int shift,int osite){

      if(dim != 0) return shift;

      int fulldim =_fdimensions[0];
      shift = (shift+fulldim)%fulldim;

      // Probably faster with table lookup;
      // or by looping over x,y,z and multiply rather than computing checkerboard.
      int ocb=CheckerBoardFromOindex(osite);
	  
      if ( (source_cb+ocb)&1 ) {
	return (shift)/2;
      } else {
	return (shift+1)/2;
      }
    }

    virtual int CheckerBoardDestination(int source_cb,int shift){
        if ((shift+_fdimensions[0])&0x1) {
            return 1-source_cb;
        } else {
            return source_cb;
        }
    };
    GridRedBlackCartesian(std::vector<int> &dimensions,
			  std::vector<int> &simd_layout,
			  std::vector<int> &processor_grid) : GridBase(processor_grid)
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
            _fdimensions[d] = dimensions[d];
            _gdimensions[d] = _fdimensions[d];
            if (d==0) _gdimensions[0] = _gdimensions[0]/2; // Remove a checkerboard
            _ldimensions[d] = _gdimensions[d]/_processors[d];
                
            // Use a reduced simd grid
            _simd_layout[d] = simd_layout[d];
            _rdimensions[d]= _ldimensions[d]/_simd_layout[d];

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
        int idx=_ostride[0]*((coor[0]/2)%_rdimensions[0]);
        for(int d=1;d<_ndimension;d++) idx+=_ostride[d]*(coor[d]%_rdimensions[d]);
        return idx;
    };
        
};

}
#endif
