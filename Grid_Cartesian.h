#include "Grid.h"
#include "Grid_vComplexD.h"

namespace dpo{
class GridCartesian: public Grid {
public:
    virtual int CheckerBoard(std::vector<int> site){
        return 0;
    }
    virtual int CheckerBoardDestination(int cb,int shift){
        return 0;
    }
    virtual int CheckerBoardShift(int source_cb,int dim,int shift){
        return shift;
    }
    GridCartesian(std::vector<int> &dimensions,std::vector<int> layout)
    {
        ///////////////////////
        // Grid information
        ///////////////////////
        _ndimension = dimensions.size();
            
        _dimensions.resize(_ndimension);
        _rdimensions.resize(_ndimension);
        _layout.resize(_ndimension);
            
        _ostride.resize(_ndimension);
        _istride.resize(_ndimension);
            
        _osites = 1;
        _isites = 1;
        for(int d=0;d<_ndimension;d++){
            _dimensions[d] = dimensions[d];
            _layout[d]     = layout[d];
                
            // Use a reduced simd grid
            _rdimensions[d]= _dimensions[d]/_layout[d];
            _osites *= _rdimensions[d];
            _isites *= _layout[d];
                
            // Addressing support
            if ( d==0 ) {
                _ostride[d] = 1;
                _istride[d] = 1;
            } else {
                _ostride[d] = _ostride[d-1]*_rdimensions[d-1];
                _istride[d] = _istride[d-1]*_layout[d-1];
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
            
        if ( _isites != vComplex::Nsimd()) {
            printf("bad layout for grid isites %d Nsimd %d\n",_isites,vComplex::Nsimd());
            exit(0);
        }
    };
};

// Specialise this for red black grids storing half the data like a chess board.
class GridRedBlackCartesian : public Grid
{
public:
    virtual int CheckerBoard(std::vector<int> site){
        return site[0]&0x1;
    }
    virtual int CheckerBoardShift(int source_cb,int dim,int shift){
        if ( dim == 0 ){
            int fulldim =2*_dimensions[0];
            shift = (shift+fulldim)%fulldim;
            if ( source_cb ) {
                // Shift 0,1 -> 0
                return (shift+1)/2;
            } else {
                // Shift 0->0, 1 -> 1, 2->1
                return (shift)/2;
            }
        } else {
            return shift;
        }
    }
    virtual int CheckerBoardDestination(int source_cb,int shift){
        if ((shift+2*_dimensions[0])&0x1) {
            return 1-source_cb;
        } else {
            return source_cb;
        }
    };
    GridRedBlackCartesian(std::vector<int> &dimensions,std::vector<int> layout)
    {
    ///////////////////////
    // Grid information
    ///////////////////////
        _ndimension = dimensions.size();
        
        _dimensions.resize(_ndimension);
        _rdimensions.resize(_ndimension);
        _layout.resize(_ndimension);
        
        _ostride.resize(_ndimension);
        _istride.resize(_ndimension);
        
        _osites = 1;
        _isites = 1;
        for(int d=0;d<_ndimension;d++){
            _dimensions[d] = dimensions[d];
            _dimensions[0] = dimensions[0]/2;
            _layout[d]     = layout[d];
                
            // Use a reduced simd grid
            _rdimensions[d]= _dimensions[d]/_layout[d];
            _osites *= _rdimensions[d];
            _isites *= _layout[d];
                
            // Addressing support
            if ( d==0 ) {
                _ostride[d] = 1;
                _istride[d] = 1;
            } else {
                _ostride[d] = _ostride[d-1]*_rdimensions[d-1];
                _istride[d] = _istride[d-1]*_layout[d-1];
            }
        }
            
        ////////////////////////////////////////////////////////////////////////////////////////////
        // subplane information
        // It may be worth the investment of generating a more general subplane "iterator",
        // and providing support for threads grabbing a unit of allocation.
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
            
        if ( _isites != vComplex::Nsimd()) {
            printf("bad layout for grid isites %d Nsimd %d\n",_isites,vComplex::Nsimd());
            exit(0);
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
