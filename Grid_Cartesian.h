#ifndef GRID_CARTESIAN_H
#define GRID_CARTESIAN_H

#include "Grid.h"

namespace dpo{

    
/////////////////////////////////////////////////////////////////////////////////////////
// Grid Support. Following will go into Grid.h.
/////////////////////////////////////////////////////////////////////////////////////////
    // Cartesian grids
    // dpo::Grid
    // dpo::GridCartesian
    // dpo::GridCartesianRedBlack
    
class Grid {
public:
    // Give Lattice access
    template<class object> friend class Lattice;

        
//protected:
        
    // Lattice wide random support. not yet fully implemented. Need seed strategy
    // and one generator per site.
    //std::default_random_engine generator;
  //    static std::mt19937  generator( 9 );

        
    // Grid information.
    unsigned long _ndimension;
    std::vector<int> _layout;     // Which dimensions get relayed out over simd lanes.
    std::vector<int> _dimensions; // Dimensions of array
    std::vector<int> _rdimensions;// Reduced dimensions with simd lane images removed
    std::vector<int> _ostride;    // Outer stride for each dimension
    std::vector<int> _istride;    // Inner stride i.e. within simd lane
    int _osites;                  // _isites*_osites = product(dimensions).
    int _isites;
        
    // subslice information
    std::vector<int> _slice_block;
    std::vector<int> _slice_stride;
    std::vector<int> _slice_nblock;
public:
    
    // These routines are key. Subdivide the linearised cartesian index into
    //      "inner" index identifying which simd lane of object<vFcomplex> is associated with coord
    //      "outer" index identifying which element of _odata in class "Lattice" is associated with coord.
    // Compared to, say, Blitz++ we simply need to store BOTH an inner stride and an outer
    // stride per dimension. The cost of evaluating the indexing information is doubled for an n-dimensional
    // coordinate. Note, however, for data parallel operations the "inner" indexing cost is not paid and all
    // lanes are operated upon simultaneously.
    
    inline int oIndexReduced(std::vector<int> &rcoor)
    {
        int idx=0;
        for(int d=0;d<_ndimension;d++) idx+=_ostride[d]*rcoor[d];
        return idx;
    }
    virtual int oIndex(std::vector<int> &coor)
    {
        int idx=0;
        for(int d=0;d<_ndimension;d++) idx+=_ostride[d]*(coor[d]%_rdimensions[d]);
        return idx;
    }
    inline int iIndex(std::vector<int> &rcoor)
    {
        int idx=0;
        for(int d=0;d<_ndimension;d++) idx+=_istride[d]*(rcoor[d]/_rdimensions[d]);
        return idx;
    }
        
    inline int oSites(void) { return _osites; };
    inline int iSites(void) { return _isites; };
    virtual int CheckerBoard(std::vector<int> site)=0;
    virtual int CheckerBoardDestination(int source_cb,int shift)=0;
    virtual int CheckerBoardShift(int source_cb,int dim,int shift)=0;
};

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
#endif
