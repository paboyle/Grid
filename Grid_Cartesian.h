#ifndef GRID_CARTESIAN_H
#define GRID_CARTESIAN_H

#include <Grid.h>
#include <Grid_Communicator.h>
namespace Grid{
    
/////////////////////////////////////////////////////////////////////////////////////////
// Grid Support. Following will go into Grid.h.
/////////////////////////////////////////////////////////////////////////////////////////
    // Cartesian grids
    // Grid::Grid
    // Grid::GridCartesian
    // Grid::GridCartesianRedBlack


class SimdGrid : public CartesianCommunicator {
public:

  
 SimdGrid(std::vector<int> & processor_grid) : CartesianCommunicator(processor_grid) {};

    // Give Lattice access
    template<class object> friend class Lattice;
        
//protected:
        
  // Lattice wide random support. not yet fully implemented. Need seed strategy
  // and one generator per site.
  //std::default_random_engine generator;
  //    static std::mt19937  generator( 9 );
        
    // Grid information.

    // Commicator provides
    //    unsigned long _ndimension;
    //    std::vector<int> _processors; // processor grid
    //    int              _processor;  // linear processor rank
    //    std::vector<int> _processor_coor;  // linear processor rank

    std::vector<int> _simd_layout;     // Which dimensions get relayed out over simd lanes.


    std::vector<int> _fdimensions;// Global dimensions of array prior to cb removal
    std::vector<int> _gdimensions;// Global dimensions of array after cb removal
    std::vector<int> _ldimensions;// local dimensions of array with processor images removed
    std::vector<int> _rdimensions;// Reduced local dimensions with simd lane images and processor images removed 

    //    std::vector<int> _lstart;     // local start of array in gcoors. _processor_coor[d]*_ldimensions[d]
    //    std::vector<int> _lend;       // local end of array in gcoors    _processor_coor[d]*_ldimensions[d]+_ldimensions_[d]-1

 
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
    inline int iCoordFromIsite(int lane,int mu)
    {
      std::vector<int> coor(_ndimension);
      for(int d=0;d<_ndimension;d++){
	coor[d] = lane % _simd_layout[d];
	lane    = lane / _simd_layout[d];
      }
      return coor[mu];
    }
    
    inline int oSites(void) { return _osites; };
    inline int iSites(void) { return _isites; };

    inline int CheckerBoardFromOsite (int Osite){
      std::vector<int> ocoor;
      CoordFromOsite(ocoor,Osite);
      int ss=0;
      for(int d=0;d<_ndimension;d++){
	ss=ss+ocoor[d];
      }      
      return ss&0x1;
    }
    inline void CoordFromOsite (std::vector<int>& coor,int Osite){
      coor.resize(_ndimension);
      for(int d=0;d<_ndimension;d++){
	coor[d] = Osite % _rdimensions[d];
	Osite   = Osite / _rdimensions[d];
      }
    }

    virtual int CheckerBoarded(int dim)=0;
    virtual int CheckerBoard(std::vector<int> site)=0;
    virtual int CheckerBoardDestination(int source_cb,int shift)=0;
    virtual int CheckerBoardShift(int source_cb,int dim,int shift,int osite)=0;
};

class GridCartesian: public SimdGrid {
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
		  ) : SimdGrid(processor_grid)
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
            
        if ( _isites != vComplex::Nsimd()) {
            printf("bad layout for grid isites %d Nsimd %d\n",_isites,vComplex::Nsimd());
            exit(0);
        }
    };
};
 
// Specialise this for red black grids storing half the data like a chess board.
class GridRedBlackCartesian : public SimdGrid
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
      int ocb=CheckerBoardFromOsite(osite);
	  
      if ( (source_cb+ocb)&1 ) {
	printf("Checkerboard shift %d\n",(shift)/2);
	return (shift)/2;
      } else {
	printf("Checkerboard shift %d\n",(shift+1)/2);
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
			  std::vector<int> &processor_grid) : SimdGrid(processor_grid)
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
