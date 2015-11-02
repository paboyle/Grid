#include "Grid.h"

namespace Grid {


  CartesianStencil::CartesianStencil(GridBase *grid,
				     int npoints,
				     int checkerboard,
				     const std::vector<int> &directions,
				     const std::vector<int> &distances) 
    :   _entries(npoints), _permute_type(npoints), _comm_buf_size(npoints)
    {
      _npoints = npoints;
      _grid    = grid;
      _directions = directions;
      _distances  = distances;
      _unified_buffer_size=0;
      _request_count =0;

      int osites  = _grid->oSites();

      for(int i=0;i<npoints;i++){

	int point = i;

	_entries[i].resize( osites);

	int dimension    = directions[i];
	int displacement = distances[i];
	int shift = displacement;
	
	int fd = _grid->_fdimensions[dimension];
	int rd = _grid->_rdimensions[dimension];
	_permute_type[point]=_grid->PermuteType(dimension);

	_checkerboard = checkerboard;

	// the permute type
	int simd_layout     = _grid->_simd_layout[dimension];
	int comm_dim        = _grid->_processors[dimension] >1 ;
	int splice_dim      = _grid->_simd_layout[dimension]>1 && (comm_dim);

	int sshift[2];
	
	// Underlying approach. For each local site build
	// up a table containing the npoint "neighbours" and whether they 
	// live in lattice or a comms buffer.
	if ( !comm_dim ) {
	  sshift[0] = _grid->CheckerBoardShiftForCB(_checkerboard,dimension,shift,Even);
	  sshift[1] = _grid->CheckerBoardShiftForCB(_checkerboard,dimension,shift,Odd);

	  if ( sshift[0] == sshift[1] ) {
	    Local(point,dimension,shift,0x3);
	  } else {
	    Local(point,dimension,shift,0x1);// if checkerboard is unfavourable take two passes
	    Local(point,dimension,shift,0x2);// both with block stride loop iteration
	  }
	} else { // All permute extract done in comms phase prior to Stencil application
	  //        So tables are the same whether comm_dim or splice_dim
	  sshift[0] = _grid->CheckerBoardShiftForCB(_checkerboard,dimension,shift,Even);
	  sshift[1] = _grid->CheckerBoardShiftForCB(_checkerboard,dimension,shift,Odd);
	  if ( sshift[0] == sshift[1] ) {
	    Comms(point,dimension,shift,0x3);
	  } else {
	    Comms(point,dimension,shift,0x1);// if checkerboard is unfavourable take two passes
	    Comms(point,dimension,shift,0x2);// both with block stride loop iteration
	  }
	}
      }
    }


    void CartesianStencil::Local     (int point, int dimension,int shiftpm,int cbmask)
    {
      int fd = _grid->_fdimensions[dimension];
      int rd = _grid->_rdimensions[dimension];
      int ld = _grid->_ldimensions[dimension];
      int gd = _grid->_gdimensions[dimension];
      
      // Map to always positive shift modulo global full dimension.
      int shift = (shiftpm+fd)%fd;
      
      // the permute type
      int permute_dim =_grid->PermuteDim(dimension);
      
      for(int x=0;x<rd;x++){       
	
	int o   = 0;
	int bo  = x * _grid->_ostride[dimension];
	
	int cb= (cbmask==0x2)? Odd : Even;
	  
	int sshift = _grid->CheckerBoardShiftForCB(_checkerboard,dimension,shift,cb);
	int sx     = (x+sshift)%rd;

	int wraparound=0;
	if ( (shiftpm==-1) && (sx>x)  ) {
	  wraparound = 1;
	}
	if ( (shiftpm== 1) && (sx<x)  ) {
	  wraparound = 1;
	}

	  
	int permute_slice=0;
	if(permute_dim){
	  int wrap = sshift/rd;
	  int  num = sshift%rd;
	  if ( x< rd-num ) permute_slice=wrap;
	  else permute_slice = 1-wrap;
	}

  	CopyPlane(point,dimension,x,sx,cbmask,permute_slice,wraparound);
  
      }
    }

    void CartesianStencil::Comms     (int point,int dimension,int shiftpm,int cbmask)
    {
      GridBase *grid=_grid;
      
      int fd              = _grid->_fdimensions[dimension];
      int ld              = _grid->_ldimensions[dimension];
      int rd              = _grid->_rdimensions[dimension];
      int pd              = _grid->_processors[dimension];
      int simd_layout     = _grid->_simd_layout[dimension];
      int comm_dim        = _grid->_processors[dimension] >1 ;
      
      //      assert(simd_layout==1); // Why?
      assert(comm_dim==1);
      int shift = (shiftpm + fd) %fd;
      assert(shift>=0);
      assert(shift<fd);
      
      int buffer_size = _grid->_slice_nblock[dimension]*_grid->_slice_block[dimension];
      _comm_buf_size[point] = buffer_size; // Size of _one_ plane. Multiple planes may be gathered and
                                           // send to one or more remote nodes.

      int cb= (cbmask==0x2)? Odd : Even;
      int sshift= _grid->CheckerBoardShiftForCB(_checkerboard,dimension,shift,cb);
      
      for(int x=0;x<rd;x++){       
	
	int offnode = (((x+sshift)%fd) >= rd ); 
	//	int comm_proc   = ((x+sshift)/ld)%pd;        
	//	int offnode     = (comm_proc!=0);
	int sx          = (x+sshift)%rd;

	int wraparound=0;
	if ( (shiftpm==-1) && (sx>x) && (grid->_processor_coor[dimension]==0) ) {
	  wraparound = 1;
	}
	if ( (shiftpm== 1) && (sx<x) && (grid->_processor_coor[dimension]==grid->_processors[dimension]-1) ) {
	  wraparound = 1;
	}
	if (!offnode) {
	  
	  int permute_slice=0;
	  CopyPlane(point,dimension,x,sx,cbmask,permute_slice,wraparound); 
	  
	} else {
	  
	  int words = buffer_size;
	  if (cbmask != 0x3) words=words>>1;
	  
	  //	  GatherPlaneSimple (point,dimension,sx,cbmask);
	  
	  int rank           = grid->_processor;
	  int recv_from_rank;
	  int xmit_to_rank;

	  int unified_buffer_offset = _unified_buffer_size;
	  _unified_buffer_size    += words;
	  ScatterPlane(point,dimension,x,cbmask,unified_buffer_offset,wraparound); // permute/extract/merge is done in comms phase
	  
	}
      }
    }
  // Routine builds up integer table for each site in _offsets, _is_local, _permute
  void CartesianStencil::CopyPlane(int point, int dimension,int lplane,int rplane,int cbmask,int permute,int wrap)
    {
      int rd = _grid->_rdimensions[dimension];
      
      if ( !_grid->CheckerBoarded(dimension) ) {
	
	int o   = 0;                                     // relative offset to base within plane
	int ro  = rplane*_grid->_ostride[dimension]; // base offset for start of plane 
	int lo  = lplane*_grid->_ostride[dimension]; // offset in buffer
	
	// Simple block stride gather of SIMD objects
	for(int n=0;n<_grid->_slice_nblock[dimension];n++){
	  for(int b=0;b<_grid->_slice_block[dimension];b++){
	    _entries[point][lo+o+b]._offset  =ro+o+b;
	    _entries[point][lo+o+b]._is_local=1;
	    _entries[point][lo+o+b]._permute=permute;
	    _entries[point][lo+o+b]._around_the_world=wrap;
	  }
	  o +=_grid->_slice_stride[dimension];
	}
	
      } else {
	
	int ro  = rplane*_grid->_ostride[dimension]; // base offset for start of plane 
	int lo  = lplane*_grid->_ostride[dimension]; // base offset for start of plane 
	int o   = 0;                                     // relative offset to base within plane
	
	for(int n=0;n<_grid->_slice_nblock[dimension];n++){
	  for(int b=0;b<_grid->_slice_block[dimension];b++){
	    
	    int ocb=1<<_grid->CheckerBoardFromOindex(o+b);
	    
	    if ( ocb&cbmask ) {
	      _entries[point][lo+o+b]._offset =ro+o+b;
	      _entries[point][lo+o+b]._is_local=1;
	      _entries[point][lo+o+b]._permute=permute;
	      _entries[point][lo+o+b]._around_the_world=wrap;
	    }
	    
	    }
	  o +=_grid->_slice_stride[dimension];
	}
	
      }
    }
  // Routine builds up integer table for each site in _offsets, _is_local, _permute
   void CartesianStencil::ScatterPlane (int point,int dimension,int plane,int cbmask,int offset, int wrap)
    {
      int rd = _grid->_rdimensions[dimension];
      
      if ( !_grid->CheckerBoarded(dimension) ) {
	
	int so  = plane*_grid->_ostride[dimension]; // base offset for start of plane 
	int o   = 0;                                    // relative offset to base within plane
	int bo  = 0;                                    // offset in buffer

	// Simple block stride gather of SIMD objects
	for(int n=0;n<_grid->_slice_nblock[dimension];n++){
	  for(int b=0;b<_grid->_slice_block[dimension];b++){
	    _entries[point][so+o+b]._offset  =offset+(bo++);
	    _entries[point][so+o+b]._is_local=0;
	    _entries[point][so+o+b]._permute=0;
	    _entries[point][so+o+b]._around_the_world=wrap;
	  }
	  o +=_grid->_slice_stride[dimension];
	}

      } else { 
	
	int so  = plane*_grid->_ostride[dimension]; // base offset for start of plane 
	int o   = 0;                                      // relative offset to base within plane
	int bo  = 0;                                      // offset in buffer
    
	for(int n=0;n<_grid->_slice_nblock[dimension];n++){
	  for(int b=0;b<_grid->_slice_block[dimension];b++){

	    int ocb=1<<_grid->CheckerBoardFromOindex(o+b);// Could easily be a table lookup
	    if ( ocb & cbmask ) {
	      _entries[point][so+o+b]._offset  =offset+(bo++);
	      _entries[point][so+o+b]._is_local=0;
	      _entries[point][so+o+b]._permute =0;
	      _entries[point][so+o+b]._around_the_world=wrap;
	    }
	  }
	  o +=_grid->_slice_stride[dimension];
	}
      }
    }
}
