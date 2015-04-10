#ifndef GRID_STENCIL_H
#define GRID_STENCIL_H
//////////////////////////////////////////////////////////////////////////////////////////
// Must not lose sight that goal is to be able to construct really efficient
// gather to a point stencil code. CSHIFT is not the best way, so probably need
// additional stencil support.
//
// Stencil based code could pre-exchange haloes and use a table lookup for neighbours
//
// Lattice <foo> could also allocate haloes which get used for stencil code.
//
// Grid could create a neighbour index table for a given stencil.
//
// Could also implement CovariantCshift, to fuse the loops and enhance performance.
//
//
// General stencil computation:
//
// Generic services
// 0) Prebuild neighbour tables
// 1) Compute sizes of all haloes/comms buffers; allocate them.
//
// 2) Gather all faces, and communicate.
// 3) Loop over result sites, giving nbr index/offnode info for each
// 
// Could take a 
// SpinProjectFaces 
// start comms
// complete comms 
// Reconstruct Umu
//
// Approach.
//
//////////////////////////////////////////////////////////////////////////////////////////

namespace Grid {

    class Stencil {
  public:

      Stencil(GridBase *grid,
	      int npoints,
	      int checkerboard,
	      std::vector<int> directions,
	      std::vector<int> distances);

      void Stencil_local     (int dimension,int shift,int cbmask);
      void Stencil_comms     (int dimension,int shift,int cbmask);
      void Stencil_comms_simd(int dimension,int shift,int cbmask);

      // Will need to implement actions for
      Copy_plane;
      Copy_plane_permute;
      Gather_plane;

      // The offsets to all neibours in stencil in each direction
      int                               _checkerboard;
      int                               _npoints; // Move to template param?
      GridBase *                        _grid;

      // Store these as SIMD Integer needed
      //
      // std::vector< iVector<Integer, Npoint> > _offsets;
      // std::vector< iVector<Integer, Npoint> > _local;
      // std::vector< iVector<Integer, Npoint> > _comm_buf_size;
      // std::vector< iVector<Integer, Npoint> > _permute;

      std::vector<std::vector<int>    > _offsets;
      std::vector<std::vector<int>    > _local;
      std::vector<int>                  _comm_buf_size;
      std::vector<int>                  _permute;

    };

    Stencil::Stencil(GridBase *grid,
		     int npoints,
		     int checkerboard,
		     std::vector<int> directions,
		     std::vector<int> distances){
      
      _npoints = npoints;
      _grid    = grid;
      
      for(int i=0;i<npoints;i++){

	int dimension    = directions[i];
	int displacement = distances[i];

	int fd = _grid->_fdimensions[dimension];
	int rd = _grid->_rdimensions[dimension];

	_checkerboard = checkerboard;

	// the permute type
	int simd_layout     = _grid->_simd_layout[dimension];
	int comm_dim        = _grid->_processors[dimension] >1 ;
	int splice_dim      = _grid->_simd_layout[dimension]>1 && (comm_dim);

	int sshift[2];

	if ( !comm_dim ) {
	  sshift[0] = _grid->CheckerBoardShift(_checkerboard,dimension,shift,0);
	  sshift[1] = _grid->CheckerBoardShift(_checkerboard,dimension,shift,1);

	  if ( sshift[0] == sshift[1] ) {
	    Stencil_local(dimension,shift,0x3);
	  } else {
	    Stencil_local(dimension,shift,0x1);// if checkerboard is unfavourable take two passes
	    Stencil_local(dimension,shift,0x2);// both with block stride loop iteration
	  }
	} else if ( splice_dim ) {
	  sshift[0] = _grid->CheckerBoardShift(_checkerboard,dimension,shift,0);
	  sshift[1] = _grid->CheckerBoardShift(_checkerboard,dimension,shift,1);
	  
	  if ( sshift[0] == sshift[1] ) {
	    Stencil_comms_simd(dimension,shift,0x3);
	  } else {
	    Stencil_comms_simd(dimension,shift,0x1);// if checkerboard is unfavourable take two passes
	    Stencil_comms_simd(dimension,shift,0x2);// both with block stride loop iteration
	  }
	} else {
	  //	  Cshift_comms(ret,rhs,dimension,shift);
	  sshift[0] = _grid->CheckerBoardShift(_checkerboard,dimension,shift,0);
	  sshift[1] = _grid->CheckerBoardShift(_checkerboard,dimension,shift,1);
	  if ( sshift[0] == sshift[1] ) {
	    Stencil_comms(dimension,shift,0x3);
	  } else {
	    Stencil_comms(dimension,shift,0x1);// if checkerboard is unfavourable take two passes
	    Stencil_comms(dimension,shift,0x2);// both with block stride loop iteration
	  }
	}
      }
    }


      void Stencil::Stencil_local     (int dimension,int shift,int cbmask)
      {
	int fd = _grid->_fdimensions[dimension];
	int rd = _grid->_rdimensions[dimension];
	int ld = _grid->_ldimensions[dimension];
	int gd = _grid->_gdimensions[dimension];

	// Map to always positive shift modulo global full dimension.
	shift = (shift+fd)%fd;
	
	// the permute type
	int permute_dim =_grid->PermuteDim(dimension);
	int permute_type=_grid->PermuteType(dimension);
	
	for(int x=0;x<rd;x++){       
	  
	  int o   = 0;
	  int bo  = x * _grid->_ostride[dimension];
	  
	  int cb= (cbmask==0x2)? 1 : 0;
	  
	  int sshift = _grid->CheckerBoardShift(_checkerboard,dimension,shift,cb);
	  int sx     = (x+sshift)%rd;
	  
	  int permute_slice=0;
	  if(permute_dim){
	    int wrap = sshift/rd;
	    int  num = sshift%rd;
	    if ( x< rd-num ) permute_slice=wrap;
	    else permute_slice = 1-wrap;
	  }

	  if ( permute_slice ) Copy_plane_permute(dimension,x,sx,cbmask,permute_type);
	  else                 Copy_plane        (dimension,x,sx,cbmask); 
  
	}
      }

      void Stencil::Stencil_comms     (int dimension,int shift,int cbmask)
      {
	typedef typename vobj::vector_type vector_type;
	typedef typename vobj::scalar_type scalar_type;
	
	GridBase *grid=_grid;

	int fd              = _grid->_fdimensions[dimension];
	int rd              = _grid->_rdimensions[dimension];
	int simd_layout     = _grid->_simd_layout[dimension];
	int comm_dim        = _grid->_processors[dimension] >1 ;

	assert(simd_layout==1);
	assert(comm_dim==1);
	assert(shift>=0);
	assert(shift<fd);
	
	int buffer_size = _grid->_slice_nblock[dimension]*rhs._grid->_slice_block[dimension];
	// FIXME: Do something with buffer_size??

	int cb= (cbmask==0x2)? 1 : 0;
	int sshift= _grid->CheckerBoardShift(_checkerboard,dimension,shift,cb);
	
	for(int x=0;x<rd;x++){       
	  
	  int offnode = ( x+sshift >= rd );
	  int sx        = (x+sshift)%rd;
	  int comm_proc = (x+sshift)/rd;
	  
	  if (!offnode) {
	    
	    Copy_plane(dimension,x,sx,cbmask); 
	    
	  } else {
	    
	    int words = send_buf.size();
	    if (cbmask != 0x3) words=words>>1;
	    
	    int bytes = words * sizeof(vobj);
	    
	    Gather_plane_simple (dimension,sx,cbmask);
	    
	    int rank           = grid->_processor;
	    int recv_from_rank;
	    int xmit_to_rank;
	    grid->ShiftedRanks(dimension,comm_proc,xmit_to_rank,recv_from_rank);
	    /*	    
	    grid->SendToRecvFrom((void *)&send_buf[0],
				 xmit_to_rank,
				 (void *)&recv_buf[0],
				 recv_from_rank,
				 bytes);
	    */
	    Scatter_plane_simple (dimension,x,cbmask);
	  }
	}
      }

      void Stencil::Stencil_comms_simd(int dimension,int shift,int cbmask)
      {
	GridBase *grid=_grid;
	const int Nsimd = _grid->Nsimd();
	typedef typename vobj::vector_type vector_type;
	typedef typename vobj::scalar_type scalar_type;
	
	int fd = _grid->_fdimensions[dimension];
	int rd = _grid->_rdimensions[dimension];
	int ld = _grid->_ldimensions[dimension];
	int simd_layout     = _grid->_simd_layout[dimension];
	int comm_dim        = _grid->_processors[dimension] >1 ;
	
	assert(comm_dim==1);
	assert(simd_layout==2);
	assert(shift>=0);
	assert(shift<fd);
	
	int permute_type=_grid->PermuteType(dimension);
	
	///////////////////////////////////////////////
	// Simd direction uses an extract/merge pair
	///////////////////////////////////////////////
	int buffer_size = _grid->_slice_nblock[dimension]*grid->_slice_block[dimension];
	// FIXME do something with buffer size
	
	std::vector<scalar_type *> pointers(Nsimd);  // 
	std::vector<scalar_type *> rpointers(Nsimd); // received pointers
	
	///////////////////////////////////////////
	// Work out what to send where
	///////////////////////////////////////////
	
	int cb    = (cbmask==0x2)? 1 : 0;
	int sshift= _grid->CheckerBoardShift(_checkerboard,dimension,shift,cb);
	
	std::vector<int> comm_offnode(simd_layout);
	std::vector<int> comm_proc   (simd_layout);  //relative processor coord in dim=dimension
	std::vector<int> icoor(grid->Nd());
	
	for(int x=0;x<rd;x++){       
	  
	  int comm_any = 0;
	  for(int s=0;s<simd_layout;s++) {
	    int shifted_x   = x+s*rd+sshift;
	    comm_offnode[s] = shifted_x >= ld; 
	    comm_any        = comm_any | comm_offnode[s];
	    comm_proc[s]    = shifted_x/ld;     
	  }
	  
	  int o    = 0;
	  int bo   = x*grid->_ostride[dimension];
	  int sx   = (x+sshift)%rd;
	  
	  if ( comm_any ) {
	    
	    for(int i=0;i<Nsimd;i++){
	      pointers[i] = (scalar_type *)&send_buf_extract[i][0];
	    }
	    Gather_plane_extract(rhs,pointers,dimension,sx,cbmask);
	    
	    for(int i=0;i<Nsimd;i++){
	      
	      int s;
	      grid->iCoorFromIindex(icoor,i);
	      s = icoor[dimension];
	      
	      if(comm_offnode[s]){
		
		int rank           = grid->_processor;
		int recv_from_rank;
		int xmit_to_rank;
		grid->ShiftedRanks(dimension,comm_proc[s],xmit_to_rank,recv_from_rank);
		
		/*		
		grid->SendToRecvFrom((void *)&send_buf_extract[i][0],
				     xmit_to_rank,
				     (void *)&recv_buf_extract[i][0],
				     recv_from_rank,
				     bytes);
		*/

		rpointers[i] = (scalar_type *)&recv_buf_extract[i][0];
		
	      } else { 
		
		rpointers[i] = (scalar_type *)&send_buf_extract[i][0];

	      }
	      
	    }
	    
	    // Permute by swizzling pointers in merge
	    int permute_slice=0;
	    int lshift=sshift%ld;
	    int wrap  =lshift/rd;
	    int  num  =lshift%rd;

	    if ( x< rd-num ) permute_slice=wrap;
	    else permute_slice = 1-wrap;
	    
	    int toggle_bit = (Nsimd>>(permute_type+1));
	    int PermuteMap;
	    for(int i=0;i<Nsimd;i++){
	      if ( permute_slice ) {
		PermuteMap=i^toggle_bit;
		pointers[i] = rpointers[PermuteMap];
	      } else {
		pointers[i] = rpointers[i];
	      }
	    }

	    Scatter_plane_merge(pointers,dimension,x,cbmask);
	    
	  } else { 

	    int permute_slice=0;
	    int wrap = sshift/rd;
	    int  num = sshift%rd;
	    if ( x< rd-num ) permute_slice=wrap;
	    else permute_slice = 1-wrap;

	    if ( permute_slice ) Copy_plane_permute(ret,rhs,dimension,x,sx,cbmask,permute_type);
	    else                 Copy_plane(ret,rhs,dimension,x,sx,cbmask); 

	  }
	}
      }


};
#endif
