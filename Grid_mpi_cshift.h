#ifndef _GRID_MPI_CSHIFT_H_
#define _GRID_MPI_CSHIFT_H_

      
//////////////////////////////////////////////
// Q. Split this into seperate sub functions?
//////////////////////////////////////////////

// CshiftCB_comms_splice
// CshiftCB_comms
// CshiftCB_local
// CshiftCB_local_permute

// Cshift_comms_splice
// Cshift_comms
// Cshift_local
// Cshift_local_permute

// Broadly I remain annoyed that the iteration is so painful
// for red black data layout, when simple block strided descriptors suffice for non-cb. 
//
// The other option is to do it table driven, or perhaps store the CB of each site in a table.
//
// Must not lose sight that goal is to be able to construct really efficient
// gather to a point stencil code. CSHIFT is not the best way, so probably need
// additional stencil support.
//
// Could still do a templated syntax tree and make CSHIFT return lattice vector.
//
// Stencil based code could pre-exchange haloes and use a table lookup for neighbours
//
// Lattice <foo> could also allocate haloes which get used for stencil code.
//
// Grid could create a neighbour index table for a given stencil.
// Could also implement CovariantCshift.

//////////////////////////////////////////////////////
//Non checkerboarded support functions
//////////////////////////////////////////////////////

friend void Gather_plane        (Lattice<vobj> &rhs,std::vector<vobj> &buffer,             int dimension,int plane)
{
  const int Nsimd = vector_type::Nsimd();
  int rd = rhs._grid->_rdimensions[dimension];

  int so  = plane*rhs._grid->_ostride[dimension]; // base offset for start of plane 
  int o   = 0;                                    // relative offset to base within plane
#pragma omp parallel for collapse(2)
  for(int n=0;n<rhs._grid->_slice_nblock[dimension];n++){
    for(int b=0;b<rhs._grid->_slice_block[dimension];b++){
	  
	  int sshift = rhs._grid->CheckerBoardShift(rhs.checkerboard,dimension,shift,o+b);
	  int sx = (x+sshift)%rd;
	  int so = sx*rhs._grid->_ostride[dimension];
	  int permute_slice=0;
	  int wrap = sshift/rd;
	  int  num = sshift%rd;
	  
	  if ( x< rd-num ) permute_slice=wrap;
	  else permute_slice = 1-wrap;
	  
	  if ( permute_slice ) {
	    permute(ret._odata[ro+o+b],rhs._odata[so+o+b],permute_type);
	  } else {
	    ret._odata[ro+o+b]=rhs._odata[so+o+b];
	  }
	}
	o +=rhs._grid->_slice_stride[dimension];
      }
  }
  
}
//friend void Gather_plane_extract(Lattice<vobj> &rhs,std::vector<scalar_type *> pointers,int dimension,int plane);
//
//friend void Scatter_plane       (Lattice<vobj> &rhs,std::vector<vobj> face,             int dimension,int plane);
//friend void Scatter_plane_merge (Lattice<vobj> &rhs,std::vector<scalar_type *> pointers,int dimension,int plane);
//
//template<int permute_type> friend void Copy_plane_permute(Lattice<vobj> &rhs,std::vector<vobj> face, int dimension,int plane);
//                           friend void Copy_plane(Lattice<vobj> &rhs,std::vector<vobj> face, int dimension,int plane);
//

//////////////////////////////////////////////////////
//Checkerboarded support functions
//////////////////////////////////////////////////////

//friend void GatherCB_plane        (Lattice<vobj> &rhs,std::vector<vobj> face,             int dimension,int plane);
//friend void GatherCB_plane_extract(Lattice<vobj> &rhs,std::vector<scalar_type *> pointers,int dimension,int plane);
//
//friend void ScatterCB_plane       (Lattice<vobj> &rhs,std::vector<vobj> face,             int dimension,int plane);
//friend void ScatterCB_plane_merge (Lattice<vobj> &rhs,std::vector<scalar_type *> pointers,int dimension,int plane);
//
//template<int permute_type> friend void CopyCB_plane_permute(Lattice<vobj> &rhs,std::vector<vobj> face, int dimension,int plane);
//                           friend void Copy_plane(Lattice<vobj> &rhs,std::vector<vobj> face, int dimension,int plane);


friend Lattice<vobj> Cshift(Lattice<vobj> &rhs,int dimension,int shift)
{
  typedef typename vobj::vector_type vector_type;
  typedef typename vobj::scalar_type scalar_type;
  const int Nsimd = vector_type::Nsimd();

  Lattice<vobj> ret(rhs._grid);
  
  int fd = rhs._grid->_fdimensions[dimension];
  int rd = rhs._grid->_rdimensions[dimension];
  //int ld = rhs._grid->_ldimensions[dimension];
  //int gd = rhs._grid->_gdimensions[dimension];
  

  // Map to always positive shift modulo global full dimension.
  shift = (shift+fd)%fd;

  ret.checkerboard = rhs._grid->CheckerBoardDestination(rhs.checkerboard,shift);
        
  // the permute type
  int simd_layout     = rhs._grid->_simd_layout[dimension];
  int comm_dim        = rhs._grid->_processors[dimension] >1 ;
  int permute_dim     = rhs._grid->_simd_layout[dimension]>1 && (!comm_dim);
  int splice_dim      = rhs._grid->_simd_layout[dimension]>1 && (comm_dim);

  int permute_type=0;
  for(int d=0;d<dimension;d++){
    if (rhs._grid->_simd_layout[d]>1 ) permute_type++;
  }

  // Logic for non-distributed dimension	  
  std::vector<int> comm_offnode(simd_layout);
  std::vector<int> comm_to     (simd_layout);
  std::vector<int> comm_from   (simd_layout);
  std::vector<int> comm_rx     (simd_layout);  // reduced coordinate of neighbour plane
  std::vector<int> comm_simd_lane(simd_layout);// simd lane of neigbour plane

  ///////////////////////////////////////////////
  // Move via a fake comms buffer
  // Simd direction uses an extract/merge pair
  ///////////////////////////////////////////////
  int buffer_size = rhs._grid->_slice_nblock[dimension]*rhs._grid->_slice_block[dimension];
  int words = sizeof(vobj)/sizeof(vector_type);

  std::vector<vobj,alignedAllocator<vobj> > comm_buf(buffer_size);
  std::vector<std::vector<scalar_type> > comm_buf_extract(Nsimd,std::vector<scalar_type>(buffer_size*words) );
  std::vector<scalar_type *> pointers(Nsimd);



  if ( permute_dim ) {

    for(int x=0;x<rd;x++){       
      int ro  = x*rhs._grid->_ostride[dimension]; // base offset for result
      int o   = 0;                                // relative offset to base
      for(int n=0;n<rhs._grid->_slice_nblock[dimension];n++){
	for(int b=0;b<rhs._grid->_slice_block[dimension];b++){
	  
	  int sshift = rhs._grid->CheckerBoardShift(rhs.checkerboard,dimension,shift,o+b);
	  int sx = (x+sshift)%rd;
	  int so = sx*rhs._grid->_ostride[dimension];
	  int permute_slice=0;
	  int wrap = sshift/rd;
	  int  num = sshift%rd;
	  
	  if ( x< rd-num ) permute_slice=wrap;
	  else permute_slice = 1-wrap;
	  
	  if ( permute_slice ) {
	    permute(ret._odata[ro+o+b],rhs._odata[so+o+b],permute_type);
	  } else {
	    ret._odata[ro+o+b]=rhs._odata[so+o+b];
	  }
	}
	o +=rhs._grid->_slice_stride[dimension];
      }
    }

  } else if ( splice_dim ) {

    if ( rhs._grid->_simd_layout[dimension] > 2 ) exit(-1); // use Cassert. Audit code for exit and replace
    if ( rhs._grid->_simd_layout[dimension] < 1 ) exit(-1);
	    

    for(int i=0;i<vobj::vector_type::Nsimd();i++){
      pointers[i] = (scalar_type *)&comm_buf_extract[i][0];
    }

    
    for(int x=0;x<rd;x++){       


      ///////////////////////////////////////////
      // Extract one orthogonal slice at a time
      ///////////////////////////////////////////
      int ro  = x*rhs._grid->_ostride[dimension]; // base offset for result

      o   = 0;                                    // relative offset to base      
      for(int n=0;n<rhs._grid->_slice_nblock[dimension];n++){
	for(int b=0;b<rhs._grid->_slice_block[dimension];b++){
	  
	  int sshift = rhs._grid->CheckerBoardShift(rhs.checkerboard,dimension,shift,o+b);
	  int sx = (x+sshift)%rd;
	  
	  // base offset for source
	  int so = sx*rhs._grid->_ostride[dimension];

	  int permute_slice=0;
	  int wrap = sshift/rd;
	  int  num = sshift%rd;
	    
	  if ( x< rd-num ) permute_slice=wrap;
	  else permute_slice = 1-wrap;

	  if ( permute_slice ) {
	    extract(rhs._odata[so+o+b],pointers);
	  }
	}
	o +=rhs._grid->_slice_stride[dimension];
      }

    
      ///////////////////////////////////////////
      // Work out what to send where
      ///////////////////////////////////////////

      for(int s=0;s<simd_layout;s++) {


	// shift to "neighbour" takes us off node
	// coordinates (rx, simd_lane) of neighbour
	// how many nodes away is this shift
	// where we should send to
	// where we should receive from
	int shifted_x = x+s*rd+shift;
	comm_offnode[s]           = shifted_x > ld; 
	comm_send_rx[s]           = shifted_x%rd;     // which slice geton the other node
	comm_send_simd_lane [s]   = shifted_x/rd;     // which slice on the other node
	comm_from[s]    = shifted_x/ld;     
	comm_to  [s] = (2*_processors[dimension]-comm_from[s]) % _processors[dimension];
	comm_from[s] = (comm_from[s]+_processors[dimension])   % _processors[dimension];
	      
      }

      ////////////////////////////////////////////////
      // Insert communication phase
      ////////////////////////////////////////////////
#if 0
	  } else if (comm_dim ) {
	    
	    // Packed gather sequence is clean
	    int buffer_size = rhs._grid->_slice_nblock[dimension]*rhs._grid->_slice_nblock[dimension];
	    std::vector<vobj,alignedAllocator<vobj> > send_buf(buffer_size);
	    std::vector<vobj,alignedAllocator<vobj> > recv_buf(buffer_size);

	    // off node; communcate slice (ld==rd)
	    if ( x+shift > rd ) { 
	      int sb=0;
	      for(int n=0;n<rhs._grid->_slice_nblock[dimension];n++){
		for(int i=0;i<rhs._grid->_slice_block[dimension];i++){
		  send_buf[sb++]=rhs._odata[so+i];
		}
		so+=rhs._grid->_slice_stride[dimension];
	      }	      

	      // Make a comm_fake them mimics comms in periodic case.
	      
	      // scatter face
	      int rb=0;
	      for(int n=0;n<rhs._grid->_slice_nblock[dimension];n++){
		for(int i=0;i<rhs._grid->_slice_block[dimension];i++){
		  ret._odata[so+i]=recv_buf[rb++];
		}
		so+=rhs._grid->_slice_stride[dimension];
	      }	      

	    } else { 

	      for(int n=0;n<rhs._grid->_slice_nblock[dimension];n++){
		for(int i=0;i<rhs._grid->_slice_block[dimension];i++){
		  ret._odata[o+i]=rhs._odata[so+i];
		}
		o+=rhs._grid->_slice_stride[dimension];
		so+=rhs._grid->_slice_stride[dimension];
	      }
	    }

#endif	
    
      ////////////////////////////////////////////////
      // Pull receive buffers and permuted buffers in
      ////////////////////////////////////////////////
      for(int i=0;i<vobj::vector_type::Nsimd();i++){
	pointers[i] = (scalar_type *)&comm_buf_extract[permute_map[permute_type][i]][0];
      }

      o   = 0;                                    // relative offset to base
      int ro  = x*rhs._grid->_ostride[dimension]; // base offset for result
      for(int n=0;n<rhs._grid->_slice_nblock[dimension];n++){
	for(int b=0;b<rhs._grid->_slice_block[dimension];b++){
	  
	  int sshift = rhs._grid->CheckerBoardShift(rhs.checkerboard,dimension,shift,o+b);
	  int sx = (x+sshift)%rd;
	  
	  // base offset for source
	  int so = sx*rhs._grid->_ostride[dimension];
	  
	  int permute_slice=0;
	  int wrap = sshift/rd;
	  int  num = sshift%rd;
	  
	  if ( x< rd-num ) permute_slice=wrap;
	  else permute_slice = 1-wrap;
	  
	  if ( permute_slice ) {
	    merge(ret._odata[ro+o+b],pointers);
	  }
	}
	o +=rhs._grid->_slice_stride[dimension];
      }
    }
  

  } else if ( comm_dim ) { 
	

    int co; // comm offset
    int o;
    
    co=0;
    for(int x=0;x<rd;x++){       
      o=0;
      int ro  = x*rhs._grid->_ostride[dimension]; // base offset for result
      for(int n=0;n<rhs._grid->_slice_nblock[dimension];n++){
	for(int b=0;b<rhs._grid->_slice_block[dimension];b++){

	  // This call in inner loop is annoying but necessary for dimension=0
	  // in the case of RedBlack grids. Could optimise away with 
	  // alternate code paths for all other cases.
	  int sshift = rhs._grid->CheckerBoardShift(rhs.checkerboard,dimension,shift,o+b);
	  int sx = (x+sshift)%rd;
	  int so = sx*rhs._grid->_ostride[dimension];

	  comm_buf[co++]=rhs._odata[so+o+b];

	}
	o +=rhs._grid->_slice_stride[dimension];
      }

      // Step through a copy into a comms buffer and pull back in.
      // Genuine fake implementation could calculate if loops back
      co=0; o=0;
      for(int n=0;n<rhs._grid->_slice_nblock[dimension];n++){
	for(int b=0;b<rhs._grid->_slice_block[dimension];b++){
	  ret._odata[ro+o+b]=comm_buf[co++];
	}
	o +=rhs._grid->_slice_stride[dimension];
      }
    }

  } else { // Local dimension, no permute required

    for(int x=0;x<rd;x++){       
      int o=0;
      int ro  = x*rhs._grid->_ostride[dimension]; // base offset for result
      for(int n=0;n<rhs._grid->_slice_nblock[dimension];n++){
	for(int b=0;b<rhs._grid->_slice_block[dimension];b++){

	  // This call in inner loop is annoying but necessary for dimension=0
	  // in the case of RedBlack grids. Could optimise away with 
	  // alternate code paths for all other cases.
	  int sshift = rhs._grid->CheckerBoardShift(rhs.checkerboard,dimension,shift,o+b);
	  
	  int sx = (x+sshift)%rd;
	  int so = sx*rhs._grid->_ostride[dimension];
	  ret._odata[bo+o+b]=rhs._odata[so+o+b];

	}
	o +=rhs._grid->_slice_stride[dimension];
      }
    }

  }

  return ret;
}


#endif
