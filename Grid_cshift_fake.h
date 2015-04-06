#ifndef _GRID_FAKE_H_
#define _GRID_FAKE_H_

      

friend Lattice<vobj> Cshift(Lattice<vobj> &rhs,int dimension,int shift)
{
  typedef typename vobj::vector_type vector_type;
  typedef typename vobj::scalar_type scalar_type;

  Lattice<vobj> ret(rhs._grid);

  GridBase *grid=rhs._grid;
  const int Nsimd = grid->Nsimd();
  
  int fd = rhs._grid->_fdimensions[dimension];
  int rd = rhs._grid->_rdimensions[dimension];
  //int ld = rhs._grid->_ldimensions[dimension];
  //int gd = rhs._grid->_gdimensions[dimension];
  

  // Map to always positive shift modulo global full dimension.
  shift = (shift+fd)%fd;

  ret.checkerboard = rhs._grid->CheckerBoardDestination(rhs.checkerboard,shift);
        
  // the permute type

  int permute_dim =rhs._grid->_simd_layout[dimension]>1 ;
  int permute_type=0;
  for(int d=0;d<dimension;d++){
    if (rhs._grid->_simd_layout[d]>1 ) permute_type++;
  }

  ///////////////////////////////////////////////
  // Move via a fake comms buffer
  // Simd direction uses an extract/merge pair
  ///////////////////////////////////////////////
  int buffer_size = rhs._grid->_slice_nblock[dimension]*rhs._grid->_slice_block[dimension];
  int words = sizeof(vobj)/sizeof(vector_type);

  std::vector<vobj,alignedAllocator<vobj> > comm_buf(buffer_size);
  std::vector<std::vector<scalar_type> > comm_buf_extract(Nsimd,std::vector<scalar_type>(buffer_size*words) );
  std::vector<scalar_type *> pointers(Nsimd);

  for(int x=0;x<rd;x++){       

    for(int i=0;i<Nsimd;i++){
      pointers[i] = (scalar_type *)&comm_buf_extract[i][0];
    }

    int ro  = x*rhs._grid->_ostride[dimension]; // base offset for result

    if ( permute_dim ) {

      int o   = 0;                                // relative offset to base
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
	  } else {
	    ret._odata[ro+o+b]=rhs._odata[so+o+b];
	  }

	}
	o +=rhs._grid->_slice_stride[dimension];
      }

      for(int i=0;i<Nsimd;i++){
	pointers[i] = (scalar_type *)&comm_buf_extract[permute_map[permute_type][i]][0];
      }

      o   = 0;                                // relative offset to base
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

    } else {

      int co; // comm offset
      int o;

      co=0; o=0;
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
  }
  return ret;
}

/*
    friend Lattice<vobj> Cshift(Lattice<vobj> &rhs,int dimension,int shift)
    {
      Lattice<vobj> ret(rhs._grid);
        
        int rd = rhs._grid->_rdimensions[dimension];
        int ld = rhs._grid->_ldimensions[dimension];
        int gd = rhs._grid->_gdimensions[dimension];
        
        // Map to always positive shift.
        shift = (shift+gd)%gd;

        ret.checkerboard = rhs._grid->CheckerBoardDestination(rhs.checkerboard,shift);
        shift = rhs._grid->CheckerBoardShift(rhs.checkerboard,dimension,shift);
        
        // Work out whether to permute and the permute type
        // ABCDEFGH ->   AE BF CG DH       permute
        // Shift 0       AE BF CG DH       0 0 0 0    ABCDEFGH
        // Shift 1       BF CG DH AE       0 0 0 1    BCDEFGHA
        // Shift 2       CG DH AE BF       0 0 1 1    CDEFGHAB
        // Shift 3       DH AE BF CG       0 1 1 1    DEFGHABC
        // Shift 4       AE BF CG DH       1 1 1 1    EFGHABCD
        // Shift 5       BF CG DH AE       1 1 1 0    FGHACBDE
        // Shift 6       CG DH AE BF       1 1 0 0    GHABCDEF
        // Shift 7       DH AE BF CG       1 0 0 0    HABCDEFG

        int permute_dim =rhs._grid->_simd_layout[dimension]>1 ;
        int permute_type=0;
        for(int d=0;d<dimension;d++)
            if (rhs._grid->_simd_layout[d]>1 ) permute_type++;
        
        
        // loop over all work
        int work =rd*rhs._grid->_slice_nblock[dimension]*rhs._grid->_slice_block[dimension];

	// Packed gather sequence is clean
	int buffer_size = rhs._grid->_slice_nblock[dimension]*rhs._grid->_slice_block[dimension];
	
	typedef typename vobj::scalar_type scalar_t;
	typedef typename vobj::vector_type vector_t;
	const int ns=sizeof(vobj)/sizeof(scalar_t);
	const int nv=sizeof(vobj)/sizeof(vector_t);
	std::vector<vobj,alignedAllocator<vobj> > comm_buf(buffer_size);

        for(int x=0;x<rd;x++){       
	  
	  int sx = (x+shift)%rd;
	  int o  = x*rhs._grid->_ostride[dimension];
	  int so =sx*rhs._grid->_ostride[dimension];

	
	  int permute_slice=0;
	  if ( permute_dim ) {
	    permute_slice = shift/rd;
	    if ( x<shift%rd ) permute_slice = 1-permute_slice;
	  }

	  if ( permute_slice ) {
	    exit(0);
	    // For fake communication ALWAYS extract and either merge one way or other
	    scalar_t * bptr = (scalar_t *) &comm_buf[0];

	    int bo=0;
	    for(int n=0;n<rhs._grid->_slice_nblock[dimension];n++){
	    
	      vector_t *optr = (vector_t *)&ret._odata[o];
	      vector_t *iptr = (vector_t *)&rhs._odata[so];
	      int skew    = buffer_size*ns/2;

	      for(int b=0;b<rhs._grid->_slice_block[dimension];b++){ 
		for(int n=0;n<nv;n++){// number of simd vecsscalars in a vector
		  extract(iptr[b*nv+n],&bptr[n],skew,permute_type);
		}
	      }
	      o+=rhs._grid->_slice_stride[dimension];
	      //	      bo+=rhs._grid->_slice_stride[dimension]*ns/2;

	    }

	  } else {
	    int bo=0;
	    for(int n=0;n<rhs._grid->_slice_nblock[dimension];n++){
	      for(int i=0;i<rhs._grid->_slice_block[dimension];i++){
		comm_buf[bo++] =rhs._odata[so+i];
	      }
	      so+=rhs._grid->_slice_stride[dimension];
	    }
	    bo=0;
	    for(int n=0;n<rhs._grid->_slice_nblock[dimension];n++){
	      for(int i=0;i<rhs._grid->_slice_block[dimension];i++){
		ret._odata[o+i]=comm_buf[bo++];
	      }
	      o+=rhs._grid->_slice_stride[dimension];
	    }
	  }
	}
        return ret;
    };
*/

#endif
