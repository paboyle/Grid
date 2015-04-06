#ifndef _GRID_CSHIFT_COMMON_H_
#define _GRID_CSHIFT_COMMON_H_
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
// Could also implement CovariantCshift.
//////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////
// Gather for when there is no need to SIMD split
//////////////////////////////////////////////////////
friend void Gather_plane_simple (Lattice<vobj> &rhs,std::vector<vobj,alignedAllocator<vobj> > &buffer,             int dimension,int plane,int cbmask)
{
  int rd = rhs._grid->_rdimensions[dimension];

  if ( !rhs._grid->CheckerBoarded(dimension) ) {

    int so  = plane*rhs._grid->_ostride[dimension]; // base offset for start of plane 
    int o   = 0;                                    // relative offset to base within plane
    int bo  = 0;                                    // offset in buffer

    // Simple block stride gather of SIMD objects
#pragma omp parallel for collapse(2)
    for(int n=0;n<rhs._grid->_slice_nblock[dimension];n++){
      for(int b=0;b<rhs._grid->_slice_block[dimension];b++){
	buffer[bo++]=rhs._odata[so+o+b];
      }
      o +=rhs._grid->_slice_stride[dimension];
    }

  } else { 

    int so  = plane*rhs._grid->_ostride[dimension]; // base offset for start of plane 
    int o   = 0;                                      // relative offset to base within plane
    int bo  = 0;                                      // offset in buffer

#pragma omp parallel for collapse(2)
    for(int n=0;n<rhs._grid->_slice_nblock[dimension];n++){
      for(int b=0;b<rhs._grid->_slice_block[dimension];b++){

	int ocb=1<<rhs._grid->CheckerBoardFromOindex(o+b);// Could easily be a table lookup
	if ( ocb &cbmask ) {
	  buffer[bo]=rhs._odata[so+o+b];
	  bo++;
	}

      }
      o +=rhs._grid->_slice_stride[dimension];
    }
  }
}


//////////////////////////////////////////////////////
// Gather for when there *is* need to SIMD split
//////////////////////////////////////////////////////
friend void Gather_plane_extract(Lattice<vobj> &rhs,std::vector<scalar_type *> pointers,int dimension,int plane,int cbmask)
{
  int rd = rhs._grid->_rdimensions[dimension];

  if ( !rhs._grid->CheckerBoarded(dimension) ) {

    int so  = plane*rhs._grid->_ostride[dimension]; // base offset for start of plane 
    int o   = 0;                                    // relative offset to base within plane
    int bo  = 0;                                    // offset in buffer

    // Simple block stride gather of SIMD objects
#pragma omp parallel for collapse(2)
    for(int n=0;n<rhs._grid->_slice_nblock[dimension];n++){
      for(int b=0;b<rhs._grid->_slice_block[dimension];b++){
	extract(rhs._odata[so+o+b],pointers);
      }
      o +=rhs._grid->_slice_stride[dimension];
    }

  } else { 

    int so  = plane*rhs._grid->_ostride[dimension]; // base offset for start of plane 
    int o   = 0;                                      // relative offset to base within plane
    int bo  = 0;                                      // offset in buffer
    
#pragma omp parallel for collapse(2)
    for(int n=0;n<rhs._grid->_slice_nblock[dimension];n++){
      for(int b=0;b<rhs._grid->_slice_block[dimension];b++){

	int ocb=1<<rhs._grid->CheckerBoardFromOindex(o+b);
	if ( ocb & cbmask ) {
	  extract(rhs._odata[so+o+b],pointers);
	}

      }
      o +=rhs._grid->_slice_stride[dimension];
    }
  }
}



//////////////////////////////////////////////////////
// Scatter for when there is no need to SIMD split
//////////////////////////////////////////////////////
friend void Scatter_plane_simple (Lattice<vobj> &rhs,std::vector<vobj,alignedAllocator<vobj> > &buffer,             int dimension,int plane,int cbmask)
{
  int rd = rhs._grid->_rdimensions[dimension];

  if ( !rhs._grid->CheckerBoarded(dimension) ) {

    int so  = plane*rhs._grid->_ostride[dimension]; // base offset for start of plane 
    int o   = 0;                                    // relative offset to base within plane
    int bo  = 0;                                    // offset in buffer

    // Simple block stride gather of SIMD objects
#pragma omp parallel for collapse(2)
    for(int n=0;n<rhs._grid->_slice_nblock[dimension];n++){
      for(int b=0;b<rhs._grid->_slice_block[dimension];b++){
	rhs._odata[so+o+b]=buffer[bo++];
      }
      o +=rhs._grid->_slice_stride[dimension];
    }

  } else { 

    int so  = plane*rhs._grid->_ostride[dimension]; // base offset for start of plane 
    int o   = 0;                                      // relative offset to base within plane
    int bo  = 0;                                      // offset in buffer
    
#pragma omp parallel for collapse(2)
    for(int n=0;n<rhs._grid->_slice_nblock[dimension];n++){
      for(int b=0;b<rhs._grid->_slice_block[dimension];b++){

	int ocb=1<<rhs._grid->CheckerBoardFromOindex(o+b);// Could easily be a table lookup
	if ( ocb & cbmask ) {
	  rhs._odata[so+o+b]=buffer[bo++];
	}

      }
      o +=rhs._grid->_slice_stride[dimension];
    }
  }
}


//////////////////////////////////////////////////////
// Scatter for when there *is* need to SIMD split
//////////////////////////////////////////////////////
friend void Scatter_plane_merge(Lattice<vobj> &rhs,std::vector<scalar_type *> pointers,int dimension,int plane,int cbmask)
{
  int rd = rhs._grid->_rdimensions[dimension];

  if ( !rhs._grid->CheckerBoarded(dimension) ) {

    int so  = plane*rhs._grid->_ostride[dimension]; // base offset for start of plane 
    int o   = 0;                                    // relative offset to base within plane
    int bo  = 0;                                    // offset in buffer

    // Simple block stride gather of SIMD objects
#pragma omp parallel for collapse(2)
    for(int n=0;n<rhs._grid->_slice_nblock[dimension];n++){
      for(int b=0;b<rhs._grid->_slice_block[dimension];b++){
	merge(rhs._odata[so+o+b],pointers);
      }
      o +=rhs._grid->_slice_stride[dimension];
    }

  } else { 

    int so  = plane*rhs._grid->_ostride[dimension]; // base offset for start of plane 
    int o   = 0;                                      // relative offset to base within plane
    int bo  = 0;                                      // offset in buffer
    
#pragma omp parallel for collapse(2)
    for(int n=0;n<rhs._grid->_slice_nblock[dimension];n++){
      for(int b=0;b<rhs._grid->_slice_block[dimension];b++){

	int ocb=1<<rhs._grid->CheckerBoardFromOindex(o+b);
	if ( ocb&cbmask ) {
	  merge(rhs._odata[so+o+b],pointers);
	}

      }
      o +=rhs._grid->_slice_stride[dimension];
    }
  }
}


//////////////////////////////////////////////////////
// local to node block strided copies
//////////////////////////////////////////////////////
// if lhs is odd, rhs even??
friend void Copy_plane(Lattice<vobj>& lhs,Lattice<vobj> &rhs, int dimension,int lplane,int rplane,int cbmask)
{
  int rd = rhs._grid->_rdimensions[dimension];

  if ( !rhs._grid->CheckerBoarded(dimension) ) {

    int o   = 0;                                     // relative offset to base within plane
    int ro  = rplane*rhs._grid->_ostride[dimension]; // base offset for start of plane 
    int lo  = lplane*lhs._grid->_ostride[dimension]; // offset in buffer

  // Simple block stride gather of SIMD objects
#pragma omp parallel for collapse(2)
    for(int n=0;n<rhs._grid->_slice_nblock[dimension];n++){
      for(int b=0;b<rhs._grid->_slice_block[dimension];b++){
	lhs._odata[lo+o+b]=rhs._odata[ro+o+b];
      }
      o +=rhs._grid->_slice_stride[dimension];
    }

  } else {

    int ro  = rplane*rhs._grid->_ostride[dimension]; // base offset for start of plane 
    int lo  = lplane*lhs._grid->_ostride[dimension]; // base offset for start of plane 
    int o   = 0;                                     // relative offset to base within plane

#pragma omp parallel for collapse(2)
    for(int n=0;n<rhs._grid->_slice_nblock[dimension];n++){
      for(int b=0;b<rhs._grid->_slice_block[dimension];b++){

	int ocb=1<<lhs._grid->CheckerBoardFromOindex(o+b);

	if ( ocb&cbmask ) {
	  lhs._odata[lo+o+b]=rhs._odata[ro+o+b];
	}

      }
      o +=rhs._grid->_slice_stride[dimension];
    }

  }
}

friend void Copy_plane_permute(Lattice<vobj>& lhs,Lattice<vobj> &rhs, int dimension,int lplane,int rplane,int cbmask,int permute_type)
{
  int rd = rhs._grid->_rdimensions[dimension];


  if ( !rhs._grid->CheckerBoarded(dimension) ) {

    int o   = 0;                                     // relative offset to base within plane
    int ro  = rplane*rhs._grid->_ostride[dimension]; // base offset for start of plane 
    int lo  = lplane*rhs._grid->_ostride[dimension]; // offset in buffer

  // Simple block stride gather of SIMD objects
#pragma omp parallel for collapse(2)
    for(int n=0;n<rhs._grid->_slice_nblock[dimension];n++){
      for(int b=0;b<rhs._grid->_slice_block[dimension];b++){
	permute(lhs._odata[lo+o+b],rhs._odata[ro+o+b],permute_type);
      }
      o +=rhs._grid->_slice_stride[dimension];
    }

  } else {

    int ro  = rplane*rhs._grid->_ostride[dimension]; // base offset for start of plane 
    int lo  = lplane*lhs._grid->_ostride[dimension]; // base offset for start of plane 
    int o   = 0;                                     // relative offset to base within plane
    
#pragma omp parallel for collapse(2)
    for(int n=0;n<rhs._grid->_slice_nblock[dimension];n++){
      for(int b=0;b<rhs._grid->_slice_block[dimension];b++){

	int ocb=1<<lhs._grid->CheckerBoardFromOindex(o+b);

	if ( ocb&cbmask ) {
	  permute(lhs._odata[lo+o+b],rhs._odata[ro+o+b],permute_type);
	}

      }
      o +=rhs._grid->_slice_stride[dimension];
    }

  }
}

//////////////////////////////////////////////////////
// Local to node Cshift
//////////////////////////////////////////////////////

  // Work out whether to permute 
  // ABCDEFGH ->   AE BF CG DH       permute              wrap num
  //
  // Shift 0       AE BF CG DH       0 0 0 0    ABCDEFGH   0   0
  // Shift 1       BF CG DH AE       0 0 0 1    BCDEFGHA   0   1
  // Shift 2       CG DH AE BF       0 0 1 1    CDEFGHAB   0   2
  // Shift 3       DH AE BF CG       0 1 1 1    DEFGHABC   0   3
  // Shift 4       AE BF CG DH       1 1 1 1    EFGHABCD   1   0 
  // Shift 5       BF CG DH AE       1 1 1 0    FGHACBDE   1   1 
  // Shift 6       CG DH AE BF       1 1 0 0    GHABCDEF   1   2
  // Shift 7       DH AE BF CG       1 0 0 0    HABCDEFG   1   3

  // Suppose 4way simd in one dim.
  // ABCDEFGH ->   AECG BFDH      permute              wrap num

  // Shift 0       AECG BFDH      0,00 0,00 ABCDEFGH         0     0
  // Shift 1       BFDH CGEA      0,00 1,01 BCDEFGHA         0     1
  // Shift 2       CGEA DHFB      1,01 1,01 CDEFGHAB         1     0
  // Shift 3       DHFB EAGC      1,01 1,11 DEFGHABC         1     1
  // Shift 4       EAGC FBHD      1,11 1,11 EFGHABCD         2     0 
  // Shift 5       FBHD GCAE      1,11 1,10 FGHABCDE         2     1
  // Shift 6       GCAE HDBF      1,10 1,10 GHABCDEF         3     0
  // Shift 7       HDBF AECG      1,10 0,00 HABCDEFG         3     1

  // Generalisation to 8 way simd, 16 way simd required.
  //
  // Need log2 Nway masks. consisting of 
  //	    1 bit  256 bit granule
  //	    2 bit  128 bit granule
  //        4 bits 64  bit granule
  //        8 bits 32  bit granules
  //
  //        15 bits....

friend void Cshift_local(Lattice<vobj>& ret,Lattice<vobj> &rhs,int dimension,int shift)
{
  int sshift[2];

  sshift[0] = rhs._grid->CheckerBoardShift(rhs.checkerboard,dimension,shift,0);
  sshift[1] = rhs._grid->CheckerBoardShift(rhs.checkerboard,dimension,shift,1);

  if ( sshift[0] == sshift[1] ) {
    Cshift_local(ret,rhs,dimension,shift,0x3);
  } else {
    Cshift_local(ret,rhs,dimension,shift,0x1);// if checkerboard is unfavourable take two passes
    Cshift_local(ret,rhs,dimension,shift,0x2);// both with block stride loop iteration
  }
}


friend Lattice<vobj> Cshift_local(Lattice<vobj> &ret,Lattice<vobj> &rhs,int dimension,int shift,int cbmask)
{
  int fd = rhs._grid->_fdimensions[dimension];
  int rd = rhs._grid->_rdimensions[dimension];
  int ld = rhs._grid->_ldimensions[dimension];
  int gd = rhs._grid->_gdimensions[dimension];
  

  // Map to always positive shift modulo global full dimension.
  shift = (shift+fd)%fd;

  ret.checkerboard = rhs._grid->CheckerBoardDestination(rhs.checkerboard,shift);
        
  // the permute type
  int permute_dim =rhs._grid->_simd_layout[dimension]>1 ;
  int permute_type=0;
  for(int d=0;d<dimension;d++){
    if (rhs._grid->_simd_layout[d]>1 ) permute_type++;
  }

  for(int x=0;x<rd;x++){       

    int o   = 0;
    int bo  = x * rhs._grid->_ostride[dimension];
    
    int cb= (cbmask==0x2)? 1 : 0;

    int sshift = rhs._grid->CheckerBoardShift(rhs.checkerboard,dimension,shift,cb);
    int sx     = (x+sshift)%rd;
	
    int permute_slice=0;
    if(permute_dim){
      int wrap = sshift/rd;
      int  num = sshift%rd;
      if ( x< rd-num ) permute_slice=wrap;
      else permute_slice = 1-wrap;
    }

    if ( permute_slice ) Copy_plane_permute(ret,rhs,dimension,x,sx,cbmask,permute_type);
    else                 Copy_plane(ret,rhs,dimension,x,sx,cbmask); 

  
  }
  return ret;
}

#endif
