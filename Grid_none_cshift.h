#ifndef _GRID_NONE_CSHIFT_H_
#define _GRID_NONE_CSHIFT_H_
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
  
// For optimisation:
//
// split into Cshift_none_rb_permute
// split into Cshift_none_rb_simple
//
// split into Cshift_none_permute
// split into Cshift_none_simple

friend Lattice<vobj> Cshift(Lattice<vobj> &rhs,int dimension,int shift)
{
  Lattice<vobj> ret(rhs._grid);
  
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

    int bo  = x*rhs._grid->_ostride[dimension];
    int o   = 0;

    if ( permute_dim ) {

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
	    permute(ret._odata[bo+o+b],rhs._odata[so+o+b],permute_type);
	  } else {
	    ret._odata[bo+o+b]=rhs._odata[so+o+b];
	  }

	}
	o +=rhs._grid->_slice_stride[dimension];

      }
    } else {
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
        
#if 0 
// Collapse doesn't appear to work the way I think it should in icpc
friend Lattice<vobj> Cshift(Lattice<vobj> &rhs,int dimension,int shift)
{
  Lattice<vobj> ret(rhs._grid);
        
  ret.checkerboard = rhs._grid->CheckerBoardDestination(rhs.checkerboard,shift);
  shift = rhs._grid->CheckerBoardShift(rhs.checkerboard,dimension,shift);
  int sx,so,o;
  int rd = rhs._grid->_rdimensions[dimension];
  int ld = rhs._grid->_dimensions[dimension];
  // Map to always positive shift.
  shift = (shift+ld)%ld;
  // Work out whether to permute and the permute type
  // ABCDEFGH ->   AE BF CG DH       permute
  // Shift 0       AE BF CG DH       0 0 0 0    ABCDEFGH
  // Shift 1       DH AE BF CG       1 0 0 0    HABCDEFG
  // Shift 2       CG DH AE BF       1 1 0 0    GHABCDEF
  // Shift 3       BF CG DH AE       1 1 1 0    FGHACBDE
  // Shift 4       AE BF CG DH       1 1 1 1    EFGHABCD
  // Shift 5       DH AE BF CG       0 1 1 1    DEFGHABC
  // Shift 6       CG DH AE BF       0 0 1 1    CDEFGHAB
  // Shift 7       BF CG DH AE       0 0 0 1    BCDEFGHA
  int permute_dim =rhs._grid->_layout[dimension]>1 ;
  int permute_type=0;
  for(int d=0;d<dimension;d++)
    if (rhs._grid->_layout[d]>1 ) permute_type++;
        
  // loop over perp slices.
  // Threading considerations:
  //   Need to map thread_num to
  //
  //               x_min,x_max for Loop-A
  //               n_min,n_max for Loop-B
  //               b_min,b_max for Loop-C
  //  In a way that maximally load balances.
  //
  //  Optimal:
  //      There are rd*n_block*block items of work.
  //      These serialise as item "w"
  //      b=w%block; w=w/block
  //      n=w%nblock; x=w/nblock. Perhaps 20 cycles?
  //
  //  Logic:
  //      x_chunk = (rd+thread)/nthreads simply divide work across nodes.
  //
  //      rd=5 , threads = 8;
  //      0 1 2 3 4 5 6 7
  //      0 0 0 1 1 1 1 1
  for(int x=0;x<rd;x++){         // Loop A
    sx = (x-shift+ld)%rd;
    o  = x*rhs._grid->_ostride[dimension];
    so =sx*rhs._grid->_ostride[dimension];
    int permute_slice=0;
    if ( permute_dim ) {
      permute_slice = shift/rd;
      if ( x<shift%rd ) permute_slice = 1-permute_slice;
    }
#if 0
    if ( permute_slice ) {
                
      int internal=sizeof(vobj)/sizeof(vComplex);
      int num =rhs._grid->_slice_block[dimension]*internal;
                
      for(int n=0;n<rhs._grid->_slice_nblock[dimension];n++){
	vComplex *optr = (vComplex *)&ret._odata[o];
	vComplex *iptr = (vComplex *)&rhs._odata[so];
	for(int b=0;b<num;b++){
	  permute(optr[b],iptr[b],permute_type);
	}
	o+=rhs._grid->_slice_stride[dimension];
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
#else
    if ( permute_slice ) {
      int internal=sizeof(vobj)/sizeof(vComplex);
      int num =rhs._grid->_slice_block[dimension]*internal;
#pragma omp parallel for collapse(2)
      for(int n=0;n<rhs._grid->_slice_nblock[dimension];n++){
	for(int b=0;b<num;b++){
	  vComplex *optr = (vComplex *)&ret._odata[o +n*rhs._grid->_slice_stride[dimension]];
	  vComplex *iptr = (vComplex *)&rhs._odata[so+n*rhs._grid->_slice_stride[dimension]];
	  permute(optr[b],iptr[b],permute_type);
	}
      }
    } else {
#pragma omp parallel for collapse(2)
      for(int n=0;n<rhs._grid->_slice_nblock[dimension];n++){
	for(int i=0;i<rhs._grid->_slice_block[dimension];i++){
	  int oo = o+ n*rhs._grid->_slice_stride[dimension];
	  int soo=so+ n*rhs._grid->_slice_stride[dimension];
	  ret._odata[oo+i]=rhs._odata[soo+i];
	}
      }
            
    }
#endif
  }
  return ret;
}
friend Lattice<vobj> Cshift(Lattice<vobj> &rhs,int dimension,int shift)
{
  Lattice<vobj> ret(rhs._grid);
        
  ret.checkerboard = rhs._grid->CheckerBoardDestination(rhs.checkerboard,shift);
  shift = rhs._grid->CheckerBoardShift(rhs.checkerboard,dimension,shift);
  int rd = rhs._grid->_rdimensions[dimension];
  int ld = rhs._grid->_dimensions[dimension];
        
  // Map to always positive shift.
  shift = (shift+ld)%ld;
        
  // Work out whether to permute and the permute type
  // ABCDEFGH ->   AE BF CG DH       permute
  // Shift 0       AE BF CG DH       0 0 0 0    ABCDEFGH
  // Shift 1       DH AE BF CG       1 0 0 0    HABCDEFG
  // Shift 2       CG DH AE BF       1 1 0 0    GHABCDEF
  // Shift 3       BF CG DH AE       1 1 1 0    FGHACBDE
  // Shift 4       AE BF CG DH       1 1 1 1    EFGHABCD
  // Shift 5       DH AE BF CG       0 1 1 1    DEFGHABC
  // Shift 6       CG DH AE BF       0 0 1 1    CDEFGHAB
  // Shift 7       BF CG DH AE       0 0 0 1    BCDEFGHA
  int permute_dim =rhs._grid->_layout[dimension]>1 ;
  int permute_type=0;
  for(int d=0;d<dimension;d++)
    if (rhs._grid->_layout[d]>1 ) permute_type++;
        
        
  // loop over all work
  int work =rd*rhs._grid->_slice_nblock[dimension]*rhs._grid->_slice_block[dimension];

#pragma omp parallel for
  for(int ww=0;ww<work;ww++){

            
    // can optimise this if know w moves congtiguously for a given thread.
    // b=(b+1);
    // if (b==_slice_block) {b=0; n=n+1;}
    // if (n==_slice_nblock) { n=0; x=x+1}
    //
    // Perhaps a five cycle iterator, or so.
    int w=ww;
    int b = w%rhs._grid->_slice_block[dimension] ; w=w/rhs._grid->_slice_block[dimension];
    int n = w%rhs._grid->_slice_nblock[dimension]; w=w/rhs._grid->_slice_nblock[dimension];
    int x = w;

    int sx,so,o;
    sx = (x-shift+ld)%rd;
    o  = x*rhs._grid->_ostride[dimension]+n*rhs._grid->_slice_stride[dimension]; // common sub expression alert.
    so =sx*rhs._grid->_ostride[dimension]+n*rhs._grid->_slice_stride[dimension];
            
    int permute_slice=0;
    if ( permute_dim ) {
      permute_slice = shift/rd;
      if ( x<shift%rd ) permute_slice = 1-permute_slice;
    }
            
    if ( permute_slice ) {
                
      int internal=sizeof(vobj)/sizeof(vComplex);
      vComplex *optr = (vComplex *)&ret._odata[o+b];
      vComplex *iptr = (vComplex *)&rhs._odata[so+b];
      const char *pf = (const char *)iptr;
      for(int i=0;i<sizeof(vobj);i+=64){
	_mm_prefetch(pf+i,_MM_HINT_T0);
      }

      for(int i=0;i<internal;i++){
	permute(optr[i],iptr[i],permute_type);
      }
    } else {
      const char *pf = (const char *) &rhs._odata[so+b];
      for(int i=0;i<sizeof(vobj);i+=64){
	_mm_prefetch(pf+i,_MM_HINT_T0);
      }
      ret._odata[o+b]=rhs._odata[so+b];
    }
  }
  return ret;
}
#endif
#endif

