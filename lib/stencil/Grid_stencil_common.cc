#include "Grid.h"

namespace Grid {



void CartesianStencil::LebesgueOrder(void) 
{
  _LebesgueReorder.resize(0);
  
  // Align up dimensions to power of two.
  const StencilInteger one=1;
  StencilInteger ND = _grid->_ndimension;
  std::vector<StencilInteger> dims(ND);
  std::vector<StencilInteger> adims(ND);
  std::vector<std::vector<StencilInteger> > bitlist(ND);


  for(StencilInteger mu=0;mu<ND;mu++){
    dims[mu] = _grid->_rdimensions[mu];
    assert ( dims[mu] != 0 );
    adims[mu] = alignup(dims[mu]);
  }

  // List which bits of padded volume coordinate contribute; this strategy 
  // i) avoids recursion 
  // ii) has loop lengths at most the width of a 32 bit word.
  int sitebit=0;
  int split=24;
  for(int mu=0;mu<ND;mu++){   // mu 0 takes bit 0; mu 1 takes bit 1 etc...
    for(int bit=0;bit<split;bit++){
    StencilInteger mask = one<<bit;
      if ( mask&(adims[mu]-1) ){
	bitlist[mu].push_back(sitebit);
	sitebit++;
      }
    }
  }
  for(int bit=split;bit<32;bit++){
    StencilInteger mask = one<<bit;
    for(int mu=0;mu<ND;mu++){   // mu 0 takes bit 0; mu 1 takes bit 1 etc...
      if ( mask&(adims[mu]-1) ){
	bitlist[mu].push_back(sitebit);
	sitebit++;
      }
    }
  }

  // Work out padded and unpadded volumes
  StencilInteger avol = 1;
  for(int mu=0;mu<ND;mu++) avol = avol * adims[mu];

  StencilInteger vol = 1;
  for(int mu=0;mu<ND;mu++) vol = vol * dims[mu];
  
  // Loop over padded volume, following Lebesgue curve
  // We interleave the bits from sequential "mu".
  std::vector<StencilInteger> ax(ND);
  
  for(StencilInteger asite=0;asite<avol;asite++){

    // Start with zero and collect bits
    for(int mu=0;mu<ND;mu++) ax[mu] = 0;

    int contained = 1;
    for(int mu=0;mu<ND;mu++){

      // Build the coordinate on the aligned volume
      for(int bit=0;bit<bitlist[mu].size();bit++){
	int sbit=bitlist[mu][bit];

	if(asite&(one<<sbit)){
	  ax[mu]|=one<<bit;
	}
      }

      // Is it contained in original box
      if ( ax[mu]>dims[mu]-1 ) contained = 0;

    }

    if ( contained ) {
      int site = ax[0]
	+        dims[0]*ax[1]
        +dims[0]*dims[1]*ax[2]
        +dims[0]*dims[1]*dims[2]*ax[3];

      _LebesgueReorder.push_back(site);
    }
  }

  assert( _LebesgueReorder.size() == vol );
}

  CartesianStencil::CartesianStencil(GridBase *grid,
				     int npoints,
				     int checkerboard,
				     const std::vector<int> &directions,
				     const std::vector<int> &distances) 
    :   _offsets(npoints), 
	_is_local(npoints), 
	_comm_buf_size(npoints), 
	_permute_type(npoints),
	_permute(npoints)
    {
      _npoints = npoints;
      _grid    = grid;
      _directions = directions;
      _distances  = distances;
      _unified_buffer_size=0;
      _request_count =0;

      LebesgueOrder();

      int osites  = _grid->oSites();

      for(int i=0;i<npoints;i++){

	int point = i;

	_offsets[i].resize( osites);
	_is_local[i].resize(osites);
	_permute[i].resize( osites);

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
	  sshift[0] = _grid->CheckerBoardShift(_checkerboard,dimension,shift,0);
	  sshift[1] = _grid->CheckerBoardShift(_checkerboard,dimension,shift,1);

	  if ( sshift[0] == sshift[1] ) {
	    Local(point,dimension,shift,0x3);
	  } else {
	    Local(point,dimension,shift,0x1);// if checkerboard is unfavourable take two passes
	    Local(point,dimension,shift,0x2);// both with block stride loop iteration
	  }
	} else { // All permute extract done in comms phase prior to Stencil application
	  //        So tables are the same whether comm_dim or splice_dim
	  sshift[0] = _grid->CheckerBoardShift(_checkerboard,dimension,shift,0);
	  sshift[1] = _grid->CheckerBoardShift(_checkerboard,dimension,shift,1);
	  if ( sshift[0] == sshift[1] ) {
	    Comms(point,dimension,shift,0x3);
	  } else {
	    Comms(point,dimension,shift,0x1);// if checkerboard is unfavourable take two passes
	    Comms(point,dimension,shift,0x2);// both with block stride loop iteration
	  }
	}
      }
    }


    void CartesianStencil::Local     (int point, int dimension,int shift,int cbmask)
    {
      int fd = _grid->_fdimensions[dimension];
      int rd = _grid->_rdimensions[dimension];
      int ld = _grid->_ldimensions[dimension];
      int gd = _grid->_gdimensions[dimension];
      
      // Map to always positive shift modulo global full dimension.
      shift = (shift+fd)%fd;
      
      // the permute type
      int permute_dim =_grid->PermuteDim(dimension);
      
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

  	CopyPlane(point,dimension,x,sx,cbmask,permute_slice);
  
      }
    }

    void CartesianStencil::Comms     (int point,int dimension,int shift,int cbmask)
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
      shift = (shift + fd) %fd;
      assert(shift>=0);
      assert(shift<fd);
      
      int buffer_size = _grid->_slice_nblock[dimension]*_grid->_slice_block[dimension];
      _comm_buf_size[point] = buffer_size; // Size of _one_ plane. Multiple planes may be gathered and
                                           // send to one or more remote nodes.

      int cb= (cbmask==0x2)? 1 : 0;
      int sshift= _grid->CheckerBoardShift(_checkerboard,dimension,shift,cb);
      
      for(int x=0;x<rd;x++){       
	
	int offnode = (((x+sshift)%fd) >= rd ); 
	//	int comm_proc   = ((x+sshift)/ld)%pd;        
	//	int offnode     = (comm_proc!=0);
	int sx          = (x+sshift)%rd;

	if (!offnode) {
	  
	  int permute_slice=0;
	  CopyPlane(point,dimension,x,sx,cbmask,permute_slice); 
	  
	} else {
	  
	  int words = buffer_size;
	  if (cbmask != 0x3) words=words>>1;
	  
	  //	  GatherPlaneSimple (point,dimension,sx,cbmask);
	  
	  int rank           = grid->_processor;
	  int recv_from_rank;
	  int xmit_to_rank;

	  int unified_buffer_offset = _unified_buffer_size;
	  _unified_buffer_size    += words;
	  ScatterPlane(point,dimension,x,cbmask,unified_buffer_offset); // permute/extract/merge is done in comms phase
	  
	}
      }
    }
  // Routine builds up integer table for each site in _offsets, _is_local, _permute
  void CartesianStencil::CopyPlane(int point, int dimension,int lplane,int rplane,int cbmask,int permute)
    {
      int rd = _grid->_rdimensions[dimension];
      
      if ( !_grid->CheckerBoarded(dimension) ) {
	
	int o   = 0;                                     // relative offset to base within plane
	int ro  = rplane*_grid->_ostride[dimension]; // base offset for start of plane 
	int lo  = lplane*_grid->_ostride[dimension]; // offset in buffer
	
	// Simple block stride gather of SIMD objects
	for(int n=0;n<_grid->_slice_nblock[dimension];n++){
	  for(int b=0;b<_grid->_slice_block[dimension];b++){
	    _offsets [point][lo+o+b]=ro+o+b;
	    _is_local[point][lo+o+b]=1;
	    _permute [point][lo+o+b]=permute;
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
	      _offsets [point][lo+o+b]=ro+o+b;
	      _is_local[point][lo+o+b]=1;
	      _permute [point][lo+o+b]=permute;
	    }
	    
	    }
	  o +=_grid->_slice_stride[dimension];
	}
	
      }
    }
  // Routine builds up integer table for each site in _offsets, _is_local, _permute
    void CartesianStencil::ScatterPlane (int point,int dimension,int plane,int cbmask,int offset)
    {
      int rd = _grid->_rdimensions[dimension];
      
      if ( !_grid->CheckerBoarded(dimension) ) {
	
	int so  = plane*_grid->_ostride[dimension]; // base offset for start of plane 
	int o   = 0;                                    // relative offset to base within plane
	int bo  = 0;                                    // offset in buffer

	// Simple block stride gather of SIMD objects
	for(int n=0;n<_grid->_slice_nblock[dimension];n++){
	  for(int b=0;b<_grid->_slice_block[dimension];b++){
	    _offsets [point][so+o+b]=offset+(bo++);
	    _is_local[point][so+o+b]=0;
	    _permute [point][so+o+b]=0;
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
	      _offsets [point][so+o+b]=offset+(bo++);
	      _is_local[point][so+o+b]=0;
	      _permute [point][so+o+b]=0;
	    }
	  }
	  o +=_grid->_slice_stride[dimension];
	}
      }
    }
}
