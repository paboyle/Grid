    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/Stencil.h

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: Peter Boyle <peterboyle@Peters-MacBook-Pro-2.local>
Author: paboyle <paboyle@ph.ed.ac.uk>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

    See the full license in the file "LICENSE" in the top level distribution directory
    *************************************************************************************/
    /*  END LEGAL */
#ifndef GRID_STENCIL_H
#define GRID_STENCIL_H

#include <stencil/Lebesgue.h>   // subdir aggregate

//////////////////////////////////////////////////////////////////////////////////////////
// Must not lose sight that goal is to be able to construct really efficient
// gather to a point stencil code. CSHIFT is not the best way, so need
// additional stencil support.
//
// Stencil based code will pre-exchange haloes and use a table lookup for neighbours.
// This will be done with generality to allow easier efficient implementations.
// Overlap of comms and compute could be semi-automated by tabulating off-node connected,
// and 
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
  
  struct StencilEntry { 
    int _offset;
    int _is_local;
    int _permute;
    int _around_the_world;
  };

  template<class vobj,class cobj, class compressor>
  class CartesianStencil { // Stencil runs along coordinate axes only; NO diagonal fill in.
  public:

      typedef uint32_t StencilInteger;
      typedef typename cobj::vector_type vector_type;
      typedef typename cobj::scalar_type scalar_type;
      typedef typename cobj::scalar_object scalar_object;

      int                               _checkerboard;
      int                               _npoints; // Move to template param?
      GridBase *                        _grid;
      
      // npoints of these
      std::vector<int>                  _directions;
      std::vector<int>                  _distances;
      std::vector<int>                  _comm_buf_size;
      std::vector<int>                  _permute_type;

      // npoints x Osites() of these
      std::vector<std::vector<StencilEntry> > _entries;

      // Comms buffers
      std::vector<std::vector<scalar_object> > send_buf_extract;
      std::vector<std::vector<scalar_object> > recv_buf_extract;
      std::vector<scalar_object *> pointers;
      std::vector<scalar_object *> rpointers;
      Vector<cobj> send_buf;

      inline StencilEntry * GetEntry(int &ptype,int point,int osite) { ptype = _permute_type[point]; return & _entries[point][osite]; }

      int _unified_buffer_size;
      int _request_count;

      double buftime;
      double gathertime;
      double commtime;
      double commstime;
      double halotime;
      double scattertime;
      double mergetime;
      double gathermtime;
      double splicetime;
      double nosplicetime;




  CartesianStencil(GridBase *grid,
				     int npoints,
				     int checkerboard,
				     const std::vector<int> &directions,
				     const std::vector<int> &distances) 
    :   _entries(npoints), _permute_type(npoints), _comm_buf_size(npoints)
    {
      gathertime=0;
      commtime=0;
      commstime=0;
      halotime=0;
      scattertime=0;
      mergetime=0;
      gathermtime=0;
      buftime=0;
      splicetime=0;
      nosplicetime=0;

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
	    //	    std::cout<<"Comms 0x3"<<std::endl;
	    Comms(point,dimension,shift,0x3);
	  } else {
	    //	    std::cout<<"Comms 0x1 ; 0x2"<<std::endl;
	    Comms(point,dimension,shift,0x1);// if checkerboard is unfavourable take two passes
	    Comms(point,dimension,shift,0x2);// both with block stride loop iteration
	  }
	}
	//	for(int ss=0;ss<osites;ss++){
	//	  std::cout << "point["<<i<<"] "<<ss<<"-> o"<<_entries[i][ss]._offset<<"; l"<<
	//	    _entries[i][ss]._is_local<<"; p"<<_entries[i][ss]._permute<<std::endl;
	//	}
      }
    }


    void Local     (int point, int dimension,int shiftpm,int cbmask)
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

    void Comms     (int point,int dimension,int shiftpm,int cbmask)
    {
      GridBase *grid=_grid;
      const int Nsimd = grid->Nsimd();
      
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

      int buffer_size = _grid->_slice_nblock[dimension]*_grid->_slice_block[dimension]; // done in reduced dims, so SIMD factored
      //      std::cout << " dim " <<dimension<<" buffersize "<<buffer_size<<std::endl;
      _comm_buf_size[point] = buffer_size; // Size of _one_ plane. Multiple planes may be gathered and
                                           // send to one or more remote nodes.

      int cb= (cbmask==0x2)? Odd : Even;
      int sshift= _grid->CheckerBoardShiftForCB(_checkerboard,dimension,shift,cb);
      

      for(int x=0;x<rd;x++){       


	int permute_type=grid->PermuteType(dimension);

	int sx        =  (x+sshift)%rd;
	
    	int offnode = 0;
	if ( simd_layout > 1 ) {

	  for(int i=0;i<Nsimd;i++){

	    int inner_bit = (Nsimd>>(permute_type+1));
	    int ic= (i&inner_bit)? 1:0;
	    int my_coor          = rd*ic + x;
	    int nbr_coor         = my_coor+sshift;
	    int nbr_proc = ((nbr_coor)/ld) % pd;// relative shift in processors

	    if ( nbr_proc ) { 
	      offnode =1;
	    }
	  }
	  
	} else { 
	  int comm_proc = ((x+sshift)/rd)%pd;
	  offnode = (comm_proc!= 0);
	  //	  std::cout << "Stencil x "<<x<<" shift "<<shift<<" sshift "<<sshift<<" fd "<<fd<<" rd " <<rd<<" offnode "<<offnode<<" sx "<<sx<< " comm_proc "<<comm_proc<<" pd "<< pd <<std::endl;
	}


	// Stencil x 1 shift 3 sshift 3 fd 8 rd 2 offnode 0 sx 0 comm_proc 0 pd 2
	// x+sshift = 4
	// x+sshift/2 = 2
	// 2%2 == 0
	// Problem: sshift is wrong in "rd" for SIMD directions. The complex logic in Cshift_mpi is needed.

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
	  //	  std::cout<< "Comms dim "<<dimension<<" offset "<<unified_buffer_offset<<" size "<<" " << _unified_buffer_size<<std::endl;
	  ScatterPlane(point,dimension,x,cbmask,unified_buffer_offset,wraparound); // permute/extract/merge is done in comms phase
	  
	}
      }
    }
  // Routine builds up integer table for each site in _offsets, _is_local, _permute
  void CopyPlane(int point, int dimension,int lplane,int rplane,int cbmask,int permute,int wrap)
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
   void ScatterPlane (int point,int dimension,int plane,int cbmask,int offset, int wrap)
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

//      CartesianStencil(GridBase *grid,
//		       int npoints,
//		       int checkerboard,
//		       const std::vector<int> &directions,
//		       const std::vector<int> &distances);


      // Add to tables for various cases;  is this mistaken. only local if 1 proc in dim
      // Can this be avoided with simpler coding of comms?
   //      void Local     (int point, int dimension,int shift,int cbmask);
   //      void Comms     (int point, int dimension,int shift,int cbmask);
   //      void CopyPlane(int point, int dimension,int lplane,int rplane,int cbmask,int permute,int wrap);
   //      void ScatterPlane (int point,int dimension,int plane,int cbmask,int offset,int wrap);

      // Could allow a functional munging of the halo to another type during the comms.
      // this could implement the 16bit/32bit/64bit compression.
      void HaloExchange(const Lattice<vobj> &source,std::vector<cobj,alignedAllocator<cobj> > &u_comm_buf,compressor &compress)
      {
	// conformable(source._grid,_grid);
	assert(source._grid==_grid);
	halotime-=usecond();
	if (u_comm_buf.size() != _unified_buffer_size ) u_comm_buf.resize(_unified_buffer_size);
	int u_comm_offset=0;

	// Gather all comms buffers
	for(int point = 0 ; point < _npoints; point++) {

	  compress.Point(point);

	  int dimension    = _directions[point];
	  int displacement = _distances[point];
	  
	  int fd = _grid->_fdimensions[dimension];
	  int rd = _grid->_rdimensions[dimension];


	  // Map to always positive shift modulo global full dimension.
	  int shift = (displacement+fd)%fd;
	  
	  //     	  int checkerboard = _grid->CheckerBoardDestination(source.checkerboard,shift);
	  assert (source.checkerboard== _checkerboard);

	  // the permute type
	  int simd_layout     = _grid->_simd_layout[dimension];
	  int comm_dim        = _grid->_processors[dimension] >1 ;
	  int splice_dim      = _grid->_simd_layout[dimension]>1 && (comm_dim);

	  // Gather phase
	  int sshift [2];
	  if ( comm_dim ) {
	    sshift[0] = _grid->CheckerBoardShiftForCB(_checkerboard,dimension,shift,Even);
	    sshift[1] = _grid->CheckerBoardShiftForCB(_checkerboard,dimension,shift,Odd);
	    if ( sshift[0] == sshift[1] ) {
	      if (splice_dim) {
		splicetime-=usecond();
		GatherStartCommsSimd(source,dimension,shift,0x3,u_comm_buf,u_comm_offset,compress);
		splicetime+=usecond();
	      } else { 
		nosplicetime-=usecond();
		GatherStartComms(source,dimension,shift,0x3,u_comm_buf,u_comm_offset,compress);
		nosplicetime+=usecond();
	      }
	    } else {
	      //	      std::cout << "dim "<<dimension<<"cb "<<_checkerboard<<"shift "<<shift<<" sshift " << sshift[0]<<" "<<sshift[1]<<std::endl;
	      if(splice_dim){
		splicetime-=usecond();
		GatherStartCommsSimd(source,dimension,shift,0x1,u_comm_buf,u_comm_offset,compress);// if checkerboard is unfavourable take two passes
		GatherStartCommsSimd(source,dimension,shift,0x2,u_comm_buf,u_comm_offset,compress);// both with block stride loop iteration
		splicetime+=usecond();
	      } else {
		nosplicetime-=usecond();
		GatherStartComms(source,dimension,shift,0x1,u_comm_buf,u_comm_offset,compress);
		GatherStartComms(source,dimension,shift,0x2,u_comm_buf,u_comm_offset,compress);
		nosplicetime+=usecond();
	      }
	    }
	  }
	}
	halotime+=usecond();
      }

        void GatherStartComms(const Lattice<vobj> &rhs,int dimension,int shift,int cbmask,
			      std::vector<cobj,alignedAllocator<cobj> > &u_comm_buf,
			      int &u_comm_offset,compressor & compress)
	{
	  typedef typename cobj::vector_type vector_type;
	  typedef typename cobj::scalar_type scalar_type;
	  
	  GridBase *grid=_grid;
	  assert(rhs._grid==_grid);
	  //	  conformable(_grid,rhs._grid);

	  int fd              = _grid->_fdimensions[dimension];
	  int rd              = _grid->_rdimensions[dimension];
	  int pd              = _grid->_processors[dimension];
	  int simd_layout     = _grid->_simd_layout[dimension];
	  int comm_dim        = _grid->_processors[dimension] >1 ;
	  assert(simd_layout==1);
	  assert(comm_dim==1);
	  assert(shift>=0);
	  assert(shift<fd);

	  int buffer_size = _grid->_slice_nblock[dimension]*_grid->_slice_block[dimension];

	  if(send_buf.size()<buffer_size) send_buf.resize(buffer_size);

	  int cb= (cbmask==0x2)? Odd : Even;
	  int sshift= _grid->CheckerBoardShiftForCB(rhs.checkerboard,dimension,shift,cb);

	  for(int x=0;x<rd;x++){       

	    int sx        = (x+sshift)%rd;
	    int comm_proc = ((x+sshift)/rd)%pd;

	    if (comm_proc) {
	      
	      int words = buffer_size;
	      if (cbmask != 0x3) words=words>>1;
	    
	      int bytes = words * sizeof(cobj);

	      gathertime-=usecond();
	      Gather_plane_simple (rhs,send_buf,dimension,sx,cbmask,compress);
	      gathertime+=usecond();

	      int rank           = _grid->_processor;
	      int recv_from_rank;
	      int xmit_to_rank;
	      _grid->ShiftedRanks(dimension,comm_proc,xmit_to_rank,recv_from_rank);
	      assert (xmit_to_rank != _grid->ThisRank());
	      assert (recv_from_rank != _grid->ThisRank());

	      //      FIXME Implement asynchronous send & also avoid buffer copy
	      commtime-=usecond();
	      _grid->SendToRecvFrom((void *)&send_buf[0],
				   xmit_to_rank,
				    (void *)&u_comm_buf[u_comm_offset],
				   recv_from_rank,
				   bytes);
	      commtime+=usecond();

	      u_comm_offset+=words;
	    }
	  }
	}


	void  GatherStartCommsSimd(const Lattice<vobj> &rhs,int dimension,int shift,int cbmask,
				   std::vector<cobj,alignedAllocator<cobj> > &u_comm_buf,
				   int &u_comm_offset,compressor &compress)
	{
	  buftime-=usecond();
	  const int Nsimd = _grid->Nsimd();

	  
	  int fd = _grid->_fdimensions[dimension];
	  int rd = _grid->_rdimensions[dimension];
	  int ld = _grid->_ldimensions[dimension];
	  int pd              = _grid->_processors[dimension];
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
	  int buffer_size = _grid->_slice_nblock[dimension]*_grid->_slice_block[dimension];
	  int words = sizeof(cobj)/sizeof(vector_type);

	  assert(cbmask==0x3); // Fixme think there is a latent bug if not true

	  //	Should grow to max size and then cost very little thereafter
	  send_buf_extract.resize(Nsimd);
	  recv_buf_extract.resize(Nsimd);
	  for(int l=0;l<Nsimd;l++){
	    if( send_buf_extract[l].size() < buffer_size) {
	      send_buf_extract[l].resize(buffer_size);
	      recv_buf_extract[l].resize(buffer_size);
	    }
	  }
	  pointers.resize(Nsimd);
	  rpointers.resize(Nsimd);

	  int bytes = buffer_size*sizeof(scalar_object);
	  
	  buftime+=usecond();
	  
	  ///////////////////////////////////////////
	  // Work out what to send where
	  ///////////////////////////////////////////

	  int cb    = (cbmask==0x2)? Odd : Even;
	  int sshift= _grid->CheckerBoardShiftForCB(rhs.checkerboard,dimension,shift,cb);
	  
	  // loop over outer coord planes orthog to dim
	  for(int x=0;x<rd;x++){       

	    int any_offnode = ( ((x+sshift)%fd) >= rd );

	    if ( any_offnode ) {

	      for(int i=0;i<Nsimd;i++){       
		pointers[i] = &send_buf_extract[i][0];
	      }
	      int sx   = (x+sshift)%rd;
	      
	      gathermtime-=usecond();
	      Gather_plane_extract<cobj>(rhs,pointers,dimension,sx,cbmask,compress);
	      gathermtime+=usecond();

	      for(int i=0;i<Nsimd;i++){

		int inner_bit = (Nsimd>>(permute_type+1));
		int ic= (i&inner_bit)? 1:0;

		int my_coor          = rd*ic + x;
		int nbr_coor         = my_coor+sshift;
		int nbr_proc = ((nbr_coor)/ld) % pd;// relative shift in processors
		int nbr_lcoor= (nbr_coor%ld);
		int nbr_ic   = (nbr_lcoor)/rd;    // inner coord of peer
		int nbr_ox   = (nbr_lcoor%rd);    // outer coord of peer
		int nbr_lane = (i&(~inner_bit));
		
		int recv_from_rank;
		int xmit_to_rank;
		
		if (nbr_ic) nbr_lane|=inner_bit;
		assert (sx == nbr_ox);

		
		if(nbr_proc){
		  
		  _grid->ShiftedRanks(dimension,nbr_proc,xmit_to_rank,recv_from_rank); 
		  
		  commstime-=usecond();
		  _grid->SendToRecvFrom((void *)&send_buf_extract[nbr_lane][0],
					xmit_to_rank,
					(void *)&recv_buf_extract[i][0],
					recv_from_rank,
					bytes);
		  commstime+=usecond();
		  
		  rpointers[i] = &recv_buf_extract[i][0];

		} else { 
		  rpointers[i] = &send_buf_extract[nbr_lane][0];
		}
	      }

	      //	      std::cout << " CommsSimd ["<<dimension<<"] offset "<<u_comm_offset<<" buffsize "<<buffer_size  <<" unified  buffer size "<<_unified_buffer_size<<std::endl;
	      mergetime-=usecond();
PARALLEL_FOR_LOOP
	      for(int i=0;i<buffer_size;i++){
		//		std::cout<<"buffer loop " << i<<" "<<u_comm_offset+i<<" / "<<_unified_buffer_size<<std::endl;
		//		assert(u_comm_offset+i<_unified_buffer_size);
		merge(u_comm_buf[u_comm_offset+i],rpointers,i);
	      }
	      mergetime+=usecond();
	      u_comm_offset+=buffer_size;
	    }
	  }
	}
  };
}
#endif
