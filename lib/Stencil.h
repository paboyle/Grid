    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/Stencil.h

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>

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

#include <thread>

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

  template<class vobj,class cobj>
  class CartesianStencil { // Stencil runs along coordinate axes only; NO diagonal fill in.
  public:

      typedef uint32_t StencilInteger;
      typedef typename cobj::vector_type vector_type;
      typedef typename cobj::scalar_type scalar_type;
      typedef typename cobj::scalar_object scalar_object;

      //////////////////////////////////////////
      // Comms packet queue for asynch thread
      //////////////////////////////////////////

      struct Packet {
	void * send_buf;
	void * recv_buf;
	Integer to_rank;
	Integer from_rank;
	Integer bytes;
	volatile Integer done;
      };

      std::vector<Packet> Packets;

#define SEND_IMMEDIATE
#define SERIAL_SENDS

      void AddPacket(void *xmit,void * rcv, Integer to,Integer from,Integer bytes){
	comms_bytes+=2.0*bytes;
#ifdef SEND_IMMEDIATE
	commtime-=usecond();
	_grid->SendToRecvFrom(xmit,to,rcv,from,bytes);
	commtime+=usecond();
#endif
	Packet p;
	p.send_buf = xmit;
	p.recv_buf = rcv;
	p.to_rank  = to;
	p.from_rank= from;
	p.bytes    = bytes;
	p.done     = 0;
	comms_bytes+=2.0*bytes;
	Packets.push_back(p);

      }

#ifdef SERIAL_SENDS
      void Communicate(void ) { 
	commtime-=usecond();
	for(int i=0;i<Packets.size();i++){
#ifndef SEND_IMMEDIATE
	  _grid->SendToRecvFrom(
				Packets[i].send_buf,
				Packets[i].to_rank,
				Packets[i].recv_buf,
				Packets[i].from_rank,
				Packets[i].bytes);
#endif
	  Packets[i].done = 1;
	}
	commtime+=usecond();
      }
#else
      void Communicate(void ) { 
	typedef CartesianCommunicator::CommsRequest_t CommsRequest_t;
	std::vector<std::vector<CommsRequest_t> > reqs(Packets.size());
	commtime-=usecond();
	const int concurrency=2;
	for(int i=0;i<Packets.size();i+=concurrency){
	  for(int ii=0;ii<concurrency;ii++){
	    int j = i+ii;
	    if ( j<Packets.size() ) {
#ifndef SEND_IMMEDIATE
	      _grid->SendToRecvFromBegin(reqs[j],
					 Packets[j].send_buf,
					 Packets[j].to_rank,
					 Packets[j].recv_buf,
					 Packets[j].from_rank,
					 Packets[j].bytes);
#endif
	    }
	  }
	  for(int ii=0;ii<concurrency;ii++){
	    int j = i+ii;
	    if ( j<Packets.size() ) {
#ifndef SEND_IMMEDIATE
	      _grid->SendToRecvFromComplete(reqs[i]);
#endif
	    }
	  }
	  for(int ii=0;ii<concurrency;ii++){
	    int j = i+ii;
	    if ( j<Packets.size() ) {
	      Packets[j].done = 1;
	    }
	  }
	}
	commtime+=usecond();
      }
#endif

      ///////////////////////////////////////////
      // Simd merge queue for asynch comms
      ///////////////////////////////////////////
      struct Merge {
        cobj * mpointer;
	std::vector<scalar_object *> rpointers;
	Integer buffer_size;
	Integer packet_id;
      };

      std::vector<Merge> Mergers;

      void AddMerge(cobj *merge_p,std::vector<scalar_object *> &rpointers,Integer buffer_size,Integer packet_id) {
	Merge m;
	m.mpointer = merge_p;
	m.rpointers= rpointers;
	m.buffer_size = buffer_size;
	m.packet_id   = packet_id;
#ifdef SEND_IMMEDIATE
	mergetime-=usecond();
PARALLEL_FOR_LOOP
        for(int o=0;o<m.buffer_size;o++){
	  merge1(m.mpointer[o],m.rpointers,o);
	}
	mergetime+=usecond();
#else
	Mergers.push_back(m);
#endif

      }

      void CommsMerge(void ) { 
	//PARALLEL_NESTED_LOOP2 
	for(int i=0;i<Mergers.size();i++){	
	  
	spintime-=usecond();
	int packet_id = Mergers[i].packet_id;
	while(! Packets[packet_id].done ); // spin for completion
	spintime+=usecond();

#ifndef SEND_IMMEDIATE
	mergetime-=usecond();
PARALLEL_FOR_LOOP
	  for(int o=0;o<Mergers[i].buffer_size;o++){
	    merge1(Mergers[i].mpointer[o],Mergers[i].rpointers,o);
	  }
	mergetime+=usecond();
#endif

	}
      }

      ////////////////////////////////////////
      // Basic Grid and stencil info
      ////////////////////////////////////////

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
      inline StencilEntry * GetEntry(int &ptype,int point,int osite) { ptype = _permute_type[point]; return & _entries[point][osite]; }

      // Comms buffers
      std::vector<Vector<scalar_object> > u_simd_send_buf;
      std::vector<Vector<scalar_object> > u_simd_recv_buf;
      Vector<cobj>          u_send_buf;
      Vector<cobj>          comm_buf;
      int u_comm_offset;
      int _unified_buffer_size;

      /////////////////////////////////////////
      // Timing info; ugly; possibly temporary
      /////////////////////////////////////////
#define TIMING_HACK
#ifdef TIMING_HACK
      double jointime;
      double gathertime;
      double commtime;
      double halogtime;
      double mergetime;
      double spintime;
      double comms_bytes;
      double gathermtime;
      double splicetime;
      double nosplicetime;
#endif

  CartesianStencil(GridBase *grid,
				     int npoints,
				     int checkerboard,
				     const std::vector<int> &directions,
				     const std::vector<int> &distances) 
    :   _entries(npoints), _permute_type(npoints), _comm_buf_size(npoints)
    {
#ifdef TIMING_HACK
      gathertime=0;
      jointime=0;
      commtime=0;
      halogtime=0;
      mergetime=0;
      spintime=0;
      gathermtime=0;
      splicetime=0;
      nosplicetime=0;
      comms_bytes=0;
#endif
      _npoints = npoints;
      _grid    = grid;
      _directions = directions;
      _distances  = distances;
      _unified_buffer_size=0;

      int osites  = _grid->oSites();

      for(int ii=0;ii<npoints;ii++){

	int i = ii; // reverse direction to get SIMD comms done first
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
	int rotate_dim      = _grid->_simd_layout[dimension]>2;

	assert ( (rotate_dim && comm_dim) == false) ; // Do not think spread out is supported

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
      u_send_buf.resize(_unified_buffer_size);
      comm_buf.resize(_unified_buffer_size);
      
      const int Nsimd = grid->Nsimd();
      u_simd_send_buf.resize(Nsimd);
      u_simd_recv_buf.resize(Nsimd);
      for(int l=0;l<Nsimd;l++){
	u_simd_send_buf[l].resize(_unified_buffer_size);
	u_simd_recv_buf[l].resize(_unified_buffer_size);
      }
    }

    void Local     (int point, int dimension,int shiftpm,int cbmask)
    {
      int fd = _grid->_fdimensions[dimension];
      int rd = _grid->_rdimensions[dimension];
      int ld = _grid->_ldimensions[dimension];
      int gd = _grid->_gdimensions[dimension];
      int ly = _grid->_simd_layout[dimension];

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
	  else permute_slice = (wrap+1)%ly;
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
      
      assert(comm_dim==1);
      int shift = (shiftpm + fd) %fd;
      assert(shift>=0);
      assert(shift<fd);

      int buffer_size = _grid->_slice_nblock[dimension]*_grid->_slice_block[dimension]; // done in reduced dims, so SIMD factored

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
	}


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


      template<class compressor>
      std::thread HaloExchangeBegin(const Lattice<vobj> &source,compressor &compress) {
	Mergers.resize(0); 
	Packets.resize(0);
	HaloGather(source,compress);
        return std::thread([&] { this->Communicate(); });
      }

      template<class compressor>
      void HaloExchange(const Lattice<vobj> &source,compressor &compress) 
      {
	Mergers.resize(0); 
	Packets.resize(0);
	HaloGather(source,compress);
	Communicate();
	CommsMerge();
      }

      void HaloExchangeComplete(std::thread &thr) 
      {
	CommsMerge(); // spins
	jointime-=usecond();
	thr.join();
	jointime+=usecond();
      }

      template<class compressor>
      void HaloGatherDir(const Lattice<vobj> &source,compressor &compress,int point)
      {
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
		GatherSimd(source,dimension,shift,0x3,compress);
		splicetime+=usecond();
	      } else { 
		nosplicetime-=usecond();
		Gather(source,dimension,shift,0x3,compress);
		nosplicetime+=usecond();
	      }
	    } else {
	      if(splice_dim){
		splicetime-=usecond();
		GatherSimd(source,dimension,shift,0x1,compress);// if checkerboard is unfavourable take two passes
		GatherSimd(source,dimension,shift,0x2,compress);// both with block stride loop iteration
		splicetime+=usecond();
	      } else {
		nosplicetime-=usecond();
		Gather(source,dimension,shift,0x1,compress);
		Gather(source,dimension,shift,0x2,compress);
		nosplicetime+=usecond();
	      }
	    }
	  }
      }

      template<class compressor>
      void HaloGather(const Lattice<vobj> &source,compressor &compress)
      {
	// conformable(source._grid,_grid);
	assert(source._grid==_grid);
	halogtime-=usecond();

	assert (comm_buf.size() == _unified_buffer_size );
	u_comm_offset=0;

	// Gather all comms buffers
	for(int point = 0 ; point < _npoints; point++) {
	  compress.Point(point);
	  HaloGatherDir(source,compress,point);
	}

	assert(u_comm_offset==_unified_buffer_size);
	halogtime+=usecond();
      }

      template<class compressor>
        void Gather(const Lattice<vobj> &rhs,int dimension,int shift,int cbmask,compressor & compress)
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
	      Gather_plane_simple (rhs,u_send_buf,dimension,sx,cbmask,compress,u_comm_offset);
	      gathertime+=usecond();

	      int rank           = _grid->_processor;
	      int recv_from_rank;
	      int xmit_to_rank;
	      _grid->ShiftedRanks(dimension,comm_proc,xmit_to_rank,recv_from_rank);
	      assert (xmit_to_rank   != _grid->ThisRank());
	      assert (recv_from_rank != _grid->ThisRank());

	      //      FIXME Implement asynchronous send & also avoid buffer copy
	      AddPacket((void *)&u_send_buf[u_comm_offset],
			(void *)  &comm_buf[u_comm_offset],
			xmit_to_rank,
			recv_from_rank,
			bytes);
			
	      u_comm_offset+=words;
	    }
	  }
	}


      template<class compressor>
	void  GatherSimd(const Lattice<vobj> &rhs,int dimension,int shift,int cbmask,compressor &compress)
	{
	  const int Nsimd = _grid->Nsimd();
	  
	  int fd = _grid->_fdimensions[dimension];
	  int rd = _grid->_rdimensions[dimension];
	  int ld = _grid->_ldimensions[dimension];
	  int pd              = _grid->_processors[dimension];
	  int simd_layout     = _grid->_simd_layout[dimension];
	  int comm_dim        = _grid->_processors[dimension] >1 ;

	  assert(comm_dim==1);
	  // This will not work with a rotate dim
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

	  int bytes = buffer_size*sizeof(scalar_object);

	  std::vector<scalar_object *> rpointers(Nsimd);
	  std::vector<scalar_object *> spointers(Nsimd);

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
		spointers[i] = &u_simd_send_buf[i][u_comm_offset];
	      }

	      int sx   = (x+sshift)%rd;
	      
	      gathermtime-=usecond();
	      Gather_plane_extract<cobj>(rhs,spointers,dimension,sx,cbmask,compress);
	      gathermtime+=usecond();

	      for(int i=0;i<Nsimd;i++){
		
		// FIXME 
		// This logic is hard coded to simd_layout ==2 and not allowing >2
		//		std::cout << "GatherSimd : lane 1st elem " << i << u_simd_send_buf[i ][u_comm_offset]<<std::endl;

		int inner_bit = (Nsimd>>(permute_type+1));
		int ic= (i&inner_bit)? 1:0;

		int my_coor          = rd*ic + x;
		int nbr_coor         = my_coor+sshift;
		int nbr_proc = ((nbr_coor)/ld) % pd;// relative shift in processors
		int nbr_lcoor= (nbr_coor%ld);
		int nbr_ic   = (nbr_lcoor)/rd;    // inner coord of peer
		int nbr_ox   = (nbr_lcoor%rd);    // outer coord of peer
		int nbr_lane = (i&(~inner_bit));
		
		if (nbr_ic) nbr_lane|=inner_bit;
		assert (sx == nbr_ox);

		auto rp = &u_simd_recv_buf[i       ][u_comm_offset];
		auto sp = &u_simd_send_buf[nbr_lane][u_comm_offset];

		void *vrp = (void *)rp;
		void *vsp = (void *)sp;


		if(nbr_proc){
		  
		  int recv_from_rank;
		  int xmit_to_rank;

		  _grid->ShiftedRanks(dimension,nbr_proc,xmit_to_rank,recv_from_rank); 
		  
		  AddPacket( vsp,vrp,xmit_to_rank,recv_from_rank,bytes);
		  
		  rpointers[i] = rp;

		} else { 

		  rpointers[i] = sp;

		}
	      }

	      AddMerge(&comm_buf[u_comm_offset],rpointers,buffer_size,Packets.size()-1);

	      u_comm_offset     +=buffer_size;
	    }
	  }
	}

  };
}
#endif
