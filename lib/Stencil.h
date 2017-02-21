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

 #include <Grid/stencil/Lebesgue.h>   // subdir aggregate

#if defined __GNUC__
 #pragma GCC diagnostic push
 #pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#endif

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

inline void Gather_plane_simple_table_compute (GridBase *grid,int dimension,int plane,int cbmask,
					       int off,std::vector<std::pair<int,int> > & table)
{
  table.resize(0);
  //int rd = grid->_rdimensions[dimension];

  if ( !grid->CheckerBoarded(dimension) ) {
    cbmask = 0x3;
  }
  //int so= plane*grid->_ostride[dimension]; // base offset for start of plane 
  int e1=grid->_slice_nblock[dimension];
  int e2=grid->_slice_block[dimension];

  int stride=grid->_slice_stride[dimension];
  if ( cbmask == 0x3 ) { 
    table.resize(e1*e2);
    for(int n=0;n<e1;n++){
      for(int b=0;b<e2;b++){
	int o  = n*stride;
	int bo = n*e2;
	table[bo+b]=std::pair<int,int>(bo+b,o+b);
      }
    }
  } else { 
     int bo=0;
     table.resize(e1*e2/2);
     for(int n=0;n<e1;n++){
       for(int b=0;b<e2;b++){
	 int o  = n*stride;
	 int ocb=1<<grid->CheckerBoardFromOindexTable(o+b);
	 if ( ocb &cbmask ) {
	   table[bo]=std::pair<int,int>(bo,o+b); bo++;
	 }
       }
     }
  }
}

template<class vobj,class cobj,class compressor> void 
Gather_plane_simple_table (std::vector<std::pair<int,int> >& table,const Lattice<vobj> &rhs,cobj *buffer,compressor &compress, int off,int so)
{
PARALLEL_FOR_LOOP     
     for(int i=0;i<table.size();i++){
       vstream(buffer[off+table[i].first],compress(rhs._odata[so+table[i].second]));
     }
}

struct StencilEntry { 
  uint64_t _offset;
  uint64_t _byte_offset;
  uint16_t _is_local;
  uint16_t _permute;
  uint32_t _around_the_world; //256 bits, 32 bytes, 1/2 cacheline
};

template<class vobj,class cobj>
class CartesianStencil { // Stencil runs along coordinate axes only; NO diagonal fill in.
 public:

  typedef CartesianCommunicator::CommsRequest_t CommsRequest_t;
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
  };

  std::vector<Packet> Packets;

  int face_table_computed;
  std::vector<std::vector<std::pair<int,int> > > face_table ;
  
  void AddPacket(void *xmit,void * rcv, Integer to,Integer from,Integer bytes){
    Packet p;
    p.send_buf = xmit;
    p.recv_buf = rcv;
    p.to_rank  = to;
    p.from_rank= from;
    p.bytes    = bytes;
    comms_bytes+=2.0*bytes;
    Packets.push_back(p);
  }

  void CommunicateBegin(std::vector<std::vector<CommsRequest_t> > &reqs)
  {
    reqs.resize(Packets.size());
    commtime-=usecond();
    for(int i=0;i<Packets.size();i++){
	_grid->StencilSendToRecvFromBegin(reqs[i],
					  Packets[i].send_buf,
					  Packets[i].to_rank,
					  Packets[i].recv_buf,
					  Packets[i].from_rank,
					  Packets[i].bytes);
	/*
      }else{
	_grid->SendToRecvFromBegin(reqs[i],
				   Packets[i].send_buf,
				   Packets[i].to_rank,
				   Packets[i].recv_buf,
				   Packets[i].from_rank,
				   Packets[i].bytes);
      }
	*/
    }
    commtime+=usecond();
  }
  void CommunicateComplete(std::vector<std::vector<CommsRequest_t> > &reqs)
  {
    commtime-=usecond();

    for(int i=0;i<Packets.size();i++){
      //      if( ShmDirectCopy ) 
	_grid->StencilSendToRecvFromComplete(reqs[i]);
	//      else 
	//	_grid->SendToRecvFromComplete(reqs[i]);
    }
    commtime+=usecond();
  }

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
    Mergers.push_back(m);
  }

  void CommsMerge(void ) { 

    for(int i=0;i<Mergers.size();i++){	
      
      mergetime-=usecond();
PARALLEL_FOR_LOOP
      for(int o=0;o<Mergers[i].buffer_size;o++){
	merge1(Mergers[i].mpointer[o],Mergers[i].rpointers,o);
      }
      mergetime+=usecond();

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
  // Flat vector, change layout for cache friendly.
  Vector<StencilEntry>  _entries;
  
  void PrecomputeByteOffsets(void){
    for(int i=0;i<_entries.size();i++){
      if( _entries[i]._is_local ) {
	_entries[i]._byte_offset = _entries[i]._offset*sizeof(vobj);
      } else { 
	_entries[i]._byte_offset = _entries[i]._offset*sizeof(cobj);
      }
    }
  };

  inline StencilEntry * GetEntry(int &ptype,int point,int osite) { ptype = _permute_type[point]; return & _entries[point+_npoints*osite]; }
  inline uint64_t GetInfo(int &ptype,int &local,int &perm,int point,int ent,uint64_t base) {
    uint64_t cbase = (uint64_t)&u_recv_buf_p[0];
    local = _entries[ent]._is_local;
    perm  = _entries[ent]._permute;
    if (perm)  ptype = _permute_type[point]; 
    if (local) {
      return  base + _entries[ent]._byte_offset;
    } else {
      return cbase + _entries[ent]._byte_offset;
    }
  }
  inline uint64_t GetPFInfo(int ent,uint64_t base) {
    uint64_t cbase = (uint64_t)&u_recv_buf_p[0];
    int local = _entries[ent]._is_local;
    if (local) return  base + _entries[ent]._byte_offset;
    else       return cbase + _entries[ent]._byte_offset;
  }
  
  ///////////////////////////////////////////////////////////
  // Unified Comms buffers for all directions
  ///////////////////////////////////////////////////////////
  // Vectors that live on the symmetric heap in case of SHMEM
  //  std::vector<commVector<scalar_object> > u_simd_send_buf_hide;
  //  std::vector<commVector<scalar_object> > u_simd_recv_buf_hide;
  //  commVector<cobj>          u_send_buf_hide;
  //  commVector<cobj>          u_recv_buf_hide;

  // These are used; either SHM objects or refs to the above symmetric heap vectors
  // depending on comms target
  cobj* u_recv_buf_p;
  cobj* u_send_buf_p;
  std::vector<scalar_object *> u_simd_send_buf;
  std::vector<scalar_object *> u_simd_recv_buf;

  int u_comm_offset;
  int _unified_buffer_size;
  
  cobj *CommBuf(void) { return u_recv_buf_p; }

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
  double t_data;
  double t_table;
  double calls;
  
  void ZeroCounters(void) {
    gathertime = 0.;
    jointime = 0.;
    commtime = 0.;
    halogtime = 0.;
    mergetime = 0.;
    spintime = 0.;
    gathermtime = 0.;
    splicetime = 0.;
    nosplicetime = 0.;
    t_data = 0.0;
    t_table= 0.0;
    comms_bytes = 0.;
    calls = 0.;
  };
  
  void Report(void) {
#define PRINTIT(A)	\
 std::cout << GridLogMessage << " Stencil " << #A << " "<< A/calls<<std::endl;
    if ( calls > 0. ) {
      std::cout << GridLogMessage << " Stencil calls "<<calls<<std::endl;
      PRINTIT(halogtime);
      PRINTIT(gathertime);
      PRINTIT(gathermtime);
      PRINTIT(mergetime);
      if(comms_bytes>1.0){
	PRINTIT(comms_bytes);
	PRINTIT(commtime);
	std::cout << GridLogMessage << " Stencil " << comms_bytes/commtime/1000. << " GB/s "<<std::endl;
      }
      PRINTIT(jointime);
      PRINTIT(spintime);
      PRINTIT(splicetime);
      PRINTIT(nosplicetime);
      PRINTIT(t_table);
	   PRINTIT(t_data);
    }
  };
 #endif

 CartesianStencil(GridBase *grid,
		  int npoints,
		  int checkerboard,
		  const std::vector<int> &directions,
		  const std::vector<int> &distances) 
   :   _permute_type(npoints), _comm_buf_size(npoints)
  {
    face_table_computed=0;
    _npoints = npoints;
    _grid    = grid;
    _directions = directions;
    _distances  = distances;
    _unified_buffer_size=0;

    int osites  = _grid->oSites();
    
    _entries.resize(_npoints* osites);
    for(int ii=0;ii<npoints;ii++){
      
      int i = ii; // reverse direction to get SIMD comms done first
      int point = i;
      
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

    /////////////////////////////////////////////////////////////////////////////////
    // Try to allocate for receiving in a shared memory region, fall back to buffer
    /////////////////////////////////////////////////////////////////////////////////
    const int Nsimd = grid->Nsimd();

    _grid->ShmBufferFreeAll();

    u_simd_send_buf.resize(Nsimd);
    u_simd_recv_buf.resize(Nsimd);

    u_send_buf_p=(cobj *)_grid->ShmBufferMalloc(_unified_buffer_size*sizeof(cobj));
    u_recv_buf_p=(cobj *)_grid->ShmBufferMalloc(_unified_buffer_size*sizeof(cobj));
    for(int l=0;l<Nsimd;l++){
      u_simd_recv_buf[l] = (scalar_object *)_grid->ShmBufferMalloc(_unified_buffer_size*sizeof(scalar_object));
      u_simd_send_buf[l] = (scalar_object *)_grid->ShmBufferMalloc(_unified_buffer_size*sizeof(scalar_object));
    }

    PrecomputeByteOffsets();
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
	  int idx=point+(lo+o+b)*_npoints;
	  _entries[idx]._offset  =ro+o+b;
	  _entries[idx]._permute=permute;
	  _entries[idx]._is_local=1;
	  _entries[idx]._around_the_world=wrap;
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
	    int idx = point+(lo+o+b)*_npoints;
	    _entries[idx]._offset =ro+o+b;
	    _entries[idx]._is_local=1;
	    _entries[idx]._permute=permute;
	    _entries[idx]._around_the_world=wrap;
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
	  int idx=point+(so+o+b)*_npoints;
	  _entries[idx]._offset  =offset+(bo++);
	  _entries[idx]._is_local=0;
	  _entries[idx]._permute=0;
	  _entries[idx]._around_the_world=wrap;
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
	    int idx = point+(so+o+b)*_npoints;
	    _entries[idx]._offset  =offset+(bo++);
	    _entries[idx]._is_local=0;
	    _entries[idx]._permute =0;
	    _entries[idx]._around_the_world=wrap;
	  }
	}
	o +=_grid->_slice_stride[dimension];
      }
    }
  }
  
  template<class compressor> void HaloExchange(const Lattice<vobj> &source,compressor &compress) 
  {
    std::vector<std::vector<CommsRequest_t> > reqs;
    calls++;
    Mergers.resize(0);
    Packets.resize(0);
    _grid->StencilBarrier();
    HaloGather(source,compress);
    this->CommunicateBegin(reqs);
    _grid->StencilBarrier();
    this->CommunicateComplete(reqs);
    _grid->StencilBarrier();
    CommsMerge(); // spins
  }
  
  template<class compressor> void HaloGatherDir(const Lattice<vobj> &source,compressor &compress,int point,int & face_idx)
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
	  GatherSimd(source,dimension,shift,0x3,compress,face_idx);
	  splicetime+=usecond();
	} else { 
	  nosplicetime-=usecond();
	  Gather(source,dimension,shift,0x3,compress,face_idx);
	  nosplicetime+=usecond();
	}
      } else {
	if(splice_dim){
	  splicetime-=usecond();
	  GatherSimd(source,dimension,shift,0x1,compress,face_idx);// if checkerboard is unfavourable take two passes
	  GatherSimd(source,dimension,shift,0x2,compress,face_idx);// both with block stride loop iteration
	  splicetime+=usecond();
	} else {
	  nosplicetime-=usecond();
	  Gather(source,dimension,shift,0x1,compress,face_idx);
	  Gather(source,dimension,shift,0x2,compress,face_idx);
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
    
    u_comm_offset=0;
    
    // Gather all comms buffers
    int face_idx=0;
    for(int point = 0 ; point < _npoints; point++) {
      compress.Point(point);
      HaloGatherDir(source,compress,point,face_idx);
    }
    face_table_computed=1;
    
    assert(u_comm_offset==_unified_buffer_size);
    halogtime+=usecond();
  }
  
  template<class compressor>
  void Gather(const Lattice<vobj> &rhs,int dimension,int shift,int cbmask,compressor & compress,int &face_idx)
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
	int so  = sx*rhs._grid->_ostride[dimension]; // base offset for start of plane 
	if ( !face_table_computed ) {
	  t_table-=usecond();
	  face_table.resize(face_idx+1);
	  Gather_plane_simple_table_compute ((GridBase *)_grid,dimension,sx,cbmask,u_comm_offset,
					     face_table[face_idx]);
	  t_table+=usecond();
	}
	
	
	int rank           = _grid->_processor;
	int recv_from_rank;
	int xmit_to_rank;
	_grid->ShiftedRanks(dimension,comm_proc,xmit_to_rank,recv_from_rank);
	
	assert (xmit_to_rank   != _grid->ThisRank());
	assert (recv_from_rank != _grid->ThisRank());
	
	/////////////////////////////////////////////////////////
	// try the direct copy if possible
	/////////////////////////////////////////////////////////


	cobj *send_buf = (cobj *)_grid->ShmBufferTranslate(xmit_to_rank,u_recv_buf_p);
	if ( send_buf==NULL ) { 
	  send_buf = u_send_buf_p;
	}
	//	std::cout << " send_bufs  "<<std::hex<< send_buf <<" ubp "<<u_send_buf_p <<std::dec<<std::endl;
	t_data-=usecond();
	assert(u_send_buf_p!=NULL);
	assert(send_buf!=NULL);
	Gather_plane_simple_table         (face_table[face_idx],rhs,send_buf,compress,u_comm_offset,so);  face_idx++;
	t_data+=usecond();
	
	AddPacket((void *)&send_buf[u_comm_offset],
		  (void *)&u_recv_buf_p[u_comm_offset],
		  xmit_to_rank,
		  recv_from_rank,
		  bytes);

	gathertime+=usecond();
	u_comm_offset+=words;
      }
    }
  }
  
  template<class compressor>
  void  GatherSimd(const Lattice<vobj> &rhs,int dimension,int shift,int cbmask,compressor &compress,int & face_idx)
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
	  
	  if(nbr_proc){
	    
	    int recv_from_rank;
	    int xmit_to_rank;
	    
	    _grid->ShiftedRanks(dimension,nbr_proc,xmit_to_rank,recv_from_rank); 
 
	    scalar_object *shm = (scalar_object *) _grid->ShmBufferTranslate(recv_from_rank,sp);
	    //	    if ((ShmDirectCopy==0)||(shm==NULL)) { 
	    if (shm==NULL) { 
	      shm = rp;
	    } 
	    
	    // if Direct, StencilSendToRecvFrom will suppress copy to a peer on node
	    // assuming above pointer flip
	    AddPacket((void *)sp,(void *)rp,xmit_to_rank,recv_from_rank,bytes);
	    
	    rpointers[i] = shm;
	    
	  } else { 
	    
	    rpointers[i] = sp;
	    
	  }
	}

	AddMerge(&u_recv_buf_p[u_comm_offset],rpointers,buffer_size,Packets.size()-1);
	
	u_comm_offset     +=buffer_size;
      }
    }
  }
  
};
}

#if defined __GNUC__
 #pragma GCC diagnostic pop
#endif

#endif
