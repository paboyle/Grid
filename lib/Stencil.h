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

  class CartesianStencil { // Stencil runs along coordinate axes only; NO diagonal fill in.
  public:

      typedef uint32_t StencilInteger;

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

      int _unified_buffer_size;
      int _request_count;


      CartesianStencil(GridBase *grid,
		       int npoints,
		       int checkerboard,
		       const std::vector<int> &directions,
		       const std::vector<int> &distances);


      // Add to tables for various cases;  is this mistaken. only local if 1 proc in dim
      // Can this be avoided with simpler coding of comms?
      void Local     (int point, int dimension,int shift,int cbmask);
      void Comms     (int point, int dimension,int shift,int cbmask);
      void CopyPlane(int point, int dimension,int lplane,int rplane,int cbmask,int permute,int wrap);
      void ScatterPlane (int point,int dimension,int plane,int cbmask,int offset,int wrap);

      // Could allow a functional munging of the halo to another type during the comms.
      // this could implement the 16bit/32bit/64bit compression.
      template<class vobj,class cobj, class compressor> void 
	HaloExchange(const Lattice<vobj> &source,std::vector<cobj,alignedAllocator<cobj> > &u_comm_buf,compressor &compress)
      {
	// conformable(source._grid,_grid);
	assert(source._grid==_grid);
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
		GatherStartCommsSimd(source,dimension,shift,0x3,u_comm_buf,u_comm_offset,compress);
	      } else { 
		GatherStartComms(source,dimension,shift,0x3,u_comm_buf,u_comm_offset,compress);
	      }
	    } else {
	      if(splice_dim){
		GatherStartCommsSimd(source,dimension,shift,0x1,u_comm_buf,u_comm_offset,compress);// if checkerboard is unfavourable take two passes
		GatherStartCommsSimd(source,dimension,shift,0x2,u_comm_buf,u_comm_offset,compress);// both with block stride loop iteration
	      } else {
		GatherStartComms(source,dimension,shift,0x1,u_comm_buf,u_comm_offset,compress);
		GatherStartComms(source,dimension,shift,0x2,u_comm_buf,u_comm_offset,compress);
	      }
	    }
	  }
	}
      }

      template<class vobj,class cobj, class compressor> 
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
	  
	  std::vector<cobj,alignedAllocator<cobj> > send_buf(buffer_size); // hmm...
	  std::vector<cobj,alignedAllocator<cobj> > recv_buf(buffer_size);
	  
	  int cb= (cbmask==0x2)? Odd : Even;
	  int sshift= _grid->CheckerBoardShiftForCB(rhs.checkerboard,dimension,shift,cb);
	  
	  for(int x=0;x<rd;x++){       
	    
	    int sx        = (x+sshift)%rd;
	    int comm_proc = ((x+sshift)/rd)%pd;

	    if (comm_proc) {
	      
	      int words = send_buf.size();
	      if (cbmask != 0x3) words=words>>1;
	    
	      int bytes = words * sizeof(cobj);

	      Gather_plane_simple (rhs,send_buf,dimension,sx,cbmask,compress);

	      int rank           = _grid->_processor;
	      int recv_from_rank;
	      int xmit_to_rank;
	      _grid->ShiftedRanks(dimension,comm_proc,xmit_to_rank,recv_from_rank);
	      assert (xmit_to_rank != _grid->ThisRank());
	      assert (recv_from_rank != _grid->ThisRank());

	      //      FIXME Implement asynchronous send & also avoid buffer copy
	      _grid->SendToRecvFrom((void *)&send_buf[0],
				   xmit_to_rank,
				   (void *)&recv_buf[0],
				   recv_from_rank,
				   bytes);

	      for(int i=0;i<buffer_size;i++){
		u_comm_buf[u_comm_offset+i]=recv_buf[i];
	      }
	      u_comm_offset+=buffer_size;
	    }
	  }
	}


      template<class vobj,class cobj, class compressor> 
	void  GatherStartCommsSimd(const Lattice<vobj> &rhs,int dimension,int shift,int cbmask,
				   std::vector<cobj,alignedAllocator<cobj> > &u_comm_buf,
				   int &u_comm_offset,compressor &compress)
	{
	  const int Nsimd = _grid->Nsimd();

	  typedef typename cobj::vector_type vector_type;
	  typedef typename cobj::scalar_type scalar_type;
	  typedef typename cobj::scalar_object scalar_object;
	  
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

	  /*
	   * possibly slow to allocate
	   * Doesn't matter in this test, but may want to preallocate in the 
	   * dirac operators
	   */
	  std::vector<std::vector<scalar_object> > send_buf_extract(Nsimd,std::vector<scalar_object>(buffer_size) ); 
	  std::vector<std::vector<scalar_object> > recv_buf_extract(Nsimd,std::vector<scalar_object>(buffer_size) );
	  int bytes = buffer_size*sizeof(scalar_object);

	  std::vector<scalar_object *> pointers(Nsimd);  //
	  std::vector<scalar_object *> rpointers(Nsimd); // received pointers
	  
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
	      
	      Gather_plane_extract<cobj>(rhs,pointers,dimension,sx,cbmask,compress);

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
		  
		  _grid->SendToRecvFrom((void *)&send_buf_extract[nbr_lane][0],
					xmit_to_rank,
					(void *)&recv_buf_extract[i][0],
					recv_from_rank,
					bytes);
		  
		  rpointers[i] = &recv_buf_extract[i][0];

		} else { 
		  rpointers[i] = &send_buf_extract[nbr_lane][0];
		}
	      }

	      // Here we don't want to scatter, just place into a buffer.
	      for(int i=0;i<buffer_size;i++){
		assert(u_comm_offset+i<_unified_buffer_size);
		merge(u_comm_buf[u_comm_offset+i],rpointers,i);
	      }

	      u_comm_offset+=buffer_size;
	    }
	  }
	}
  };
}
#endif
