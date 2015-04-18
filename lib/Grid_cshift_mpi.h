#ifndef _GRID_CSHIFT_MPI_H_
#define _GRID_CSHIFT_MPI_H_

#ifndef MAX
#define MAX(x,y) ((x)>(y)?(x):(y))
#define MIN(x,y) ((x)>(y)?(y):(x))
#endif

namespace Grid { 

template<class vobj> Lattice<vobj> Cshift(Lattice<vobj> &rhs,int dimension,int shift)
{
  typedef typename vobj::vector_type vector_type;
  typedef typename vobj::scalar_type scalar_type;

  Lattice<vobj> ret(rhs._grid);
  
  int fd = rhs._grid->_fdimensions[dimension];
  int rd = rhs._grid->_rdimensions[dimension];

  // Map to always positive shift modulo global full dimension.
  shift = (shift+fd)%fd;

  ret.checkerboard = rhs._grid->CheckerBoardDestination(rhs.checkerboard,shift);
        
  // the permute type
  int simd_layout     = rhs._grid->_simd_layout[dimension];
  int comm_dim        = rhs._grid->_processors[dimension] >1 ;
  int splice_dim      = rhs._grid->_simd_layout[dimension]>1 && (comm_dim);


  if ( !comm_dim ) {
    Cshift_local(ret,rhs,dimension,shift); // Handles checkerboarding
  } else if ( splice_dim ) {
    Cshift_comms_simd(ret,rhs,dimension,shift);
  } else {
    Cshift_comms(ret,rhs,dimension,shift);
  }
  return ret;
}

template<class vobj> void Cshift_comms(Lattice<vobj>& ret,Lattice<vobj> &rhs,int dimension,int shift)
{
  int sshift[2];

  sshift[0] = rhs._grid->CheckerBoardShift(rhs.checkerboard,dimension,shift,0);
  sshift[1] = rhs._grid->CheckerBoardShift(rhs.checkerboard,dimension,shift,1);

  if ( sshift[0] == sshift[1] ) {
    Cshift_comms(ret,rhs,dimension,shift,0x3);
  } else {
    Cshift_comms(ret,rhs,dimension,shift,0x1);// if checkerboard is unfavourable take two passes
    Cshift_comms(ret,rhs,dimension,shift,0x2);// both with block stride loop iteration
  }
}

template<class vobj> void Cshift_comms_simd(Lattice<vobj>& ret,Lattice<vobj> &rhs,int dimension,int shift)
{
  int sshift[2];

  sshift[0] = rhs._grid->CheckerBoardShift(rhs.checkerboard,dimension,shift,0);
  sshift[1] = rhs._grid->CheckerBoardShift(rhs.checkerboard,dimension,shift,1);

  if ( sshift[0] == sshift[1] ) {
    Cshift_comms_simd(ret,rhs,dimension,shift,0x3);
  } else {
    Cshift_comms_simd(ret,rhs,dimension,shift,0x1);// if checkerboard is unfavourable take two passes
    Cshift_comms_simd(ret,rhs,dimension,shift,0x2);// both with block stride loop iteration
  }
}

template<class vobj> void Cshift_comms(Lattice<vobj> &ret,Lattice<vobj> &rhs,int dimension,int shift,int cbmask)
{
  typedef typename vobj::vector_type vector_type;
  typedef typename vobj::scalar_type scalar_type;

  GridBase *grid=rhs._grid;
  Lattice<vobj> temp(rhs._grid);

  int fd              = rhs._grid->_fdimensions[dimension];
  int rd              = rhs._grid->_rdimensions[dimension];
  int simd_layout     = rhs._grid->_simd_layout[dimension];
  int comm_dim        = rhs._grid->_processors[dimension] >1 ;
  assert(simd_layout==1);
  assert(comm_dim==1);
  assert(shift>=0);
  assert(shift<fd);
  
  int buffer_size = rhs._grid->_slice_nblock[dimension]*rhs._grid->_slice_block[dimension];
  std::vector<vobj,alignedAllocator<vobj> > send_buf(buffer_size);
  std::vector<vobj,alignedAllocator<vobj> > recv_buf(buffer_size);

  int cb= (cbmask==0x2)? 1 : 0;
  int sshift= rhs._grid->CheckerBoardShift(rhs.checkerboard,dimension,shift,cb);

  for(int x=0;x<rd;x++){       

    int offnode = ( x+sshift >= rd );
    int sx        = (x+sshift)%rd;
    int comm_proc = (x+sshift)/rd;
    
    if (!offnode) {

      Copy_plane(ret,rhs,dimension,x,sx,cbmask); 

    } else {

      int words = send_buf.size();
      if (cbmask != 0x3) words=words>>1;

      int bytes = words * sizeof(vobj);

      Gather_plane_simple (rhs,send_buf,dimension,sx,cbmask);

      int rank           = grid->_processor;
      int recv_from_rank;
      int xmit_to_rank;
      grid->ShiftedRanks(dimension,comm_proc,xmit_to_rank,recv_from_rank);

      grid->SendToRecvFrom((void *)&send_buf[0],
			   xmit_to_rank,
			   (void *)&recv_buf[0],
			   recv_from_rank,
			   bytes);

      Scatter_plane_simple (ret,recv_buf,dimension,x,cbmask);
    }
  }
}

template<class vobj> void  Cshift_comms_simd(Lattice<vobj> &ret,Lattice<vobj> &rhs,int dimension,int shift,int cbmask)
{
  GridBase *grid=rhs._grid;
  const int Nsimd = grid->Nsimd();
  typedef typename vobj::vector_type vector_type;
  typedef typename vobj::scalar_type scalar_type;
   
  int fd = grid->_fdimensions[dimension];
  int rd = grid->_rdimensions[dimension];
  int ld = grid->_ldimensions[dimension];
  int simd_layout     = grid->_simd_layout[dimension];
  int comm_dim        = grid->_processors[dimension] >1 ;

  assert(comm_dim==1);
  assert(simd_layout==2);
  assert(shift>=0);
  assert(shift<fd);

  int permute_type=grid->PermuteType(dimension);

  ///////////////////////////////////////////////
  // Simd direction uses an extract/merge pair
  ///////////////////////////////////////////////
  int buffer_size = grid->_slice_nblock[dimension]*grid->_slice_block[dimension];
  int words = sizeof(vobj)/sizeof(vector_type);

  std::vector<std::vector<scalar_type> > send_buf_extract(Nsimd,std::vector<scalar_type>(buffer_size*words) );
  std::vector<std::vector<scalar_type> > recv_buf_extract(Nsimd,std::vector<scalar_type>(buffer_size*words) );
  int bytes = buffer_size*words*sizeof(scalar_type);

  std::vector<scalar_type *> pointers(Nsimd);  // 
  std::vector<scalar_type *> rpointers(Nsimd); // received pointers

  ///////////////////////////////////////////
  // Work out what to send where
  ///////////////////////////////////////////

  int cb    = (cbmask==0x2)? 1 : 0;
  int sshift= grid->CheckerBoardShift(rhs.checkerboard,dimension,shift,cb);
  
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
	  

	  grid->SendToRecvFrom((void *)&send_buf_extract[i][0],
			    xmit_to_rank,
			    (void *)&recv_buf_extract[i][0],
			    recv_from_rank,
			    bytes);

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

      Scatter_plane_merge(ret,pointers,dimension,x,cbmask);

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
}
#endif
