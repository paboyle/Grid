#ifndef _GRID_MPI_CSHIFT_H_
#define _GRID_MPI_CSHIFT_H_

#define MAX(x,y) ((x)>(y)?(x):(y))
#define MIN(x,y) ((x)>(y)?(y):(x))
//////////////////////////////////////////////////////////////////////////////////////////
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
//////////////////////////////////////////////////////////////////////////////////////////

      
/////////////////////////////////////////////////////////////
// Q. Further split this into separate sub functions?
/////////////////////////////////////////////////////////////

// CshiftCB_local
// CshiftCB_local_permute

// Cshift_comms_splice
// Cshift_comms
// Cshift_local
// Cshift_local_permute


friend Lattice<vobj> Cshift(Lattice<vobj> &rhs,int dimension,int shift)
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

friend void Cshift_comms(Lattice<vobj>& ret,Lattice<vobj> &rhs,int dimension,int shift)
{
  int sshift[2];

  sshift[0] = rhs._grid->CheckerBoardShift(rhs.checkerboard,dimension,shift,0);
  sshift[1] = rhs._grid->CheckerBoardShift(rhs.checkerboard,dimension,shift,1);

  if ( sshift[0] == sshift[1] ) {
    //    printf("Cshift_comms : single pass\n");
    Cshift_comms(ret,rhs,dimension,shift,0x3);
  } else {
    //    printf("Cshift_comms : two pass\n");
    //    printf("call1\n");
    Cshift_comms(ret,rhs,dimension,shift,0x1);// if checkerboard is unfavourable take two passes
    //    printf("call2\n");
    Cshift_comms(ret,rhs,dimension,shift,0x2);// both with block stride loop iteration
    //    printf("done\n");

  }
}

friend void Cshift_comms_simd(Lattice<vobj>& ret,Lattice<vobj> &rhs,int dimension,int shift)
{
  int sshift[2];

  sshift[0] = rhs._grid->CheckerBoardShift(rhs.checkerboard,dimension,shift,0);
  sshift[1] = rhs._grid->CheckerBoardShift(rhs.checkerboard,dimension,shift,1);

  if ( sshift[0] == sshift[1] ) {
    Cshift_comms_simd(ret,rhs,dimension,shift,0x3);
  } else {
    //    printf("call1 0x1 cb=even\n");
    Cshift_comms_simd(ret,rhs,dimension,shift,0x1);// if checkerboard is unfavourable take two passes
    //    printf("call2 0x2 cb=odd\n");
    Cshift_comms_simd(ret,rhs,dimension,shift,0x2);// both with block stride loop iteration
    //    printf("done\n");
  }
}


friend void Cshift_comms(Lattice<vobj> &ret,Lattice<vobj> &rhs,int dimension,int shift,int cbmask)
{
  typedef typename vobj::vector_type vector_type;
  typedef typename vobj::scalar_type scalar_type;

  SimdGrid *grid=rhs._grid;
  Lattice<vobj> temp(rhs._grid);

  int fd              = rhs._grid->_fdimensions[dimension];
  int rd              = rhs._grid->_rdimensions[dimension];
  int simd_layout     = rhs._grid->_simd_layout[dimension];
  int comm_dim        = rhs._grid->_processors[dimension] >1 ;
  assert(simd_layout==1);
  assert(comm_dim==1);
  assert(shift>=0);
  assert(shift<fd);
  
  // Packed gather sequence is clean
  int buffer_size = rhs._grid->_slice_nblock[dimension]*rhs._grid->_slice_block[dimension];
  std::vector<vobj,alignedAllocator<vobj> > send_buf(buffer_size);
  std::vector<vobj,alignedAllocator<vobj> > recv_buf(buffer_size);

  // This code could be simplified by multiple calls to single routine with extra params to
  // encapsulate the difference in the code paths.
  int cb= (cbmask==0x2)? 1 : 0;
  int sshift= rhs._grid->CheckerBoardShift(rhs.checkerboard,dimension,shift,cb);

  for(int x=0;x<rd;x++){       

    int offnode = ( x+sshift >= rd );
    int sx        = (x+sshift)%rd;
    int comm_proc = (x+sshift)/rd;
    
    if (!offnode) {
      //      printf("local x %d sshift %d offnode %d rd %d cb %d\n",x,sshift,offnode,rd,cb);
      Copy_plane(ret,rhs,dimension,x,sx,cbmask); 
    } else {

      int words = send_buf.size();
      if (cbmask != 0x3) words=words>>1;

      int bytes = words * sizeof(vobj);

      //      printf("nonlocal x %d sx %d sshift %d offnode %d rd %d cb %d cbmask %d rhscb %d comm_proc %d\n",
      //	     x,sx,sshift,offnode,rd,cb,cbmask,rhs.checkerboard,comm_proc);
      //      Copy_plane(temp,rhs,dimension,x,sx,cbmask); 

      // Bug found; cbmask may differ between sx plan and rx plane.
      Gather_plane_simple (rhs,send_buf,dimension,sx,cbmask);
      //      for(int i=0;i<MIN(words,8);i++){
      //	float *ptr = (float *)&send_buf[i];
      //	printf("send buf shift %d cbmask %d i %d %le\n",sshift,cbmask,i,*ptr);
      //      }
      //      Gather_plane_simple (rhs,send_buf,dimension,sx,cbmask^0x3);
      //      for(int i=0;i<MIN(words,8);i++){
      //	float *ptr = (float *)&send_buf[i];
      //	printf("send buf shift %d cbmask %d i %d %le\n",sshift,cbmask,i,*ptr);
      //      }
      //      recv_buf=send_buf;

      int rank           = grid->_processor;
      int recv_from_rank;
      int xmit_to_rank;
      grid->ShiftedRanks(dimension,comm_proc,xmit_to_rank,recv_from_rank);

      //      printf("bytes %d node %d sending to %d receiving from %d\n",bytes,rank,xmit_to_rank,recv_from_rank );
      grid->SendToRecvFrom((void *)&send_buf[0],
			   xmit_to_rank,
			   (void *)&recv_buf[0],
			   recv_from_rank,
			   bytes);

      Scatter_plane_simple (ret,recv_buf,dimension,x,cbmask);
    }
  }
}


friend void  Cshift_comms_simd(Lattice<vobj> &ret,Lattice<vobj> &rhs,int dimension,int shift,int cbmask)
{
  const int Nsimd = vector_type::Nsimd();
  SimdGrid *grid=rhs._grid;
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

  int permute_type=0;
  for(int d=0;d<dimension;d++){
    if (grid->_simd_layout[d]>1 ) permute_type++;
  }

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
  
  //  printf("cshift-comms-simd: shift = %d ; sshift = %d ; cbmask %d ; simd_layout %d\n",shift,sshift,cbmask,simd_layout);
  std::vector<int> comm_offnode(simd_layout);
  std::vector<int> comm_proc   (simd_layout);  //relative processor coord in dim=dimension

  // Strategy
  //
  //*  Loop over source planes
  //*    if any communication needed extract and send
  //*    if communication needed extract and send

  for(int x=0;x<rd;x++){       

    int comm_any = 0;
    for(int s=0;s<simd_layout;s++) {
      // does shift to "neighbour" takes us off node?
      // coordinates (reduce plane, simd_lane) of neighbour?
      // how many nodes away is this shift?
      // where we should send to?
      // where we should receive from?
      int shifted_x   = x+s*rd+sshift;
      comm_offnode[s] = shifted_x >= ld; 
      comm_any        = comm_any | comm_offnode[s];
      comm_proc[s]    = shifted_x/ld;     
      //      printf("rd %d x %d shifted %d s=%d comm_any %d\n",rd, x,shifted_x,s,comm_any);
    }
    
    int o    = 0;
    int bo   = x*grid->_ostride[dimension];
    int sx   = (x+sshift)%rd;

    // Need Convenience function in _grid. Move this in
    if ( comm_any ) {

      for(int i=0;i<Nsimd;i++){
	pointers[i] = (scalar_type *)&send_buf_extract[i][0];
      }
      Gather_plane_extract(rhs,pointers,dimension,sx,cbmask);
      //      for(int i=0;i<Nsimd;i++){
      //	printf("extracted %d %le\n",i,real(send_buf_extract[i][0]));
      //      }

      for(int i=0;i<Nsimd;i++){

	int s = grid->iCoordFromIsite(i,dimension);

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

	  //	  printf("Cshift_simd comms %d %le %le\n",i,real(recv_buf_extract[i][0]),real(send_buf_extract[i][0]));

	  rpointers[i] = (scalar_type *)&recv_buf_extract[i][0];

	} else { 

	  rpointers[i] = (scalar_type *)&send_buf_extract[i][0];
	  //	  printf("Cshift_simd local %d %le \n",i,real(send_buf_extract[i][0]));

	}

      }

      // Permute by swizzling pointers in merge
      int permute_slice=0;
      int lshift=sshift%ld;
      int wrap  =lshift/rd;
      int  num  =lshift%rd;

      if ( x< rd-num ) permute_slice=wrap;
      else permute_slice = 1-wrap;

      for(int i=0;i<vobj::vector_type::Nsimd();i++){
	if ( permute_slice ) {
	  pointers[i] = rpointers[permute_map[permute_type][i]];
	} else {
	  pointers[i] = rpointers[i];
	}
	//	printf("Cshift_simd perm %d num %d wrap %d swiz %d %le unswiz %le\n",permute_slice,num,wrap,i,real(pointers[i][0]),real(rpointers[i][0]));
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




#endif
