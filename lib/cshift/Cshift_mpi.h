    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/cshift/Cshift_mpi.h

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
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
#ifndef _GRID_CSHIFT_MPI_H_
#define _GRID_CSHIFT_MPI_H_


namespace Grid { 

template<class vobj> Lattice<vobj> Cshift(const Lattice<vobj> &rhs,int dimension,int shift)
{
  typedef typename vobj::vector_type vector_type;
  typedef typename vobj::scalar_type scalar_type;

  Lattice<vobj> ret(rhs._grid); 
  
  int fd = rhs._grid->_fdimensions[dimension];
  int rd = rhs._grid->_rdimensions[dimension];

  // Map to always positive shift modulo global full dimension.
  shift = (shift+fd)%fd;

  ret.checkerboard = rhs._grid->CheckerBoardDestination(rhs.checkerboard,shift,dimension);
        
  // the permute type
  int simd_layout     = rhs._grid->_simd_layout[dimension];
  int comm_dim        = rhs._grid->_processors[dimension] >1 ;
  int splice_dim      = rhs._grid->_simd_layout[dimension]>1 && (comm_dim);


  if ( !comm_dim ) {
    //    std::cout << "Cshift_local" <<std::endl;
    Cshift_local(ret,rhs,dimension,shift); // Handles checkerboarding
  } else if ( splice_dim ) {
    //    std::cout << "Cshift_comms_simd" <<std::endl;
    Cshift_comms_simd(ret,rhs,dimension,shift);
  } else {
    //    std::cout << "Cshift_comms" <<std::endl;
    Cshift_comms(ret,rhs,dimension,shift);
  }
  return ret;
}

template<class vobj> void Cshift_comms(Lattice<vobj>& ret,const Lattice<vobj> &rhs,int dimension,int shift)
{
  int sshift[2];

  sshift[0] = rhs._grid->CheckerBoardShiftForCB(rhs.checkerboard,dimension,shift,Even);
  sshift[1] = rhs._grid->CheckerBoardShiftForCB(rhs.checkerboard,dimension,shift,Odd);

  //  std::cout << "Cshift_comms dim "<<dimension<<"cb "<<rhs.checkerboard<<"shift "<<shift<<" sshift " << sshift[0]<<" "<<sshift[1]<<std::endl;

  if ( sshift[0] == sshift[1] ) {
    //    std::cout << "Single pass Cshift_comms" <<std::endl;
    Cshift_comms(ret,rhs,dimension,shift,0x3);
  } else {
    //    std::cout << "Two pass Cshift_comms" <<std::endl;
    Cshift_comms(ret,rhs,dimension,shift,0x1);// if checkerboard is unfavourable take two passes
    Cshift_comms(ret,rhs,dimension,shift,0x2);// both with block stride loop iteration
  }
}

template<class vobj> void Cshift_comms_simd(Lattice<vobj>& ret,const Lattice<vobj> &rhs,int dimension,int shift)
{
  int sshift[2];

  sshift[0] = rhs._grid->CheckerBoardShiftForCB(rhs.checkerboard,dimension,shift,Even);
  sshift[1] = rhs._grid->CheckerBoardShiftForCB(rhs.checkerboard,dimension,shift,Odd);

  if ( sshift[0] == sshift[1] ) {
    Cshift_comms_simd(ret,rhs,dimension,shift,0x3);
  } else {
    Cshift_comms_simd(ret,rhs,dimension,shift,0x1);// if checkerboard is unfavourable take two passes
    Cshift_comms_simd(ret,rhs,dimension,shift,0x2);// both with block stride loop iteration
  }
}

template<class vobj> void Cshift_comms(Lattice<vobj> &ret,const Lattice<vobj> &rhs,int dimension,int shift,int cbmask)
{
  typedef typename vobj::vector_type vector_type;
  typedef typename vobj::scalar_type scalar_type;

  GridBase *grid=rhs._grid;
  Lattice<vobj> temp(rhs._grid);

  int fd              = rhs._grid->_fdimensions[dimension];
  int rd              = rhs._grid->_rdimensions[dimension];
  int pd              = rhs._grid->_processors[dimension];
  int simd_layout     = rhs._grid->_simd_layout[dimension];
  int comm_dim        = rhs._grid->_processors[dimension] >1 ;
  assert(simd_layout==1);
  assert(comm_dim==1);
  assert(shift>=0);
  assert(shift<fd);
  
  int buffer_size = rhs._grid->_slice_nblock[dimension]*rhs._grid->_slice_block[dimension];
  std::vector<vobj,alignedAllocator<vobj> > send_buf(buffer_size);
  std::vector<vobj,alignedAllocator<vobj> > recv_buf(buffer_size);

  int cb= (cbmask==0x2)? Odd : Even;
  int sshift= rhs._grid->CheckerBoardShiftForCB(rhs.checkerboard,dimension,shift,cb);

  for(int x=0;x<rd;x++){       

    int sx        =  (x+sshift)%rd;
    int comm_proc = ((x+sshift)/rd)%pd;
    
    if (comm_proc==0) {

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

      //      for(int i=0;i<words;i++){
      //	std::cout << "SendRecv ["<<i<<"] snd "<<send_buf[i]<<" rcv " << recv_buf[i] << "  0x" << cbmask<<std::endl;
      //      }
      Scatter_plane_simple (ret,recv_buf,dimension,x,cbmask);
    }
  }
}

template<class vobj> void  Cshift_comms_simd(Lattice<vobj> &ret,const Lattice<vobj> &rhs,int dimension,int shift,int cbmask)
{
  GridBase *grid=rhs._grid;
  const int Nsimd = grid->Nsimd();
  typedef typename vobj::vector_type vector_type;
  typedef typename vobj::scalar_object scalar_object;
  typedef typename vobj::scalar_type scalar_type;
   
  int fd = grid->_fdimensions[dimension];
  int rd = grid->_rdimensions[dimension];
  int ld = grid->_ldimensions[dimension];
  int pd = grid->_processors[dimension];
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

  std::vector<std::vector<scalar_object> > send_buf_extract(Nsimd,std::vector<scalar_object>(buffer_size) );
  std::vector<std::vector<scalar_object> > recv_buf_extract(Nsimd,std::vector<scalar_object>(buffer_size) );
  int bytes = buffer_size*sizeof(scalar_object);

  std::vector<scalar_object *>  pointers(Nsimd); // 
  std::vector<scalar_object *> rpointers(Nsimd); // received pointers

  ///////////////////////////////////////////
  // Work out what to send where
  ///////////////////////////////////////////
  int cb    = (cbmask==0x2)? Odd : Even;
  int sshift= grid->CheckerBoardShiftForCB(rhs.checkerboard,dimension,shift,cb);

  // loop over outer coord planes orthog to dim
  for(int x=0;x<rd;x++){       

    // FIXME call local permute copy if none are offnode.
    for(int i=0;i<Nsimd;i++){       
      pointers[i] = &send_buf_extract[i][0];
    }
    int sx   = (x+sshift)%rd;
    Gather_plane_extract(rhs,pointers,dimension,sx,cbmask);

    for(int i=0;i<Nsimd;i++){
      
      int inner_bit = (Nsimd>>(permute_type+1));
      int ic= (i&inner_bit)? 1:0;

      int my_coor          = rd*ic + x;
      int nbr_coor         = my_coor+sshift;
      int nbr_proc = ((nbr_coor)/ld) % pd;// relative shift in processors

      int nbr_ic   = (nbr_coor%ld)/rd;    // inner coord of peer
      int nbr_ox   = (nbr_coor%rd);       // outer coord of peer
      int nbr_lane = (i&(~inner_bit));

      int recv_from_rank;
      int xmit_to_rank;

      if (nbr_ic) nbr_lane|=inner_bit;

      assert (sx == nbr_ox);

      if(nbr_proc){
	grid->ShiftedRanks(dimension,nbr_proc,xmit_to_rank,recv_from_rank); 

	grid->SendToRecvFrom((void *)&send_buf_extract[nbr_lane][0],
			     xmit_to_rank,
			     (void *)&recv_buf_extract[i][0],
			     recv_from_rank,
			     bytes);

	rpointers[i] = &recv_buf_extract[i][0];
      } else { 
	rpointers[i] = &send_buf_extract[nbr_lane][0];
      }

    }
    Scatter_plane_merge(ret,rpointers,dimension,x,cbmask);
  }

 }
}
#endif
