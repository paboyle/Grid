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


NAMESPACE_BEGIN(Grid); 
const int Cshift_verbose=0;
template<class vobj> Lattice<vobj> Cshift(const Lattice<vobj> &rhs,int dimension,int shift)
{
  typedef typename vobj::vector_type vector_type;
  typedef typename vobj::scalar_type scalar_type;

  Lattice<vobj> ret(rhs.Grid()); 
  
  int fd = rhs.Grid()->_fdimensions[dimension];
  int rd = rhs.Grid()->_rdimensions[dimension];

  // Map to always positive shift modulo global full dimension.
  shift = (shift+fd)%fd;

  ret.Checkerboard() = rhs.Grid()->CheckerBoardDestination(rhs.Checkerboard(),shift,dimension);
        
  // the permute type
  int simd_layout     = rhs.Grid()->_simd_layout[dimension];
  int comm_dim        = rhs.Grid()->_processors[dimension] >1 ;
  int splice_dim      = rhs.Grid()->_simd_layout[dimension]>1 && (comm_dim);

  RealD t1,t0;
  t0=usecond();
  if ( !comm_dim ) {
    //    std::cout << "CSHIFT: Cshift_local" <<std::endl;
    Cshift_local(ret,rhs,dimension,shift); // Handles checkerboarding
  } else if ( splice_dim ) {
    //    std::cout << "CSHIFT: Cshift_comms_simd call - splice_dim = " << splice_dim << " shift " << shift << " dimension = " << dimension << std::endl;
    Cshift_comms_simd(ret,rhs,dimension,shift);
  } else {
    //    std::cout << "CSHIFT: Cshift_comms" <<std::endl;
    Cshift_comms(ret,rhs,dimension,shift);
  }
  t1=usecond();
  if(Cshift_verbose) std::cout << GridLogPerformance << "Cshift took "<< (t1-t0)/1e3 << " ms"<<std::endl;
  return ret;
}
#if 1
template<class vobj> void Cshift_comms(Lattice<vobj>& ret,const Lattice<vobj> &rhs,int dimension,int shift)
{
  int sshift[2];

  sshift[0] = rhs.Grid()->CheckerBoardShiftForCB(rhs.Checkerboard(),dimension,shift,Even);
  sshift[1] = rhs.Grid()->CheckerBoardShiftForCB(rhs.Checkerboard(),dimension,shift,Odd);

  //  std::cout << "Cshift_comms dim "<<dimension<<"cb "<<rhs.Checkerboard()<<"shift "<<shift<<" sshift " << sshift[0]<<" "<<sshift[1]<<std::endl;
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

  sshift[0] = rhs.Grid()->CheckerBoardShiftForCB(rhs.Checkerboard(),dimension,shift,Even);
  sshift[1] = rhs.Grid()->CheckerBoardShiftForCB(rhs.Checkerboard(),dimension,shift,Odd);

  //  std::cout << "Cshift_comms_simd dim "<<dimension<<"cb "<<rhs.Checkerboard()<<"shift "<<shift<<" sshift " << sshift[0]<<" "<<sshift[1]<<std::endl;
  if ( sshift[0] == sshift[1] ) {
    //    std::cout << "Single pass Cshift_comms" <<std::endl;
    Cshift_comms_simd(ret,rhs,dimension,shift,0x3);
  } else {
    //    std::cout << "Two pass Cshift_comms" <<std::endl;
    Cshift_comms_simd(ret,rhs,dimension,shift,0x1);// if checkerboard is unfavourable take two passes
    Cshift_comms_simd(ret,rhs,dimension,shift,0x2);// both with block stride loop iteration
  }
}
template<class vobj> void Cshift_comms(Lattice<vobj> &ret,const Lattice<vobj> &rhs,int dimension,int shift,int cbmask)
{
  typedef typename vobj::vector_type vector_type;
  typedef typename vobj::scalar_type scalar_type;

  GridBase *grid=rhs.Grid();
  Lattice<vobj> temp(rhs.Grid());

  int fd              = rhs.Grid()->_fdimensions[dimension];
  int rd              = rhs.Grid()->_rdimensions[dimension];
  int pd              = rhs.Grid()->_processors[dimension];
  int simd_layout     = rhs.Grid()->_simd_layout[dimension];
  int comm_dim        = rhs.Grid()->_processors[dimension] >1 ;
  assert(simd_layout==1);
  assert(comm_dim==1);
  assert(shift>=0);
  assert(shift<fd);
  
  int buffer_size = rhs.Grid()->_slice_nblock[dimension]*rhs.Grid()->_slice_block[dimension];
  static deviceVector<vobj> send_buf; send_buf.resize(buffer_size);
  static deviceVector<vobj> recv_buf; recv_buf.resize(buffer_size);
    
  int cb= (cbmask==0x2)? Odd : Even;
  int sshift= rhs.Grid()->CheckerBoardShiftForCB(rhs.Checkerboard(),dimension,shift,cb);
  RealD tcopy=0.0;
  RealD tgather=0.0;
  RealD tscatter=0.0;
  RealD tcomms=0.0;
  uint64_t xbytes=0;
  for(int x=0;x<rd;x++){       

    int sx        =  (x+sshift)%rd;
    int comm_proc = ((x+sshift)/rd)%pd;
    
    if (comm_proc==0) {
      tcopy-=usecond();
      Copy_plane(ret,rhs,dimension,x,sx,cbmask); 
      tcopy+=usecond();
    } else {

      int words = buffer_size;
      if (cbmask != 0x3) words=words>>1;

      int bytes = words * sizeof(vobj);

      tgather-=usecond();
      Gather_plane_simple (rhs,send_buf,dimension,sx,cbmask);
      tgather+=usecond();

      //      int rank           = grid->_processor;
      int recv_from_rank;
      int xmit_to_rank;
      grid->ShiftedRanks(dimension,comm_proc,xmit_to_rank,recv_from_rank);
      
      tcomms-=usecond();
      grid->Barrier();

      grid->SendToRecvFrom((void *)&send_buf[0],
			   xmit_to_rank,
			   (void *)&recv_buf[0],
			   recv_from_rank,
			   bytes);
      xbytes+=bytes;
      grid->Barrier();
      tcomms+=usecond();

      tscatter-=usecond();
      Scatter_plane_simple (ret,recv_buf,dimension,x,cbmask);
      tscatter+=usecond();
    }
  }
  if (Cshift_verbose){
    std::cout << GridLogPerformance << " Cshift copy    "<<tcopy/1e3<<" ms"<<std::endl;
    std::cout << GridLogPerformance << " Cshift gather  "<<tgather/1e3<<" ms"<<std::endl;
    std::cout << GridLogPerformance << " Cshift scatter "<<tscatter/1e3<<" ms"<<std::endl;
    std::cout << GridLogPerformance << " Cshift comm    "<<tcomms/1e3<<" ms"<<std::endl;
    std::cout << GridLogPerformance << " Cshift BW      "<<(2.0*xbytes)/tcomms<<" MB/s "<<2*xbytes<< " Bytes "<<std::endl;
  }
}

template<class vobj> void  Cshift_comms_simd(Lattice<vobj> &ret,const Lattice<vobj> &rhs,int dimension,int shift,int cbmask)
{
  GridBase *grid=rhs.Grid();
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

  //  std::cout << "Cshift_comms_simd dim "<< dimension << " fd "<<fd<<" rd "<<rd
  //	    << " ld "<<ld<<" pd " << pd<<" simd_layout "<<simd_layout 
  //	    << " comm_dim " << comm_dim << " cbmask " << cbmask <<std::endl;

  assert(comm_dim==1);
  assert(simd_layout==2);
  assert(shift>=0);
  assert(shift<fd);

  RealD tcopy=0.0;
  RealD tgather=0.0;
  RealD tscatter=0.0;
  RealD tcomms=0.0;
  uint64_t xbytes=0;
  
  int permute_type=grid->PermuteType(dimension);

  ///////////////////////////////////////////////
  // Simd direction uses an extract/merge pair
  ///////////////////////////////////////////////
  int buffer_size = grid->_slice_nblock[dimension]*grid->_slice_block[dimension];
  //  int words = sizeof(vobj)/sizeof(vector_type);

  static std::vector<deviceVector<scalar_object> >  send_buf_extract; send_buf_extract.resize(Nsimd);
  static std::vector<deviceVector<scalar_object> >  recv_buf_extract; recv_buf_extract.resize(Nsimd);
  scalar_object *  recv_buf_extract_mpi;
  scalar_object *  send_buf_extract_mpi;
 
  for(int s=0;s<Nsimd;s++){
    send_buf_extract[s].resize(buffer_size);
    recv_buf_extract[s].resize(buffer_size);
  }

  int bytes = buffer_size*sizeof(scalar_object);

  ExtractPointerArray<scalar_object>  pointers(Nsimd); // 
  ExtractPointerArray<scalar_object> rpointers(Nsimd); // received pointers

  ///////////////////////////////////////////
  // Work out what to send where
  ///////////////////////////////////////////
  int cb    = (cbmask==0x2)? Odd : Even;
  int sshift= grid->CheckerBoardShiftForCB(rhs.Checkerboard(),dimension,shift,cb);

  // loop over outer coord planes orthog to dim
  for(int x=0;x<rd;x++){       

    // FIXME call local permute copy if none are offnode.
    for(int i=0;i<Nsimd;i++){       
      pointers[i] = &send_buf_extract[i][0];
    }
    int sx   = (x+sshift)%rd;
    tgather-=usecond();
    Gather_plane_extract(rhs,pointers,dimension,sx,cbmask);
    tgather+=usecond();

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

	tcomms-=usecond();
	grid->Barrier();

	send_buf_extract_mpi = &send_buf_extract[nbr_lane][0];
	recv_buf_extract_mpi = &recv_buf_extract[i][0];
	grid->SendToRecvFrom((void *)send_buf_extract_mpi,
			     xmit_to_rank,
			     (void *)recv_buf_extract_mpi,
			     recv_from_rank,
			     bytes);

	xbytes+=bytes;
	grid->Barrier();
	tcomms+=usecond();

	rpointers[i] = &recv_buf_extract[i][0];
      } else { 
	rpointers[i] = &send_buf_extract[nbr_lane][0];
      }

    }
    tscatter-=usecond();
    Scatter_plane_merge(ret,rpointers,dimension,x,cbmask);
    tscatter+=usecond();
  }
  if(Cshift_verbose){
    std::cout << GridLogPerformance << " Cshift (s) copy    "<<tcopy/1e3<<" ms"<<std::endl;
    std::cout << GridLogPerformance << " Cshift (s) gather  "<<tgather/1e3<<" ms"<<std::endl;
    std::cout << GridLogPerformance << " Cshift (s) scatter "<<tscatter/1e3<<" ms"<<std::endl;
    std::cout << GridLogPerformance << " Cshift (s) comm    "<<tcomms/1e3<<" ms"<<std::endl;
    std::cout << GridLogPerformance << " Cshift BW      "<<(2.0*xbytes)/tcomms<<" MB/s "<<2*xbytes<< " Bytes "<<std::endl;
  }
}
#else 
template<class vobj> void Cshift_comms(Lattice<vobj> &ret,const Lattice<vobj> &rhs,int dimension,int shift,int cbmask)
{
  typedef typename vobj::vector_type vector_type;
  typedef typename vobj::scalar_type scalar_type;

  GridBase *grid=rhs.Grid();
  Lattice<vobj> temp(rhs.Grid());

  int fd              = rhs.Grid()->_fdimensions[dimension];
  int rd              = rhs.Grid()->_rdimensions[dimension];
  int pd              = rhs.Grid()->_processors[dimension];
  int simd_layout     = rhs.Grid()->_simd_layout[dimension];
  int comm_dim        = rhs.Grid()->_processors[dimension] >1 ;
  assert(simd_layout==1);
  assert(comm_dim==1);
  assert(shift>=0);
  assert(shift<fd);
  RealD tcopy=0.0;
  RealD tgather=0.0;
  RealD tscatter=0.0;
  RealD tcomms=0.0;
  uint64_t xbytes=0;
  
  int buffer_size = rhs.Grid()->_slice_nblock[dimension]*rhs.Grid()->_slice_block[dimension];
  static cshiftVector<vobj> send_buf_v; send_buf_v.resize(buffer_size);
  static cshiftVector<vobj> recv_buf_v; recv_buf_v.resize(buffer_size);
  vobj *send_buf;
  vobj *recv_buf;
  {
    grid->ShmBufferFreeAll();
    size_t bytes = buffer_size*sizeof(vobj);
    send_buf=(vobj *)grid->ShmBufferMalloc(bytes);
    recv_buf=(vobj *)grid->ShmBufferMalloc(bytes);
  }
    
  int cb= (cbmask==0x2)? Odd : Even;
  int sshift= rhs.Grid()->CheckerBoardShiftForCB(rhs.Checkerboard(),dimension,shift,cb);

  for(int x=0;x<rd;x++){       

    int sx        =  (x+sshift)%rd;
    int comm_proc = ((x+sshift)/rd)%pd;
    
    if (comm_proc==0) {

      tcopy-=usecond();
      Copy_plane(ret,rhs,dimension,x,sx,cbmask); 
      tcopy+=usecond();

    } else {

      int words = buffer_size;
      if (cbmask != 0x3) words=words>>1;

      int bytes = words * sizeof(vobj);

      tgather-=usecond();
      Gather_plane_simple (rhs,send_buf_v,dimension,sx,cbmask);
      tgather+=usecond();

      //      int rank           = grid->_processor;
      int recv_from_rank;
      int xmit_to_rank;
      grid->ShiftedRanks(dimension,comm_proc,xmit_to_rank,recv_from_rank);


      tcomms-=usecond();
      //      grid->Barrier();

      acceleratorCopyDeviceToDevice((void *)&send_buf_v[0],(void *)&send_buf[0],bytes);
      grid->SendToRecvFrom((void *)&send_buf[0],
			   xmit_to_rank,
			   (void *)&recv_buf[0],
			   recv_from_rank,
			   bytes);
      xbytes+=bytes;
      acceleratorCopyDeviceToDevice((void *)&recv_buf[0],(void *)&recv_buf_v[0],bytes);

      //      grid->Barrier();
      tcomms+=usecond();

      tscatter-=usecond();
      Scatter_plane_simple (ret,recv_buf_v,dimension,x,cbmask);
      tscatter+=usecond();
    }
  }
  if(Cshift_verbose){
    std::cout << GridLogPerformance << " Cshift copy    "<<tcopy/1e3<<" ms"<<std::endl;
    std::cout << GridLogPerformance << " Cshift gather  "<<tgather/1e3<<" ms"<<std::endl;
    std::cout << GridLogPerformance << " Cshift scatter "<<tscatter/1e3<<" ms"<<std::endl;
    std::cout << GridLogPerformance << " Cshift comm    "<<tcomms/1e3<<" ms"<<std::endl;
    std::cout << GridLogPerformance << " Cshift BW      "<<(2.0*xbytes)/tcomms<<" MB/s "<<2*xbytes<< " Bytes "<<std::endl;
  }
}

template<class vobj> void  Cshift_comms_simd(Lattice<vobj> &ret,const Lattice<vobj> &rhs,int dimension,int shift,int cbmask)
{
  GridBase *grid=rhs.Grid();
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

  //std::cout << "Cshift_comms_simd dim "<< dimension << " fd "<<fd<<" rd "<<rd
  //    << " ld "<<ld<<" pd " << pd<<" simd_layout "<<simd_layout 
  //    << " comm_dim " << comm_dim << " cbmask " << cbmask <<std::endl;

  assert(comm_dim==1);
  assert(simd_layout==2);
  assert(shift>=0);
  assert(shift<fd);
  RealD tcopy=0.0;
  RealD tgather=0.0;
  RealD tscatter=0.0;
  RealD tcomms=0.0;
  uint64_t xbytes=0;

  int permute_type=grid->PermuteType(dimension);

  ///////////////////////////////////////////////
  // Simd direction uses an extract/merge pair
  ///////////////////////////////////////////////
  int buffer_size = grid->_slice_nblock[dimension]*grid->_slice_block[dimension];
  //  int words = sizeof(vobj)/sizeof(vector_type);

  static std::vector<cshiftVector<scalar_object> >  send_buf_extract; send_buf_extract.resize(Nsimd);
  static std::vector<cshiftVector<scalar_object> >  recv_buf_extract; recv_buf_extract.resize(Nsimd);
  scalar_object *  recv_buf_extract_mpi;
  scalar_object *  send_buf_extract_mpi;
  {
    size_t bytes = sizeof(scalar_object)*buffer_size;
    grid->ShmBufferFreeAll();
    send_buf_extract_mpi = (scalar_object *)grid->ShmBufferMalloc(bytes);
    recv_buf_extract_mpi = (scalar_object *)grid->ShmBufferMalloc(bytes);
  }
  for(int s=0;s<Nsimd;s++){
    send_buf_extract[s].resize(buffer_size);
    recv_buf_extract[s].resize(buffer_size);
  }

  int bytes = buffer_size*sizeof(scalar_object);

  ExtractPointerArray<scalar_object>  pointers(Nsimd); // 
  ExtractPointerArray<scalar_object> rpointers(Nsimd); // received pointers

  ///////////////////////////////////////////
  // Work out what to send where
  ///////////////////////////////////////////
  int cb    = (cbmask==0x2)? Odd : Even;
  int sshift= grid->CheckerBoardShiftForCB(rhs.Checkerboard(),dimension,shift,cb);

  // loop over outer coord planes orthog to dim
  for(int x=0;x<rd;x++){       

    // FIXME call local permute copy if none are offnode.
    for(int i=0;i<Nsimd;i++){       
      pointers[i] = &send_buf_extract[i][0];
    }
    tgather-=usecond();
    int sx   = (x+sshift)%rd;
    Gather_plane_extract(rhs,pointers,dimension,sx,cbmask);
    tgather+=usecond();

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

	tcomms-=usecond();
	//	grid->Barrier();

	acceleratorCopyDeviceToDevice((void *)&send_buf_extract[nbr_lane][0],(void *)send_buf_extract_mpi,bytes);
	grid->SendToRecvFrom((void *)send_buf_extract_mpi,
			     xmit_to_rank,
			     (void *)recv_buf_extract_mpi,
			     recv_from_rank,
			     bytes);
	acceleratorCopyDeviceToDevice((void *)recv_buf_extract_mpi,(void *)&recv_buf_extract[i][0],bytes);
	xbytes+=bytes;

	//	grid->Barrier();
	tcomms+=usecond();
	rpointers[i] = &recv_buf_extract[i][0];
      } else { 
	rpointers[i] = &send_buf_extract[nbr_lane][0];
      }

    }
    tscatter-=usecond();
    Scatter_plane_merge(ret,rpointers,dimension,x,cbmask);
    tscatter+=usecond();

  }
  if(Cshift_verbose){
    std::cout << GridLogPerformance << " Cshift (s) copy    "<<tcopy/1e3<<" ms"<<std::endl;
    std::cout << GridLogPerformance << " Cshift (s) gather  "<<tgather/1e3<<" ms"<<std::endl;
    std::cout << GridLogPerformance << " Cshift (s) scatter "<<tscatter/1e3<<" ms"<<std::endl;
    std::cout << GridLogPerformance << " Cshift (s) comm    "<<tcomms/1e3<<" ms"<<std::endl;
    std::cout << GridLogPerformance << " Cshift BW      "<<(2.0*xbytes)/tcomms<<" MB/s"<<std::endl;
  }
}
#endif

NAMESPACE_END(Grid); 

#endif
