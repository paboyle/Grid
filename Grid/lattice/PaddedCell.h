/*************************************************************************************
    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/lattice/PaddedCell.h

    Copyright (C) 2019

Author: Peter Boyle pboyle@bnl.gov

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
#pragma once

#include<Grid/cshift/Cshift.h>

NAMESPACE_BEGIN(Grid);

//Allow the user to specify how the C-shift is performed, e.g. to respect the appropriate boundary conditions
template<typename vobj>
struct CshiftImplBase{
  virtual Lattice<vobj> Cshift(const Lattice<vobj> &in, int dir, int shift) const = 0;
  virtual ~CshiftImplBase(){}
};
template<typename vobj>
struct CshiftImplDefault: public CshiftImplBase<vobj>{
  Lattice<vobj> Cshift(const Lattice<vobj> &in, int dir, int shift) const override{ return Grid::Cshift(in,dir,shift); }
};
template<typename Gimpl>
struct CshiftImplGauge: public CshiftImplBase<typename Gimpl::GaugeLinkField::vector_object>{
  typename Gimpl::GaugeLinkField Cshift(const typename Gimpl::GaugeLinkField &in, int dir, int shift) const override{ return Gimpl::CshiftLink(in,dir,shift); }
};  


/*
 *
 * TODO: 
 *  -- address elementsof vobj via thread block in Scatter/Gather
 *  -- overlap comms with motion in Face_exchange
 *
 */

template<class vobj> inline void ScatterSlice(const cshiftVector<vobj> &buf,
					      Lattice<vobj> &lat,
					      int x,
					      int dim,
					      int offset=0)
{
  const int Nsimd=vobj::Nsimd();
  typedef typename vobj::scalar_object sobj;
  typedef typename vobj::scalar_type scalar_type;
  typedef typename vobj::vector_type vector_type;

  GridBase *grid = lat.Grid();
  Coordinate simd = grid->_simd_layout;
  int Nd          = grid->Nd();
  int block       = grid->_slice_block[dim];
  int stride      = grid->_slice_stride[dim];
  int nblock      = grid->_slice_nblock[dim];
  int rd          = grid->_rdimensions[dim];

  int ox = x%rd;
  int ix = x/rd;

  int isites = 1; for(int d=0;d<Nd;d++) if( d!=dim) isites*=simd[d];

  Coordinate rsimd= simd;  rsimd[dim]=1; // maybe reduce Nsimd

  int rNsimd = 1; for(int d=0;d<Nd;d++) rNsimd*=rsimd[d];
  int rNsimda= Nsimd/simd[dim]; // should be equal
  assert(rNsimda==rNsimd);
  int face_ovol=block*nblock;

  //  assert(buf.size()==face_ovol*rNsimd);

  /*This will work GPU ONLY unless rNsimd is put in the lexico index*/
  //Let's make it work on GPU and then make a special accelerator_for that
  //doesn't hide the SIMD direction and keeps explicit in the threadIdx
  //for cross platform
  // FIXME -- can put internal indices into thread loop
  auto buf_p = & buf[0];
  autoView(lat_v, lat, AcceleratorWrite);
  accelerator_for(ss, face_ovol/simd[dim],Nsimd,{

    // scalar layout won't coalesce
#ifdef GRID_SIMT
      {
	int blane=acceleratorSIMTlane(Nsimd); // buffer lane
#else
      for(int blane=0;blane<Nsimd;blane++) {
#endif
	int olane=blane%rNsimd;               // reduced lattice lane
	int obit =blane/rNsimd;

	///////////////////////////////////////////////////////////////
	// osite -- potentially one bit from simd in the buffer: (ss<<1)|obit
	///////////////////////////////////////////////////////////////
	int ssp = ss*simd[dim]+obit;
	int b    = ssp%block;
	int n    = ssp/block;
	int osite= b+n*stride + ox*block;
	
	////////////////////////////////////////////
	// isite -- map lane within buffer to lane within lattice
	////////////////////////////////////////////
	Coordinate icoor;
	int lane;
	Lexicographic::CoorFromIndex(icoor,olane,rsimd);
	icoor[dim]=ix;
	Lexicographic::IndexFromCoor(icoor,lane,simd);
	
	///////////////////////////////////////////
	// Transfer into lattice - will coalesce
	///////////////////////////////////////////
	//	sobj obj = extractLane(blane,buf_p[ss+offset]);
	//	insertLane(lane,lat_v[osite],obj);
	const int words=sizeof(vobj)/sizeof(vector_type);
	vector_type * from = (vector_type *)&buf_p[ss+offset];
	vector_type * to   = (vector_type *)&lat_v[osite];
	scalar_type stmp;
	for(int w=0;w<words;w++){
	  stmp = getlane(from[w], blane);
	  putlane(to[w], stmp, lane);
	}
      }
  });
}

template<class vobj> inline void GatherSlice(cshiftVector<vobj> &buf,
					     const Lattice<vobj> &lat,
					     int x,
					     int dim,
					     int offset=0)
{
  const int Nsimd=vobj::Nsimd();
  typedef typename vobj::scalar_object sobj;
  typedef typename vobj::scalar_type scalar_type;
  typedef typename vobj::vector_type vector_type;

  autoView(lat_v, lat, AcceleratorRead);

  GridBase *grid = lat.Grid();
  Coordinate simd = grid->_simd_layout;
  int Nd          = grid->Nd();
  int block       = grid->_slice_block[dim];
  int stride      = grid->_slice_stride[dim];
  int nblock      = grid->_slice_nblock[dim];
  int rd          = grid->_rdimensions[dim];

  int ox = x%rd;
  int ix = x/rd;

  int isites = 1; for(int d=0;d<Nd;d++) if( d!=dim) isites*=simd[d];

  Coordinate rsimd= simd;  rsimd[dim]=1; // maybe reduce Nsimd

  int rNsimd = 1; for(int d=0;d<Nd;d++) rNsimd*=rsimd[d];
  
  int face_ovol=block*nblock;

  //  assert(buf.size()==face_ovol*rNsimd);

  /*This will work GPU ONLY unless rNsimd is put in the lexico index*/
  //Let's make it work on GPU and then make a special accelerator_for that
  //doesn't hide the SIMD direction and keeps explicit in the threadIdx
  //for cross platform
  //For CPU perhaps just run a loop over Nsimd
  auto buf_p = & buf[0];
  accelerator_for(ss, face_ovol/simd[dim],Nsimd,{

    // scalar layout won't coalesce
#ifdef GRID_SIMT
      {
	int blane=acceleratorSIMTlane(Nsimd); // buffer lane
#else
      for(int blane=0;blane<Nsimd;blane++) {
#endif
	int olane=blane%rNsimd;               // reduced lattice lane
	int obit =blane/rNsimd;
	
	////////////////////////////////////////////
	// osite
	////////////////////////////////////////////
	int ssp = ss*simd[dim]+obit;
	int b    = ssp%block;
	int n    = ssp/block;
	int osite= b+n*stride + ox*block;

	////////////////////////////////////////////
	// isite -- map lane within buffer to lane within lattice
	////////////////////////////////////////////
	Coordinate icoor;
	int lane;
	Lexicographic::CoorFromIndex(icoor,olane,rsimd);
	icoor[dim]=ix;
	Lexicographic::IndexFromCoor(icoor,lane,simd);
	
	///////////////////////////////////////////
	// Take out of lattice
	///////////////////////////////////////////
	//	sobj obj = extractLane(lane,lat_v[osite]);
	//	insertLane(blane,buf_p[ss+offset],obj);
	const int words=sizeof(vobj)/sizeof(vector_type);
	vector_type * to    = (vector_type *)&buf_p[ss+offset];
	vector_type * from  = (vector_type *)&lat_v[osite];
	scalar_type stmp;
	for(int w=0;w<words;w++){
	  stmp = getlane(from[w], lane);
	  putlane(to[w], stmp, blane);
	}
      }
  });
}


class PaddedCell {
public:
  GridCartesian * unpadded_grid;
  int dims;
  int depth;
  std::vector<GridCartesian *> grids;

  ~PaddedCell()
  {
    DeleteGrids();
  }
  PaddedCell(int _depth,GridCartesian *_grid)
  {
    unpadded_grid = _grid;
    depth=_depth;
    dims=_grid->Nd();
    AllocateGrids();
    Coordinate local     =unpadded_grid->LocalDimensions();
    Coordinate procs     =unpadded_grid->ProcessorGrid();
    for(int d=0;d<dims;d++){
      if ( procs[d] > 1 ) assert(local[d]>=depth);
    }
  }
  void DeleteGrids(void)
  {
    Coordinate processors=unpadded_grid->_processors;
    for(int d=0;d<grids.size();d++){
      if ( processors[d] > 1 ) { 
	delete grids[d];
      }
    }
    grids.resize(0);
  };
  void AllocateGrids(void)
  {
    Coordinate local     =unpadded_grid->LocalDimensions();
    Coordinate simd      =unpadded_grid->_simd_layout;
    Coordinate processors=unpadded_grid->_processors;
    Coordinate plocal    =unpadded_grid->LocalDimensions();
    Coordinate global(dims);
    GridCartesian *old_grid = unpadded_grid;
    // expand up one dim at a time
    for(int d=0;d<dims;d++){

      if ( processors[d] > 1 ) { 
	plocal[d] += 2*depth; 
      
	for(int d=0;d<dims;d++){
	  global[d] = plocal[d]*processors[d];
	}

	old_grid = new GridCartesian(global,simd,processors);
      }
      grids.push_back(old_grid);
    }
  };
  template<class vobj>
  inline Lattice<vobj> Extract(const Lattice<vobj> &in) const
  {
    Coordinate processors=unpadded_grid->_processors;

    Lattice<vobj> out(unpadded_grid);

    Coordinate local     =unpadded_grid->LocalDimensions();
    // depends on the MPI spread      
    Coordinate fll(dims,depth);
    Coordinate tll(dims,0); // depends on the MPI spread
    for(int d=0;d<dims;d++){
      if( processors[d]==1 ) fll[d]=0;
    }
    localCopyRegion(in,out,fll,tll,local);
    return out;
  }
  template<class vobj>
  inline Lattice<vobj> Exchange(const Lattice<vobj> &in, const CshiftImplBase<vobj> &cshift = CshiftImplDefault<vobj>()) const
  {
    GridBase *old_grid = in.Grid();
    int dims = old_grid->Nd();
    Lattice<vobj> tmp = in;
    for(int d=0;d<dims;d++){
      tmp = Expand(d,tmp,cshift); // rvalue && assignment
    }
    return tmp;
  }
  template<class vobj>
  inline Lattice<vobj> ExchangePeriodic(const Lattice<vobj> &in) const
  {
    GridBase *old_grid = in.Grid();
    int dims = old_grid->Nd();
    Lattice<vobj> tmp = in;
    for(int d=0;d<dims;d++){
      tmp = ExpandPeriodic(d,tmp); // rvalue && assignment
    }
    return tmp;
  }
  // expand up one dim at a time
  template<class vobj>
  inline Lattice<vobj> Expand(int dim, const Lattice<vobj> &in, const CshiftImplBase<vobj> &cshift = CshiftImplDefault<vobj>()) const
  {
    Coordinate processors=unpadded_grid->_processors;
    GridBase *old_grid = in.Grid();
    GridCartesian *new_grid = grids[dim];//These are new grids
    Lattice<vobj>  padded(new_grid);
    Lattice<vobj> shifted(old_grid);    
    Coordinate local     =old_grid->LocalDimensions();
    Coordinate plocal    =new_grid->LocalDimensions();
    if(dim==0) conformable(old_grid,unpadded_grid);
    else       conformable(old_grid,grids[dim-1]);

    double tins=0, tshift=0;

    int islocal = 0 ;
    if ( processors[dim] == 1 ) islocal = 1;

    if ( islocal ) {

      // replace with a copy and maybe grid swizzle
      // return in;??
      double t = usecond();
      padded = in;
      tins += usecond() - t;
      
    } else {

      //////////////////////////////////////////////
      // Replace sequence with
      // ---------------------
      // (i) Gather high face(s); start comms
      // (ii) Gather low  face(s); start comms
      // (iii) Copy middle bit with localCopyRegion
      // (iv) Complete high face(s), insert slice(s)
      // (iv) Complete low  face(s), insert slice(s)
      //////////////////////////////////////////////
      // Middle bit
      double t = usecond();
      for(int x=0;x<local[dim];x++){
	InsertSliceLocal(in,padded,x,depth+x,dim);
      }
      tins += usecond() - t;
    
      // High bit
      t = usecond();
      shifted = cshift.Cshift(in,dim,depth);
      tshift += usecond() - t;

      t=usecond();
      for(int x=0;x<depth;x++){
	InsertSliceLocal(shifted,padded,local[dim]-depth+x,depth+local[dim]+x,dim);
      }
      tins += usecond() - t;
    
      // Low bit
      t = usecond();
      shifted = cshift.Cshift(in,dim,-depth);
      tshift += usecond() - t;
    
      t = usecond();
      for(int x=0;x<depth;x++){
	InsertSliceLocal(shifted,padded,x,x,dim);
      }
      tins += usecond() - t;

    }
    std::cout << GridLogPerformance << "PaddedCell::Expand timings: cshift:" << tshift/1000 << "ms, insert-slice:" << tins/1000 << "ms" << std::endl;
    
    return padded;
  }

  template<class vobj>
  inline Lattice<vobj> ExpandPeriodic(int dim, const Lattice<vobj> &in) const
  {
    Coordinate processors=unpadded_grid->_processors;
    GridBase *old_grid = in.Grid();
    GridCartesian *new_grid = grids[dim];//These are new grids
    Lattice<vobj>  padded(new_grid);
    //    Lattice<vobj> shifted(old_grid);    
    Coordinate local     =old_grid->LocalDimensions();
    Coordinate plocal    =new_grid->LocalDimensions();
    if(dim==0) conformable(old_grid,unpadded_grid);
    else       conformable(old_grid,grids[dim-1]);

    //    std::cout << " dim "<<dim<<" local "<<local << " padding to "<<plocal<<std::endl;
    double tins=0, tshift=0;

    int islocal = 0 ;
    if ( processors[dim] == 1 ) islocal = 1;

    if ( islocal ) {
      padded=in; // slightly different interface could avoid a copy operation
    } else {
      Face_exchange(in,padded,dim,depth);
      return padded;
    }
    return padded;
  }
  template<class vobj>
  void Face_exchange(const Lattice<vobj> &from,
		     Lattice<vobj> &to,
		     int dimension,int depth) const
  {
    typedef typename vobj::vector_type vector_type;
    typedef typename vobj::scalar_type scalar_type;
    typedef typename vobj::scalar_object sobj;

    RealD t_gather=0.0;
    RealD t_scatter=0.0;
    RealD t_comms=0.0;
    RealD t_copy=0.0;
    
    //    std::cout << GridLogMessage << "dimension " <<dimension<<std::endl;
    //    DumpSliceNorm(std::string("Face_exchange from"),from,dimension);
    GridBase *grid=from.Grid();
    GridBase *new_grid=to.Grid();

    Coordinate lds = from.Grid()->_ldimensions;
    Coordinate nlds=   to.Grid()->_ldimensions;
    Coordinate simd= from.Grid()->_simd_layout;
    int ld    = lds[dimension];
    int nld   = to.Grid()->_ldimensions[dimension];
    const int Nsimd = vobj::Nsimd();

    assert(depth<=lds[dimension]); // A must be on neighbouring node
    assert(depth>0);   // A caller bug if zero
    assert(ld+2*depth==nld);
    ////////////////////////////////////////////////////////////////////////////
    // Face size and byte calculations
    ////////////////////////////////////////////////////////////////////////////
    int buffer_size = 1;
    for(int d=0;d<lds.size();d++){
      if ( d!= dimension) buffer_size=buffer_size*lds[d];
    }
    buffer_size = buffer_size  / Nsimd;
    int rNsimd = Nsimd / simd[dimension];
    assert( buffer_size == from.Grid()->_slice_nblock[dimension]*from.Grid()->_slice_block[dimension] / simd[dimension]);

    static cshiftVector<vobj> send_buf; 
    static cshiftVector<vobj> recv_buf;
    send_buf.resize(buffer_size*2*depth);    
    recv_buf.resize(buffer_size*2*depth);

    std::vector<CommsRequest_t> fwd_req;   
    std::vector<CommsRequest_t> bwd_req;   

    int words = buffer_size;
    int bytes = words * sizeof(vobj);

    ////////////////////////////////////////////////////////////////////////////
    // Communication coords
    ////////////////////////////////////////////////////////////////////////////
    int comm_proc = 1;
    int xmit_to_rank;
    int recv_from_rank;
    grid->ShiftedRanks(dimension,comm_proc,xmit_to_rank,recv_from_rank);

    ////////////////////////////////////////////////////////////////////////////
    // Gather all surface terms up to depth "d"
    ////////////////////////////////////////////////////////////////////////////
    RealD t;
    RealD t_tot=-usecond();
    int plane=0;
    for ( int d=0;d < depth ; d ++ ) {
      int tag = d*1024 + dimension*2+0;

      t=usecond();
      GatherSlice(send_buf,from,d,dimension,plane*buffer_size); plane++;
      t_gather+=usecond()-t;

      t=usecond();
      grid->SendToRecvFromBegin(fwd_req,
				(void *)&send_buf[d*buffer_size], xmit_to_rank,
				(void *)&recv_buf[d*buffer_size], recv_from_rank, bytes, tag);
      t_comms+=usecond()-t;
     }
    for ( int d=0;d < depth ; d ++ ) {
      int tag = d*1024 + dimension*2+1;

      t=usecond();
      GatherSlice(send_buf,from,ld-depth+d,dimension,plane*buffer_size); plane++;
      t_gather+= usecond() - t;

      t=usecond();
      grid->SendToRecvFromBegin(bwd_req,
				(void *)&send_buf[(d+depth)*buffer_size], recv_from_rank,
				(void *)&recv_buf[(d+depth)*buffer_size], xmit_to_rank, bytes,tag);
      t_comms+=usecond()-t;
    }

    ////////////////////////////////////////////////////////////////////////////
    // Copy interior -- overlap this with comms
    ////////////////////////////////////////////////////////////////////////////
    int Nd = new_grid->Nd();
    Coordinate LL(Nd,0);
    Coordinate sz = grid->_ldimensions;
    Coordinate toLL(Nd,0);
    toLL[dimension]=depth;
    t=usecond();
    localCopyRegion(from,to,LL,toLL,sz);
    t_copy= usecond() - t;
    
    ////////////////////////////////////////////////////////////////////////////
    // Scatter all faces
    ////////////////////////////////////////////////////////////////////////////
    plane=0;

    t=usecond();
    grid->CommsComplete(fwd_req);
    t_comms+= usecond() - t;

    t=usecond();
    for ( int d=0;d < depth ; d ++ ) {
      ScatterSlice(recv_buf,to,nld-depth+d,dimension,plane*buffer_size); plane++;
    }
    t_scatter= usecond() - t;

    t=usecond();
    grid->CommsComplete(bwd_req);
    t_comms+= usecond() - t;
    
    t=usecond();
    for ( int d=0;d < depth ; d ++ ) {
      ScatterSlice(recv_buf,to,d,dimension,plane*buffer_size); plane++;
    }
    t_scatter+= usecond() - t;
    t_tot+=usecond();

    std::cout << GridLogPerformance << "PaddedCell::Expand new timings: gather :" << t_gather/1000  << "ms"<<std::endl;
    std::cout << GridLogPerformance << "PaddedCell::Expand new timings: scatter:" << t_scatter/1000   << "ms"<<std::endl;
    std::cout << GridLogPerformance << "PaddedCell::Expand new timings: copy   :" << t_copy/1000      << "ms"<<std::endl;
    std::cout << GridLogPerformance << "PaddedCell::Expand new timings: comms  :" << t_comms/1000     << "ms"<<std::endl;
    std::cout << GridLogPerformance << "PaddedCell::Expand new timings: total  :" << t_tot/1000     << "ms"<<std::endl;
    std::cout << GridLogPerformance << "PaddedCell::Expand new timings: gather :" << depth*4.0*bytes/t_gather << "MB/s"<<std::endl;
    std::cout << GridLogPerformance << "PaddedCell::Expand new timings: scatter:" << depth*4.0*bytes/t_scatter<< "MB/s"<<std::endl;
    std::cout << GridLogPerformance << "PaddedCell::Expand new timings: comms  :" << (RealD)4.0*bytes/t_comms   << "MB/s"<<std::endl;
    std::cout << GridLogPerformance << "PaddedCell::Expand new timings: face bytes  :" << depth*bytes/1e6 << "MB"<<std::endl;
  }
  
};
 

NAMESPACE_END(Grid);


