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

#define STENCIL_MAX (16)

#include <Grid/stencil/SimpleCompressor.h>   // subdir aggregate
#include <Grid/stencil/Lebesgue.h>   // subdir aggregate

//////////////////////////////////////////////////////////////////////////////////////////
// Must not lose sight that goal is to be able to construct really efficient
// gather to a point stencil code. CSHIFT is not the best way, so need
// additional stencil support.
//
// Stencil based code will exchange haloes and use a table lookup for neighbours.
// This will be done with generality to allow easier efficient implementations.
// Overlap of comms and compute is enabled by tabulating off-node connected,
//
// Generic services
// 0) Prebuild neighbour tables
// 1) Compute sizes of all haloes/comms buffers; allocate them.
// 2) Gather all faces, and communicate.
// 3) Loop over result sites, giving nbr index/offnode info for each
//
//////////////////////////////////////////////////////////////////////////////////////////

NAMESPACE_BEGIN(Grid);

// These can move into a params header and be given MacroMagic serialisation
struct DefaultImplParams {
  Coordinate dirichlet; // Blocksize of dirichlet BCs
  int  partialDirichlet;
  DefaultImplParams()  {
    dirichlet.resize(0);
    partialDirichlet=0;
  };
};

///////////////////////////////////////////////////////////////////
// Gather for when there *is* need to SIMD split with compression
///////////////////////////////////////////////////////////////////

void Gather_plane_table_compute (GridBase *grid,int dimension,int plane,int cbmask,
				 int off,std::vector<std::pair<int,int> > & table);

/*
template<class vobj,class cobj,class compressor>
void Gather_plane_simple_table (commVector<std::pair<int,int> >& table,const Lattice<vobj> &rhs,cobj *buffer,compressor &compress, int off,int so)   __attribute__((noinline));

template<class vobj,class cobj,class compressor>
void Gather_plane_simple_table (commVector<std::pair<int,int> >& table,const Lattice<vobj> &rhs,cobj *buffer,compressor &compress, int off,int so)
{
  int num=table.size();
  std::pair<int,int> *table_v = & table[0];

  auto rhs_v = rhs.View(AcceleratorRead);
  accelerator_forNB( i,num, vobj::Nsimd(), {
    compress.Compress(buffer[off+table_v[i].first],rhs_v[so+table_v[i].second]);
  });
  rhs_v.ViewClose();
}

///////////////////////////////////////////////////////////////////
// Gather for when there *is* need to SIMD split with compression
///////////////////////////////////////////////////////////////////
template<class cobj,class vobj,class compressor>
void Gather_plane_exchange_table(const Lattice<vobj> &rhs,
				 commVector<cobj *> pointers,
				 int dimension,int plane,
				 int cbmask,compressor &compress,int type) __attribute__((noinline));

template<class cobj,class vobj,class compressor>
void Gather_plane_exchange_table(commVector<std::pair<int,int> >& table,
				 const Lattice<vobj> &rhs,
				 std::vector<cobj *> &pointers,int dimension,int plane,int cbmask,
				 compressor &compress,int type)
{
  assert( (table.size()&0x1)==0);
  int num=table.size()/2;
  int so  = plane*rhs.Grid()->_ostride[dimension]; // base offset for start of plane

  auto rhs_v = rhs.View(AcceleratorRead);
  auto rhs_p = &rhs_v[0];
  auto p0=&pointers[0][0];
  auto p1=&pointers[1][0];
  auto tp=&table[0];
  accelerator_forNB(j, num, vobj::Nsimd(), {
      compress.CompressExchange(p0,p1, rhs_p, j,
				so+tp[2*j  ].second,
				so+tp[2*j+1].second,
				type);
  });
  rhs_v.ViewClose();
}
*/

void DslashResetCounts(void);
void DslashGetCounts(uint64_t &dirichlet,uint64_t &partial,uint64_t &full);
void DslashLogFull(void);
void DslashLogPartial(void);
void DslashLogDirichlet(void);

struct StencilEntry {
#ifdef GRID_CUDA
  uint64_t _byte_offset;       // 8 bytes
  uint32_t _offset;            // 4 bytes
#else
  uint64_t _byte_offset;       // 8 bytes
  uint64_t _offset;            // 8 bytes (8 ever required?)
#endif
  uint8_t _is_local;           // 1 bytes
  uint8_t _permute;            // 1 bytes
  uint8_t _around_the_world;   // 1 bytes
  uint8_t _pad;   // 1 bytes
};
// Could pack to 8 + 4 + 4 = 128 bit and use

template<class vobj,class cobj,class Parameters>
class CartesianStencilAccelerator {
 public:
  typedef AcceleratorVector<int,STENCIL_MAX> StencilVector;

  // Stencil runs along coordinate axes only; NO diagonal fill in.
  ////////////////////////////////////////
  // Basic Grid and stencil info
  ////////////////////////////////////////
  int           _checkerboard;
  int           _npoints; // Move to template param?
  int           _osites;
  StencilVector _directions;
  StencilVector _distances;
  ///////////////////////////////////////////////////
  // If true, this is FULLY communicated per face
  // Otherwise will either be full or partial dirichlet
  ///////////////////////////////////////////////////
  StencilVector _comms_send; 
  StencilVector _comms_recv; // this is FULLY communicated per face
  ///////////////////////////////////////////////////
  // If true, this is partially communicated per face
  ///////////////////////////////////////////////////
  StencilVector _comms_partial_send; 
  StencilVector _comms_partial_recv;
  //
  StencilVector _comm_buf_size;
  StencilVector _permute_type;
  StencilVector same_node;
  Coordinate    _simd_layout;
  Parameters    parameters;
  StencilEntry*  _entries_p;
  cobj* u_recv_buf_p;
  cobj* u_send_buf_p;

  accelerator_inline cobj *CommBuf(void) const { return u_recv_buf_p; }

  accelerator_inline int GetNodeLocal(int osite,int point) const {
    return this->_entries_p[point+this->_npoints*osite]._is_local;
  }
  accelerator_inline StencilEntry * GetEntry(int &ptype,int point,int osite) const {
    ptype = this->_permute_type[point]; return & this->_entries_p[point+this->_npoints*osite];
  }

  accelerator_inline uint64_t GetInfo(int &ptype,int &local,int &perm,int point,int ent,uint64_t base) const {
    uint64_t cbase = (uint64_t)&u_recv_buf_p[0];
    local = this->_entries_p[ent]._is_local;
    perm  = this->_entries_p[ent]._permute;
    if (perm)  ptype = this->_permute_type[point];
    if (local) {
      return  base + this->_entries_p[ent]._byte_offset;
    } else {
      return cbase + this->_entries_p[ent]._byte_offset;
    }
  }

  accelerator_inline uint64_t GetPFInfo(int ent,uint64_t base) const {
    uint64_t cbase = (uint64_t)&u_recv_buf_p[0];
    int local = this->_entries_p[ent]._is_local;
    if (local) return  base + this->_entries_p[ent]._byte_offset;
    else       return cbase + this->_entries_p[ent]._byte_offset;
  }

  accelerator_inline void iCoorFromIindex(Coordinate &coor,int lane) const
  {
    Lexicographic::CoorFromIndex(coor,lane,this->_simd_layout);
  }
};

template<class vobj,class cobj,class Parameters>
class CartesianStencilView : public CartesianStencilAccelerator<vobj,cobj,Parameters>
{
public:
  int *closed;
  StencilEntry *cpu_ptr;
  ViewMode      mode;
 public:
  // default copy constructor
  CartesianStencilView (const CartesianStencilView &refer_to_me) = default;

  CartesianStencilView (const CartesianStencilAccelerator<vobj,cobj,Parameters> &refer_to_me,ViewMode _mode)
    : CartesianStencilAccelerator<vobj,cobj,Parameters>(refer_to_me),
    cpu_ptr(this->_entries_p),
    mode(_mode)
  {
    this->_entries_p =(StencilEntry *)
      MemoryManager::ViewOpen(this->_entries_p,
			      this->_npoints*this->_osites*sizeof(StencilEntry),
			      mode,
			      AdviseDefault);
  }

  void ViewClose(void)
  {
    MemoryManager::ViewClose(this->cpu_ptr,this->mode);
  }

};

////////////////////////////////////////
// The Stencil Class itself
////////////////////////////////////////
template<class vobj,class cobj,class Parameters>
class CartesianStencil : public CartesianStencilAccelerator<vobj,cobj,Parameters> { // Stencil runs along coordinate axes only; NO diagonal fill in.
public:

  typedef typename cobj::vector_type vector_type;
  typedef typename cobj::scalar_object scalar_object;
  typedef const CartesianStencilView<vobj,cobj,Parameters> View_type;
  typedef typename View_type::StencilVector StencilVector;
  ///////////////////////////////////////////
  // Helper structs
  ///////////////////////////////////////////
  struct Packet {
    void * send_buf;
    void * recv_buf;
    Integer to_rank;
    Integer from_rank;
    Integer do_send;
    Integer do_recv;
    Integer xbytes;
    Integer rbytes;
  };
  struct Merge {
    static constexpr int Nsimd = vobj::Nsimd();
    cobj * mpointer;
    //    std::vector<scalar_object *> rpointers;
    std::vector<cobj *> vpointers;
    Integer buffer_size;
    Integer type;
    Integer partial; // partial dirichlet BCs
    Coordinate dims;
  };
  struct Decompress {
    static constexpr int Nsimd = vobj::Nsimd();
    cobj * kernel_p;
    cobj * mpi_p;
    Integer buffer_size;
    Integer partial; // partial dirichlet BCs
    Coordinate dims;
  };
  struct CopyReceiveBuffer {
    void * from_p;
    void * to_p;
    Integer bytes;
  };
  struct CachedTransfer {
    Integer direction;
    Integer OrthogPlane;
    Integer DestProc;
    Integer xbytes;
    Integer rbytes;
    Integer lane;
    Integer cb;
    void *recv_buf;
  };

protected:
  GridBase *                        _grid;
public:
  GridBase *Grid(void) const { return _grid; }
  LebesgueOrder *lo;

  ////////////////////////////////////////////////////////////////////////
  // Needed to conveniently communicate gparity parameters into GPU memory
  // without adding parameters. Perhaps a template parameter to StenciView is
  // required to pass general parameters.
  // Generalise as required later if needed
  ////////////////////////////////////////////////////////////////////////

  View_type View(ViewMode mode) const {
    View_type accessor(*( (View_type *) this),mode);
    return accessor;
  }

  int face_table_computed;
  int partialDirichlet;
  int fullDirichlet;
  std::vector<commVector<std::pair<int,int> > > face_table ;
  Vector<int> surface_list;

  stencilVector<StencilEntry>  _entries; // Resident in managed memory
  commVector<StencilEntry>     _entries_device; // Resident in managed memory
  std::vector<Packet> Packets;
  std::vector<Merge> Mergers;
  std::vector<Merge> MergersSHM;
  std::vector<Decompress> Decompressions;
  std::vector<Decompress> DecompressionsSHM;
  std::vector<CopyReceiveBuffer> CopyReceiveBuffers ;
  std::vector<CachedTransfer> CachedTransfers;
  std::vector<CommsRequest_t> MpiReqs;
  
  ///////////////////////////////////////////////////////////
  // Unified Comms buffers for all directions
  ///////////////////////////////////////////////////////////
  // Vectors that live on the symmetric heap in case of SHMEM
  // These are used; either SHM objects or refs to the above symmetric heap vectors
  // depending on comms target
  std::vector<cobj *> u_simd_send_buf;
  std::vector<cobj *> u_simd_recv_buf;

  int u_comm_offset;
  int _unified_buffer_size;

  ////////////////////////////////////////
  // Stencil query
  ////////////////////////////////////////
#if 1
  inline int SameNode(int point) {

    int dimension    = this->_directions[point];
    int displacement = this->_distances[point];

    int pd              = _grid->_processors[dimension];
    int fd              = _grid->_fdimensions[dimension];
    int ld              = _grid->_ldimensions[dimension];
    int rd              = _grid->_rdimensions[dimension];
    int simd_layout     = _grid->_simd_layout[dimension];
    int comm_dim        = _grid->_processors[dimension] >1 ;

    //    int recv_from_rank;
    //    int xmit_to_rank;

    if ( ! comm_dim ) return 1;
    if ( displacement == 0 ) return 1;
    return 0;
  }
#else
  // fancy calculation for shm code
  inline int SameNode(int point) {

    int dimension    = this->_directions[point];
    int displacement = this->_distances[point];

    int pd              = _grid->_processors[dimension];
    int fd              = _grid->_fdimensions[dimension];
    int ld              = _grid->_ldimensions[dimension];
    int rd              = _grid->_rdimensions[dimension];
    int simd_layout     = _grid->_simd_layout[dimension];
    int comm_dim        = _grid->_processors[dimension] >1 ;
 
    int recv_from_rank;
    int xmit_to_rank;

    if ( ! comm_dim ) return 1;

    int nbr_proc;
    if (displacement>0) nbr_proc = 1;
    else                nbr_proc = pd-1;

    // FIXME  this logic needs to be sorted for three link term
    //    assert( (displacement==1) || (displacement==-1));
    // Present hack only works for >= 4^4 subvol per node
    _grid->ShiftedRanks(dimension,nbr_proc,xmit_to_rank,recv_from_rank);

    void *shm = (void *) _grid->ShmBufferTranslate(recv_from_rank,this->u_recv_buf_p);

    if ( shm==NULL ) return 0;
    return 1;
  }
#endif
  //////////////////////////////////////////
  // Comms packet queue for asynch thread
  // Use OpenMP Tasks for cleaner ???
  // must be called *inside* parallel region
  //////////////////////////////////////////
  /*
  void CommunicateThreaded()
  {
#ifdef GRID_OMP
    int mythread = omp_get_thread_num();
    int nthreads = CartesianCommunicator::nCommThreads;
#else
    int mythread = 0;
    int nthreads = 1;
#endif
    if (nthreads == -1) nthreads = 1;
    if (mythread < nthreads) {
      for (int i = mythread; i < Packets.size(); i += nthreads) {
	uint64_t bytes = _grid->StencilSendToRecvFrom(Packets[i].send_buf,
						      Packets[i].to_rank,
						      Packets[i].recv_buf,
						      Packets[i].from_rank,
						      Packets[i].bytes,i);
      }
    }
  }
  */
  ////////////////////////////////////////////////////////////////////////
  // Non blocking send and receive. Necessarily parallel.
  ////////////////////////////////////////////////////////////////////////
  void CommunicateBegin(std::vector<std::vector<CommsRequest_t> > &reqs)
  {
    for(int i=0;i<Packets.size();i++){
      _grid->StencilSendToRecvFromBegin(MpiReqs,
					Packets[i].send_buf,
					Packets[i].to_rank,Packets[i].do_send,
					Packets[i].recv_buf,
					Packets[i].from_rank,Packets[i].do_recv,
					Packets[i].xbytes,Packets[i].rbytes,i);
    }
  }

  void CommunicateComplete(std::vector<std::vector<CommsRequest_t> > &reqs)
  {
    _grid->StencilSendToRecvFromComplete(MpiReqs,0);
    if   ( this->partialDirichlet ) DslashLogPartial();
    else if ( this->fullDirichlet ) DslashLogDirichlet();
    else DslashLogFull();
    acceleratorCopySynchronise();
    _grid->StencilBarrier(); 
  }
  ////////////////////////////////////////////////////////////////////////
  // Blocking send and receive. Either sequential or parallel.
  ////////////////////////////////////////////////////////////////////////
  void Communicate(void)
  {
    /////////////////////////////////////////////////////////
    // Concurrent and non-threaded asynch calls to MPI
    /////////////////////////////////////////////////////////
    std::vector<std::vector<CommsRequest_t> > reqs;
    this->CommunicateBegin(reqs);
    this->CommunicateComplete(reqs);
  }

  template<class compressor> void HaloExchange(const Lattice<vobj> &source,compressor &compress)
  {
    Prepare();
    HaloGather(source,compress);
    Communicate();
    CommsMergeSHM(compress);
    CommsMerge(compress);
  }

  template<class compressor> int HaloGatherDir(const Lattice<vobj> &source,compressor &compress,int point,int & face_idx)
  {
    int dimension    = this->_directions[point];
    int displacement = this->_distances[point];

    int fd = _grid->_fdimensions[dimension];
    int rd = _grid->_rdimensions[dimension];

    // Map to always positive shift modulo global full dimension.
    int shift = (displacement+fd)%fd;

    assert (source.Checkerboard()== this->_checkerboard);

    // the permute type
    int simd_layout     = _grid->_simd_layout[dimension];
    int comm_dim        = _grid->_processors[dimension] >1 ;
    int splice_dim      = _grid->_simd_layout[dimension]>1 && (comm_dim);

    int is_same_node = 1;
    // Gather phase
    int sshift [2];
    if ( comm_dim ) {
      sshift[0] = _grid->CheckerBoardShiftForCB(this->_checkerboard,dimension,shift,Even);
      sshift[1] = _grid->CheckerBoardShiftForCB(this->_checkerboard,dimension,shift,Odd);
      if ( sshift[0] == sshift[1] ) {
	if (splice_dim) {
	  auto tmp  = GatherSimd(source,dimension,shift,0x3,compress,face_idx,point);
	  is_same_node = is_same_node && tmp;
	} else {
	  auto tmp  = Gather(source,dimension,shift,0x3,compress,face_idx,point);
	  is_same_node = is_same_node && tmp;
	}
      } else {
	if(splice_dim){
	  // if checkerboard is unfavourable take two passes
	  // both with block stride loop iteration
	  auto tmp1 =  GatherSimd(source,dimension,shift,0x1,compress,face_idx,point);
	  auto tmp2 =  GatherSimd(source,dimension,shift,0x2,compress,face_idx,point);
	  is_same_node = is_same_node && tmp1 && tmp2;
	} else {
	  auto tmp1 = Gather(source,dimension,shift,0x1,compress,face_idx,point);
	  auto tmp2 = Gather(source,dimension,shift,0x2,compress,face_idx,point);
	  is_same_node = is_same_node && tmp1 && tmp2;
	}
      }
    }
    return is_same_node;
  }

  template<class compressor>
  void HaloGather(const Lattice<vobj> &source,compressor &compress)
  {
    _grid->StencilBarrier();// Synch shared memory on a single nodes

    assert(source.Grid()==_grid);

    u_comm_offset=0;

    // Gather all comms buffers
    int face_idx=0;
    for(int point = 0 ; point < this->_npoints; point++) {
      compress.Point(point);
      HaloGatherDir(source,compress,point,face_idx);
    }
    accelerator_barrier();
    face_table_computed=1;
    assert(u_comm_offset==_unified_buffer_size);

  }

  /////////////////////////
  // Implementation
  /////////////////////////
  void Prepare(void)
  {
    Decompressions.resize(0);
    DecompressionsSHM.resize(0);
    Mergers.resize(0);
    MergersSHM.resize(0);
    Packets.resize(0);
    CopyReceiveBuffers.resize(0);
    CachedTransfers.resize(0);
    MpiReqs.resize(0);
  }
  void AddCopy(void *from,void * to, Integer bytes)
  {
    CopyReceiveBuffer obj;
    obj.from_p = from;
    obj.to_p = to;
    obj.bytes= bytes;
    CopyReceiveBuffers.push_back(obj);
  }
  void CommsCopy()
  {
    //    These are device resident MPI buffers.
    for(int i=0;i<CopyReceiveBuffers.size();i++){
      cobj *from=(cobj *)CopyReceiveBuffers[i].from_p;
      cobj *to  =(cobj *)CopyReceiveBuffers[i].to_p;
      Integer words = CopyReceiveBuffers[i].bytes/sizeof(cobj);

      accelerator_forNB(j, words, cobj::Nsimd(), {
	  coalescedWrite(to[j] ,coalescedRead(from [j]));
      });
    }
  }
  
  Integer CheckForDuplicate(Integer direction, Integer OrthogPlane, Integer DestProc, void *recv_buf,Integer lane,
			    Integer xbytes,Integer rbytes,
			    Integer cb)
  {
    CachedTransfer obj;
    obj.direction   = direction;
    obj.OrthogPlane = OrthogPlane;
    obj.DestProc    = DestProc;
    obj.recv_buf    = recv_buf;
    obj.lane        = lane;
    obj.xbytes      = xbytes;
    obj.rbytes      = rbytes;
    obj.cb          = cb;

    for(int i=0;i<CachedTransfers.size();i++){
      if (   (CachedTransfers[i].direction  ==direction)
	   &&(CachedTransfers[i].OrthogPlane==OrthogPlane)
	   &&(CachedTransfers[i].DestProc   ==DestProc)
	   &&(CachedTransfers[i].xbytes      ==xbytes)
	   &&(CachedTransfers[i].rbytes      ==rbytes)
	   &&(CachedTransfers[i].lane       ==lane)
	   &&(CachedTransfers[i].cb         ==cb)
	     ){
	// FIXME worry about duplicate with partial compression
	// Wont happen as DWF has no duplicates, but...
	AddCopy(CachedTransfers[i].recv_buf,recv_buf,rbytes);
	return 1;
      }
    }

    CachedTransfers.push_back(obj);
    return 0;
  }
  void AddPacket(void *xmit,void * rcv,
		 Integer to, Integer do_send,
		 Integer from, Integer do_recv,
		 Integer xbytes,Integer rbytes){
    Packet p;
    p.send_buf = xmit;
    p.recv_buf = rcv;
    p.to_rank  = to;
    p.from_rank= from;
    p.do_send  = do_send;
    p.do_recv  = do_recv;
    p.xbytes    = xbytes;
    p.rbytes    = rbytes;
    //    if (do_send) std::cout << GridLogMessage << " MPI packet to   "<<to<< " of size "<<xbytes<<std::endl;
    //    if (do_recv) std::cout << GridLogMessage << " MPI packet from "<<from<< " of size "<<xbytes<<std::endl;
    Packets.push_back(p);
  }
  void AddDecompress(cobj *k_p,cobj *m_p,Integer buffer_size,std::vector<Decompress> &dv) {
    Decompress d;
    d.partial  = this->partialDirichlet;
    d.dims     = _grid->_fdimensions;
    d.kernel_p = k_p;
    d.mpi_p    = m_p;
    d.buffer_size = buffer_size;
    dv.push_back(d);
  }
  void AddMerge(cobj *merge_p,std::vector<cobj *> &rpointers,Integer buffer_size,Integer type,std::vector<Merge> &mv) {
    Merge m;
    m.partial  = this->partialDirichlet;
    m.dims     = _grid->_fdimensions;
    m.type     = type;
    m.mpointer = merge_p;
    m.vpointers= rpointers;
    m.buffer_size = buffer_size;
    mv.push_back(m);
  }
  template<class decompressor>  void CommsMerge(decompressor decompress)    {
    CommsCopy();
    CommsMerge(decompress,Mergers,Decompressions);
  }
  template<class decompressor>  void CommsMergeSHM(decompressor decompress) {
    assert(MergersSHM.size()==0);
    assert(DecompressionsSHM.size()==0);
  }

  template<class decompressor>
  void CommsMerge(decompressor decompress,std::vector<Merge> &mm,std::vector<Decompress> &dd)
  {
    for(int i=0;i<mm.size();i++){
      decompressor::MergeFace(decompress,mm[i]);
    }
    for(int i=0;i<dd.size();i++){
      decompressor::DecompressFace(decompress,dd[i]);
    }
  }
  ////////////////////////////////////////
  // Set up routines
  ////////////////////////////////////////
  void PrecomputeByteOffsets(void){
    for(int i=0;i<_entries.size();i++){
      if( _entries[i]._is_local ) {
	_entries[i]._byte_offset = _entries[i]._offset*sizeof(vobj);
      } else {
	_entries[i]._byte_offset = _entries[i]._offset*sizeof(cobj);
      }
    }
  };

  // Move interior/exterior split into the generic stencil
  // FIXME Explicit Ls in interface is a pain. Should just use a vol
  void BuildSurfaceList(int Ls,int vol4){

    // find same node for SHM
    // Here we know the distance is 1 for WilsonStencil
    for(int point=0;point<this->_npoints;point++){
      this->same_node[point] = this->SameNode(point);
    }

    for(int site = 0 ;site< vol4;site++){
      int local = 1;
      for(int point=0;point<this->_npoints;point++){
	if( (!this->GetNodeLocal(site*Ls,point)) && (!this->same_node[point]) ){
	  local = 0;
	}
      }
      if(local == 0) {
	for(int s=0;s<Ls;s++){
	  surface_list.push_back(site*Ls+s);
	}
      }
    }
    std::cout << GridLogDebug << "BuildSurfaceList size is "<<surface_list.size()<<std::endl;
  }
  /// Introduce a block structure and switch off comms on boundaries
  void DirichletBlock(const Coordinate &dirichlet_block)
  {
    for(int ii=0;ii<this->_npoints;ii++){
      int dimension    = this->_directions[ii];
      int displacement = this->_distances[ii];
      int gd = _grid->_gdimensions[dimension];
      int fd = _grid->_fdimensions[dimension];
      int pd = _grid->_processors [dimension];
      int pc = _grid->_processor_coor[dimension];
      int ld = fd/pd;
      ///////////////////////////////////////////
      // Figure out dirichlet send and receive
      // on this leg of stencil.
      ///////////////////////////////////////////
      int comm_dim        = _grid->_processors[dimension] >1 ;
      int block = dirichlet_block[dimension];
      this->_comms_send[ii] = comm_dim;
      this->_comms_recv[ii] = comm_dim;
      this->_comms_partial_send[ii] = 0;
      this->_comms_partial_recv[ii] = 0;
      if ( block && comm_dim ) {
	assert(abs(displacement) < ld );
	// Quiesce communication across block boundaries
	if( displacement > 0 ) {
	  // High side, low side
	  // | <--B--->|
	  // |    |    |
	  //           noR
	  // noS
	  if ( ( (ld*(pc+1) ) % block ) == 0 ) this->_comms_recv[ii] = 0;
	  if ( ( (ld*pc     ) % block ) == 0 ) this->_comms_send[ii] = 0;
	} else {
	  // High side, low side
	  // | <--B--->|
	  // |    |    |
	  //           noS
	  // noR
	  if ( ( (ld*(pc+1) ) % block ) == 0 ) this->_comms_send[ii] = 0;
	  if ( ( (ld*pc     ) % block ) == 0 ) this->_comms_recv[ii] = 0;
	}
	if ( partialDirichlet ) {
	  this->_comms_partial_send[ii] = !this->_comms_send[ii];
	  this->_comms_partial_recv[ii] = !this->_comms_recv[ii];
	}
      }
    }
  }
  CartesianStencil(GridBase *grid,
		   int npoints,
		   int checkerboard,
		   const std::vector<int> &directions,
		   const std::vector<int> &distances,
		   Parameters p=Parameters())
  {
    face_table_computed=0;
    _grid    = grid;
    this->parameters=p;
    /////////////////////////////////////
    // Initialise the base
    /////////////////////////////////////
    this->_npoints = npoints;
    this->_comm_buf_size.resize(npoints),
    this->_permute_type.resize(npoints),
    this->_simd_layout = _grid->_simd_layout; // copy simd_layout to give access to Accelerator Kernels
    this->_directions = StencilVector(directions);
    this->_distances  = StencilVector(distances);
    this->_comms_send.resize(npoints); 
    this->_comms_recv.resize(npoints); 
    this->same_node.resize(npoints);

    if ( p.dirichlet.size() ==0 ) p.dirichlet.resize(grid->Nd(),0);
    partialDirichlet = p.partialDirichlet;
    DirichletBlock(p.dirichlet); // comms send/recv set up
    fullDirichlet=0;
    for(int d=0;d<p.dirichlet.size();d++){
      if (p.dirichlet[d]) fullDirichlet=1;
    }

    _unified_buffer_size=0;
    surface_list.resize(0);

    this->_osites  = _grid->oSites();

    _entries.resize(this->_npoints* this->_osites);
    this->_entries_p = &_entries[0];
    for(int ii=0;ii<npoints;ii++){

      int i = ii; // reverse direction to get SIMD comms done first
      int point = i;

      int dimension    = directions[i];
      int displacement = distances[i];
      int shift = displacement;

      int gd = _grid->_gdimensions[dimension];
      int fd = _grid->_fdimensions[dimension];
      int pd = _grid->_processors [dimension];
      //      int ld = gd/pd;
      int rd = _grid->_rdimensions[dimension];
      int pc = _grid->_processor_coor[dimension];
      this->_permute_type[point]=_grid->PermuteType(dimension);

      this->_checkerboard = checkerboard;

      int simd_layout     = _grid->_simd_layout[dimension];
      int comm_dim        = _grid->_processors[dimension] >1 ;
      int splice_dim      = _grid->_simd_layout[dimension]>1 && (comm_dim);
      int rotate_dim      = _grid->_simd_layout[dimension]>2;

      assert ( (rotate_dim && comm_dim) == false) ; // Do not think spread out is supported

      int sshift[2];
      //////////////////////////
      // Underlying approach. For each local site build
      // up a table containing the npoint "neighbours" and whether they
      // live in lattice or a comms buffer.
      //////////////////////////
      if ( !comm_dim ) {
	sshift[0] = _grid->CheckerBoardShiftForCB(this->_checkerboard,dimension,shift,Even);
	sshift[1] = _grid->CheckerBoardShiftForCB(this->_checkerboard,dimension,shift,Odd);

	if ( sshift[0] == sshift[1] ) {
	  Local(point,dimension,shift,0x3);
	} else {
	  Local(point,dimension,shift,0x1);// if checkerboard is unfavourable take two passes
	  Local(point,dimension,shift,0x2);// both with block stride loop iteration
	}
      } else {
	// All permute extract done in comms phase prior to Stencil application
	//        So tables are the same whether comm_dim or splice_dim
	sshift[0] = _grid->CheckerBoardShiftForCB(this->_checkerboard,dimension,shift,Even);
	sshift[1] = _grid->CheckerBoardShiftForCB(this->_checkerboard,dimension,shift,Odd);
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

    int maxl=2;
    u_simd_send_buf.resize(maxl);
    u_simd_recv_buf.resize(maxl);
    this->u_send_buf_p=(cobj *)_grid->ShmBufferMalloc(_unified_buffer_size*sizeof(cobj));
    this->u_recv_buf_p=(cobj *)_grid->ShmBufferMalloc(_unified_buffer_size*sizeof(cobj));

    for(int l=0;l<maxl;l++){
      u_simd_recv_buf[l] = (cobj *)_grid->ShmBufferMalloc(_unified_buffer_size*sizeof(cobj));
      u_simd_send_buf[l] = (cobj *)_grid->ShmBufferMalloc(_unified_buffer_size*sizeof(cobj));
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

      //      int o   = 0;
      int bo  = x * _grid->_ostride[dimension];

      int cb= (cbmask==0x2)? Odd : Even;

      int sshift = _grid->CheckerBoardShiftForCB(this->_checkerboard,dimension,shift,cb);
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
	int wrap = sshift/rd; wrap=wrap % ly; // but it is local anyway
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

    int comms_recv      = this->_comms_recv[point] || this->_comms_partial_recv[point] ;
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

    // done in reduced dims, so SIMD factored
    int buffer_size = _grid->_slice_nblock[dimension]*_grid->_slice_block[dimension];

    this->_comm_buf_size[point] = buffer_size; // Size of _one_ plane. Multiple planes may be gathered and

    // send to one or more remote nodes.

    int cb= (cbmask==0x2)? Odd : Even;
    int sshift= _grid->CheckerBoardShiftForCB(this->_checkerboard,dimension,shift,cb);

    for(int x=0;x<rd;x++){

      int permute_type=grid->PermuteType(dimension);
      int permute_slice;

      int sx        =  (x+sshift)%rd;

      int offnode = 0;
      if ( simd_layout > 1 ) {

	permute_slice=1;
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
	permute_slice=0;
      }

      int wraparound=0;
      if ( (shiftpm==-1) && (sx>x) && (grid->_processor_coor[dimension]==0) ) {
	wraparound = 1;
      }
      if ( (shiftpm== 1) && (sx<x) && (grid->_processor_coor[dimension]==grid->_processors[dimension]-1) ) {
	wraparound = 1;
      }

      // Wrap locally dirichlet support case OR node local
      if ( offnode==0 ) {

	permute_slice=0;
	CopyPlane(point,dimension,x,sx,cbmask,permute_slice,wraparound);
	
      } else {

	if ( comms_recv ) {

	  ScatterPlane(point,dimension,x,cbmask,_unified_buffer_size,wraparound); // permute/extract/merge is done in comms phase

	} else { 

	  CopyPlane(point,dimension,x,sx,cbmask,permute_slice,wraparound);

	}

      }
      
      if ( offnode ) {
	int words = buffer_size;
	if (cbmask != 0x3) words=words>>1;
	_unified_buffer_size    += words;
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
	  int idx=point+(lo+o+b)*this->_npoints;
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
	    int idx = point+(lo+o+b)*this->_npoints;
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
	  int idx=point+(so+o+b)*this->_npoints;
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
	    int idx = point+(so+o+b)*this->_npoints;
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

  template<class compressor>
  int Gather(const Lattice<vobj> &rhs,int dimension,int shift,int cbmask,compressor & compress,int &face_idx, int point)
  {
    typedef typename cobj::vector_type vector_type;

    int comms_send   = this->_comms_send[point];
    int comms_recv   = this->_comms_recv[point];
    int comms_partial_send   = this->_comms_partial_send[point] ;
    int comms_partial_recv   = this->_comms_partial_recv[point] ;
    
    assert(rhs.Grid()==_grid);
    //	  conformable(_grid,rhs.Grid());

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
    int sshift= _grid->CheckerBoardShiftForCB(rhs.Checkerboard(),dimension,shift,cb);

    for(int x=0;x<rd;x++){

      int sx        = (x+sshift)%rd;
      int comm_proc = ((x+sshift)/rd)%pd;
      
      if (comm_proc) {
	
	int words = buffer_size;
	if (cbmask != 0x3) words=words>>1;

	int bytes =  words * compress.CommDatumSize();
	int xbytes;
	int rbytes;

	if ( comms_send ) xbytes = bytes; // Full send
	else if ( comms_partial_send ) xbytes = bytes/compressor::PartialCompressionFactor(_grid);
	else xbytes = 0; // full dirichlet

	if ( comms_recv ) rbytes = bytes;
	else if ( comms_partial_recv ) rbytes = bytes/compressor::PartialCompressionFactor(_grid);
	else rbytes = 0;
	
	int so  = sx*rhs.Grid()->_ostride[dimension]; // base offset for start of plane
	int comm_off = u_comm_offset;

	int recv_from_rank;
	int xmit_to_rank;
	cobj *recv_buf;
	cobj *send_buf;
	_grid->ShiftedRanks(dimension,comm_proc,xmit_to_rank,recv_from_rank);

	assert (xmit_to_rank   != _grid->ThisRank());
	assert (recv_from_rank != _grid->ThisRank());

	if ( !face_table_computed ) {
	  face_table.resize(face_idx+1);
	  std::vector<std::pair<int,int> >  face_table_host ;
	  Gather_plane_table_compute ((GridBase *)_grid,dimension,sx,cbmask,comm_off,face_table_host);
	  //	  std::cout << "bytes expect "<< bytes << " " << face_table_host.size()* compress.CommDatumSize()<<std::endl;
	    face_table[face_idx].resize(face_table_host.size());
	    acceleratorCopyToDevice(&face_table_host[0],
				    &face_table[face_idx][0],
				    face_table[face_idx].size()*sizeof(face_table_host[0]));
	}


	if ( (compress.DecompressionStep()&&comms_recv) || comms_partial_recv ) {
	  recv_buf=u_simd_recv_buf[0];
	} else {
	  recv_buf=this->u_recv_buf_p;
	}

	// potential SHM fast path for intranode
	int shm_send=0;
	int shm_recv=0;
#ifdef SHM_FAST_PATH
	// Put directly in place if we can
	send_buf = (cobj *)_grid->ShmBufferTranslate(xmit_to_rank,recv_buf);
	if ( (send_buf==NULL) ) {
	  shm_send=0;
	  send_buf = this->u_send_buf_p;
	} else {
	  shm_send=1;
	}
	void *test_ptr = _grid->ShmBufferTranslate(recv_from_rank,recv_buf);
	if ( test_ptr != NULL ) shm_recv = 1;
	//	static int printed;
	//	if (!printed){
	  //	  std::cout << " GATHER FAST PATH SHM "<<shm_send<< " "<<shm_recv<<std::endl;
	//	  printed = 1;
	//	}
#else
	////////////////////////////////////////////////////////
	// Gather locally
	////////////////////////////////////////////////////////
	send_buf = this->u_send_buf_p; // Gather locally, must send
	assert(send_buf!=NULL);
#endif

	//	std::cout << " GatherPlaneSimple partial send "<< comms_partial_send<<std::endl;
	compressor::Gather_plane_simple(face_table[face_idx],rhs,send_buf,compress,comm_off,so,comms_partial_send);

        int duplicate = CheckForDuplicate(dimension,sx,comm_proc,(void *)&recv_buf[comm_off],0,xbytes,rbytes,cbmask);
	if ( !duplicate ) { // Force comms for now
	  
	  ///////////////////////////////////////////////////////////
	  // Build a list of things to do after we synchronise GPUs
	  // Start comms now???
	  ///////////////////////////////////////////////////////////
	  int do_send = (comms_send|comms_partial_send) && (!shm_send );
	  int do_recv = (comms_send|comms_partial_send) && (!shm_recv );
	  
	  AddPacket((void *)&send_buf[comm_off],
		    (void *)&recv_buf[comm_off],
		    xmit_to_rank, do_send,
		    recv_from_rank, do_recv,
		    xbytes,rbytes);
	}

	if ( (compress.DecompressionStep() && comms_recv) || comms_partial_recv ) {
	  AddDecompress(&this->u_recv_buf_p[comm_off],
			&recv_buf[comm_off],
			words,Decompressions);
	}
	
	u_comm_offset+=words;
	face_idx++;
      }
    }
    return 0;
  }

  template<class compressor>
  int  GatherSimd(const Lattice<vobj> &rhs,int dimension,int shift,int cbmask,compressor &compress,int & face_idx,int point)
  {
    const int Nsimd = _grid->Nsimd();

    const int maxl =2;// max layout in a direction

    int comms_send   = this->_comms_send[point];
    int comms_recv   = this->_comms_recv[point];
    int comms_partial_send   = this->_comms_partial_send[point] ;
    int comms_partial_recv   = this->_comms_partial_recv[point] ;

    int fd = _grid->_fdimensions[dimension];
    int rd = _grid->_rdimensions[dimension];
    int ld = _grid->_ldimensions[dimension];
    int pd              = _grid->_processors[dimension];
    int simd_layout     = _grid->_simd_layout[dimension];
    int comm_dim        = _grid->_processors[dimension] >1 ;
    assert(comm_dim==1);
    // This will not work with a rotate dim
    assert(simd_layout==maxl);
    assert(shift>=0);
    assert(shift<fd);


    int permute_type=_grid->PermuteType(dimension);

    ///////////////////////////////////////////////
    // Simd direction uses an extract/merge pair
    ///////////////////////////////////////////////
    int buffer_size = _grid->_slice_nblock[dimension]*_grid->_slice_block[dimension];
    //    int words = sizeof(cobj)/sizeof(vector_type);

    assert(cbmask==0x3); // Fixme think there is a latent bug if not true
                         // This assert will trap it if ever hit. Not hit normally so far
    int reduced_buffer_size = buffer_size;
    if (cbmask != 0x3) reduced_buffer_size=buffer_size>>1;

    int datum_bytes = compress.CommDatumSize();
    int bytes = (reduced_buffer_size*datum_bytes)/simd_layout;

    // how many bytes on wire : partial dirichlet or dirichlet may set to < bytes
    int xbytes; 
    int rbytes; 
    
    assert(bytes*simd_layout == reduced_buffer_size*datum_bytes);

    std::vector<cobj *> rpointers(maxl);
    std::vector<cobj *> spointers(maxl);

    ///////////////////////////////////////////
    // Work out what to send where
    ///////////////////////////////////////////

    int cb    = (cbmask==0x2)? Odd : Even;
    int sshift= _grid->CheckerBoardShiftForCB(rhs.Checkerboard(),dimension,shift,cb);

    // loop over outer coord planes orthog to dim
    for(int x=0;x<rd;x++){

      int any_offnode = ( ((x+sshift)%fd) >= rd );

      if ( any_offnode ) {

	int comm_off = u_comm_offset;
	for(int i=0;i<maxl;i++){
	  spointers[i] = (cobj *) &u_simd_send_buf[i][comm_off];
	}

	int sx   = (x+sshift)%rd;

	if ( !face_table_computed ) {
	  face_table.resize(face_idx+1);
	  std::vector<std::pair<int,int> >  face_table_host ;

	  Gather_plane_table_compute ((GridBase *)_grid,dimension,sx,cbmask,comm_off,face_table_host);
	  face_table[face_idx].resize(face_table_host.size());
	  acceleratorCopyToDevice(&face_table_host[0],
				  &face_table[face_idx][0],
				  face_table[face_idx].size()*sizeof(face_table_host[0]));

	}

	
	if ( comms_send ) xbytes = bytes;
	else if ( comms_partial_send ) xbytes = bytes/compressor::PartialCompressionFactor(_grid);
	else xbytes = 0;

	if ( comms_recv ) rbytes = bytes;
	else if ( comms_partial_recv ) rbytes = bytes/compressor::PartialCompressionFactor(_grid);
	else rbytes = 0;

	// Gathers SIMD lanes for send and merge
	// Different faces can be full comms or partial comms with  multiple ranks per node
	if ( comms_send || comms_recv||comms_partial_send||comms_partial_recv ) {

	  int partial = partialDirichlet;
	  compressor::Gather_plane_exchange(face_table[face_idx],rhs,
					    spointers,dimension,sx,cbmask,
					    compress,permute_type,partial );
	}
	face_idx++;

	//spointers[0] -- low simd coor
	//spointers[1] -- high simd coor
	for(int i=0;i<maxl;i++){

	  int my_coor  = rd*i + x;            // self explanatory
	  int nbr_coor = my_coor+sshift;      // self explanatory

	  int nbr_proc = ((nbr_coor)/ld) % pd;// relative shift in processors
	  int nbr_lcoor= (nbr_coor%ld);       // local plane coor on neighbour node
	  int nbr_ic   = (nbr_lcoor)/rd;      // inner coord of peer simd lane "i"
	  int nbr_ox   = (nbr_lcoor%rd);      // outer coord of peer "x"

	  int nbr_plane = nbr_ic;
	  assert (sx == nbr_ox);

	  auto rp = &u_simd_recv_buf[i        ][comm_off];
	  auto sp = &u_simd_send_buf[nbr_plane][comm_off];

	  if(nbr_proc){

	    int recv_from_rank;
	    int xmit_to_rank;
	    int shm_send=0;

	    _grid->ShiftedRanks(dimension,nbr_proc,xmit_to_rank,recv_from_rank);
#ifdef SHM_FAST_PATH
  #warning STENCIL SHM FAST PATH SELECTED
  	  int shm_recv=0;
	    // shm == receive pointer         if offnode
	    // shm == Translate[send pointer] if on node -- my view of his send pointer
	    cobj *shm = (cobj *) _grid->ShmBufferTranslate(recv_from_rank,sp);
	    if (shm==NULL) {
	      shm = rp;
	      // we found a packet that comes from MPI and contributes to this shift.
	      // is_same_node is only used in the WilsonStencil, and gets set for this point in the stencil.
	      // Kernel will add the exterior_terms except if is_same_node.
	      // leg of stencil
	      shm_recv=0;
	    } else {
	      shm_recv=1;
	    }
	    rpointers[i] = shm;
	    // Test send side
	    void *test_ptr = (void *) _grid->ShmBufferTranslate(xmit_to_rank,sp);
	    if ( test_ptr != NULL ) shm_send = 1;
	    //	    static int printed;
	    //	    if (!printed){
	    //	      std::cout << " GATHERSIMD FAST PATH SHM "<<shm_send<< " "<<shm_recv<<std::endl;
	    //	      printed = 1;
	    //	    }
#else
	    rpointers[i] = rp;
#endif
	    
	    int duplicate = CheckForDuplicate(dimension,sx,nbr_proc,(void *)rp,i,xbytes,rbytes,cbmask);
	    if ( !duplicate  ) { 
	      if ( (bytes != rbytes) && (rbytes!=0) ){
		acceleratorMemSet(rp,0,bytes); // Zero prefill comms buffer to zero
	      }
	      int do_send = (comms_send|comms_partial_send) && (!shm_send );
	      AddPacket((void *)sp,(void *)rp,
			xmit_to_rank,do_send,
			recv_from_rank,do_send,
			xbytes,rbytes);
	    }

	  } else {

	    rpointers[i] = sp;

	  }
	}
	// rpointer may be doing a remote read in the gather over SHM
	if ( comms_recv|comms_partial_recv ) {
	  AddMerge(&this->u_recv_buf_p[comm_off],rpointers,reduced_buffer_size,permute_type,Mergers);
	}

	u_comm_offset     +=buffer_size;

      }
    }
    return 0;
  }

};
NAMESPACE_END(Grid);

#endif
