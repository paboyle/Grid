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

///////////////////////////////////////////////////////////////////
// Gather for when there *is* need to SIMD split with compression
///////////////////////////////////////////////////////////////////

void Gather_plane_table_compute (GridBase *grid,int dimension,int plane,int cbmask,
				 int off,Vector<std::pair<int,int> > & table);

template<class vobj,class cobj,class compressor>
void Gather_plane_simple_table (Vector<std::pair<int,int> >& table,const Lattice<vobj> &rhs,cobj *buffer,compressor &compress, int off,int so)   __attribute__((noinline));

template<class vobj,class cobj,class compressor>
void Gather_plane_simple_table (Vector<std::pair<int,int> >& table,const Lattice<vobj> &rhs,cobj *buffer,compressor &compress, int off,int so)
{
  int num=table.size();
  std::pair<int,int> *table_v = & table[0];

  auto rhs_v = rhs.View(AcceleratorRead);
  accelerator_forNB( i,num, vobj::Nsimd(), {
    typedef decltype(coalescedRead(buffer[0])) compressed_t;
    compressed_t   tmp_c;
    uint64_t o = table_v[i].first;
    compress.Compress(&tmp_c,0,rhs_v(so+table_v[i].second));
    coalescedWrite(buffer[off+o],tmp_c);
  });
  rhs_v.ViewClose();
// Further optimisatoin: i) software prefetch the first element of the next table entry, prefetch the table
}

///////////////////////////////////////////////////////////////////
// Gather for when there *is* need to SIMD split with compression
///////////////////////////////////////////////////////////////////
template<class cobj,class vobj,class compressor>
void Gather_plane_exchange_table(const Lattice<vobj> &rhs,
				 Vector<cobj *> pointers,int dimension,int plane,int cbmask,compressor &compress,int type) __attribute__((noinline));

template<class cobj,class vobj,class compressor>
void Gather_plane_exchange_table(Vector<std::pair<int,int> >& table,const Lattice<vobj> &rhs,
				 Vector<cobj *> pointers,int dimension,int plane,int cbmask,
				 compressor &compress,int type)
{
  assert( (table.size()&0x1)==0);
  int num=table.size()/2;
  int so  = plane*rhs.Grid()->_ostride[dimension]; // base offset for start of plane

  auto rhs_v = rhs.View(AcceleratorRead);
  auto p0=&pointers[0][0];
  auto p1=&pointers[1][0];
  auto tp=&table[0];
  accelerator_forNB(j, num, 1, {
      compress.CompressExchange(p0,p1, &rhs_v[0], j,
			      so+tp[2*j  ].second,
			      so+tp[2*j+1].second,
			      type);
  });
  rhs_v.ViewClose();
}

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
 private:
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
  typedef typename cobj::scalar_type scalar_type;
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
    Integer bytes;
  };
  struct Merge {
    cobj * mpointer;
    Vector<scalar_object *> rpointers;
    Vector<cobj *> vpointers;
    Integer buffer_size;
    Integer type;
  };
  struct Decompress {
    cobj * kernel_p;
    cobj * mpi_p;
    Integer buffer_size;
  };


protected:
  GridBase *                        _grid;

public:
  GridBase *Grid(void) const { return _grid; }

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
  std::vector<Vector<std::pair<int,int> > > face_table ;
  Vector<int> surface_list;

  stencilVector<StencilEntry>  _entries; // Resident in managed memory
  std::vector<Packet> Packets;
  std::vector<Merge> Mergers;
  std::vector<Merge> MergersSHM;
  std::vector<Decompress> Decompressions;
  std::vector<Decompress> DecompressionsSHM;

  ///////////////////////////////////////////////////////////
  // Unified Comms buffers for all directions
  ///////////////////////////////////////////////////////////
  // Vectors that live on the symmetric heap in case of SHMEM
  // These are used; either SHM objects or refs to the above symmetric heap vectors
  // depending on comms target
  Vector<cobj *> u_simd_send_buf;
  Vector<cobj *> u_simd_recv_buf;

  int u_comm_offset;
  int _unified_buffer_size;

  /////////////////////////////////////////
  // Timing info; ugly; possibly temporary
  /////////////////////////////////////////
  double commtime;
  double mpi3synctime;
  double mpi3synctime_g;
  double shmmergetime;
  double gathertime;
  double gathermtime;
  double halogtime;
  double mergetime;
  double decompresstime;
  double comms_bytes;
  double shm_bytes;
  double splicetime;
  double nosplicetime;
  double calls;
  std::vector<double> comm_bytes_thr;
  std::vector<double> shm_bytes_thr;
  std::vector<double> comm_time_thr;
  std::vector<double> comm_enter_thr;
  std::vector<double> comm_leave_thr;

  ////////////////////////////////////////
  // Stencil query
  ////////////////////////////////////////
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
    else                 nbr_proc = pd-1;

    // FIXME  this logic needs to be sorted for three link term
    //    assert( (displacement==1) || (displacement==-1));
    // Present hack only works for >= 4^4 subvol per node
    _grid->ShiftedRanks(dimension,nbr_proc,xmit_to_rank,recv_from_rank);

    void *shm = (void *) _grid->ShmBufferTranslate(recv_from_rank,this->u_recv_buf_p);

    if ( shm==NULL ) return 0;

    return 1;
  }

  //////////////////////////////////////////
  // Comms packet queue for asynch thread
  // Use OpenMP Tasks for cleaner ???
  //////////////////////////////////////////
  void CommunicateThreaded()
  {
#ifdef GRID_OMP
    // must be called in parallel region
    int mythread = omp_get_thread_num();
    int nthreads = CartesianCommunicator::nCommThreads;
#else
    int mythread = 0;
    int nthreads = 1;
#endif
    if (nthreads == -1) nthreads = 1;
    if (mythread < nthreads) {
      comm_enter_thr[mythread] = usecond();
      for (int i = mythread; i < Packets.size(); i += nthreads) {
	uint64_t bytes = _grid->StencilSendToRecvFrom(Packets[i].send_buf,
						      Packets[i].to_rank,
						      Packets[i].recv_buf,
						      Packets[i].from_rank,
						      Packets[i].bytes,i);
	comm_bytes_thr[mythread] += bytes;
	shm_bytes_thr[mythread] += 2*Packets[i].bytes-bytes; // Send + Recv.

      }
      comm_leave_thr[mythread]= usecond();
      comm_time_thr[mythread] += comm_leave_thr[mythread] - comm_enter_thr[mythread];
    }
  }

  void CollateThreads(void)
  {
    int nthreads = CartesianCommunicator::nCommThreads;
    double first=0.0;
    double last =0.0;

    for(int t=0;t<nthreads;t++) {

      double t0 = comm_enter_thr[t];
      double t1 = comm_leave_thr[t];
      comms_bytes+=comm_bytes_thr[t];
      shm_bytes  +=shm_bytes_thr[t];

      comm_enter_thr[t] = 0.0;
      comm_leave_thr[t] = 0.0;
      comm_time_thr[t]   = 0.0;
      comm_bytes_thr[t]=0;
      shm_bytes_thr[t]=0;

      if ( first == 0.0 ) first = t0;                   // first is t0
      if ( (t0 > 0.0) && ( t0 < first ) ) first = t0;   // min time seen

      if ( t1 > last ) last = t1;                       // max time seen

    }
    commtime+= last-first;
  }
  ////////////////////////////////////////////////////////////////////////
  // Non blocking send and receive. Necessarily parallel.
  ////////////////////////////////////////////////////////////////////////
  void CommunicateBegin(std::vector<std::vector<CommsRequest_t> > &reqs)
  {
    reqs.resize(Packets.size());
    commtime-=usecond();
    for(int i=0;i<Packets.size();i++){
      uint64_t bytes=_grid->StencilSendToRecvFromBegin(reqs[i],
						     Packets[i].send_buf,
						     Packets[i].to_rank,
						     Packets[i].recv_buf,
						     Packets[i].from_rank,
						     Packets[i].bytes,i);
      comms_bytes+=bytes;
      shm_bytes  +=2*Packets[i].bytes-bytes;
    }
  }

  void CommunicateComplete(std::vector<std::vector<CommsRequest_t> > &reqs)
  {
    for(int i=0;i<Packets.size();i++){
      _grid->StencilSendToRecvFromComplete(reqs[i],i);
    }
    commtime+=usecond();
  }
  ////////////////////////////////////////////////////////////////////////
  // Blocking send and receive. Either sequential or parallel.
  ////////////////////////////////////////////////////////////////////////
  void Communicate(void)
  {
    if ( CartesianCommunicator::CommunicatorPolicy == CartesianCommunicator::CommunicatorPolicySequential ){
      thread_region {
	// must be called in parallel region
	int mythread  = thread_num();
	int maxthreads= thread_max();
	int nthreads = CartesianCommunicator::nCommThreads;
	assert(nthreads <= maxthreads);
	if (nthreads == -1) nthreads = 1;
	if (mythread < nthreads) {
	  for (int i = mythread; i < Packets.size(); i += nthreads) {
	    double start = usecond();
	    uint64_t bytes= _grid->StencilSendToRecvFrom(Packets[i].send_buf,
							 Packets[i].to_rank,
							 Packets[i].recv_buf,
							 Packets[i].from_rank,
							 Packets[i].bytes,i);
	    comm_bytes_thr[mythread] += bytes;
	    shm_bytes_thr[mythread]  += Packets[i].bytes - bytes;
	    comm_time_thr[mythread]  += usecond() - start;
	  }
	}
      }
    } else { // Concurrent and non-threaded asynch calls to MPI
      std::vector<std::vector<CommsRequest_t> > reqs;
      this->CommunicateBegin(reqs);
      this->CommunicateComplete(reqs);
    }
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
	  splicetime-=usecond();
	  auto tmp  = GatherSimd(source,dimension,shift,0x3,compress,face_idx);
	  is_same_node = is_same_node && tmp;
	  splicetime+=usecond();
	} else {
	  nosplicetime-=usecond();
	  auto tmp  = Gather(source,dimension,shift,0x3,compress,face_idx);
	  is_same_node = is_same_node && tmp;
	  nosplicetime+=usecond();
	}
      } else {
	if(splice_dim){
	  splicetime-=usecond();
	  // if checkerboard is unfavourable take two passes
	  // both with block stride loop iteration
	  auto tmp1 =  GatherSimd(source,dimension,shift,0x1,compress,face_idx);
	  auto tmp2 =  GatherSimd(source,dimension,shift,0x2,compress,face_idx);
	  is_same_node = is_same_node && tmp1 && tmp2;
	  splicetime+=usecond();
	} else {
	  nosplicetime-=usecond();
	  auto tmp1 = Gather(source,dimension,shift,0x1,compress,face_idx);
	  auto tmp2 = Gather(source,dimension,shift,0x2,compress,face_idx);
	  is_same_node = is_same_node && tmp1 && tmp2;
	  nosplicetime+=usecond();
	}
      }
    }
    return is_same_node;
  }

  template<class compressor>
  void HaloGather(const Lattice<vobj> &source,compressor &compress)
  {
    mpi3synctime_g-=usecond();
    _grid->StencilBarrier();// Synch shared memory on a single nodes
    mpi3synctime_g+=usecond();

    // conformable(source.Grid(),_grid);
    assert(source.Grid()==_grid);
    halogtime-=usecond();

    u_comm_offset=0;

    // Gather all comms buffers
    int face_idx=0;
    for(int point = 0 ; point < this->_npoints; point++) {
      compress.Point(point);
      HaloGatherDir(source,compress,point,face_idx);
    }
    face_table_computed=1;
    assert(u_comm_offset==_unified_buffer_size);

    accelerator_barrier();
    halogtime+=usecond();
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
    calls++;
  }
  void AddPacket(void *xmit,void * rcv, Integer to,Integer from,Integer bytes){
    Packet p;
    p.send_buf = xmit;
    p.recv_buf = rcv;
    p.to_rank  = to;
    p.from_rank= from;
    p.bytes    = bytes;
    Packets.push_back(p);
  }
  void AddDecompress(cobj *k_p,cobj *m_p,Integer buffer_size,std::vector<Decompress> &dv) {
    Decompress d;
    d.kernel_p = k_p;
    d.mpi_p    = m_p;
    d.buffer_size = buffer_size;
    dv.push_back(d);
  }
  void AddMerge(cobj *merge_p,Vector<cobj *> &rpointers,Integer buffer_size,Integer type,std::vector<Merge> &mv) {
    Merge m;
    m.type     = type;
    m.mpointer = merge_p;
    m.vpointers= rpointers;
    m.buffer_size = buffer_size;
    mv.push_back(m);
  }
  template<class decompressor>  void CommsMerge(decompressor decompress)    {
    CommsMerge(decompress,Mergers,Decompressions);
  }
  template<class decompressor>  void CommsMergeSHM(decompressor decompress) {
    mpi3synctime-=usecond();
    _grid->StencilBarrier();// Synch shared memory on a single nodes
    mpi3synctime+=usecond();
    shmmergetime-=usecond();
    CommsMerge(decompress,MergersSHM,DecompressionsSHM);
    shmmergetime+=usecond();
  }

  template<class decompressor>
  void CommsMerge(decompressor decompress,std::vector<Merge> &mm,std::vector<Decompress> &dd) {

    mergetime-=usecond();
    for(int i=0;i<mm.size();i++){
      auto mp = &mm[i].mpointer[0];
      auto vp0= &mm[i].vpointers[0][0];
      auto vp1= &mm[i].vpointers[1][0];
      auto type= mm[i].type;
      accelerator_forNB(o,mm[i].buffer_size/2,1,{
	  decompress.Exchange(mp,vp0,vp1,type,o);
      });
    }
    mergetime+=usecond();

    decompresstime-=usecond();
    for(int i=0;i<dd.size();i++){
      auto kp = dd[i].kernel_p;
      auto mp = dd[i].mpi_p;
      accelerator_forNB(o,dd[i].buffer_size,1,{
	decompress.Decompress(kp,mp,o);
      });
    }
    decompresstime+=usecond();
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
	surface_list.push_back(site);
      }
    }
  }

  CartesianStencil(GridBase *grid,
		   int npoints,
		   int checkerboard,
		   const std::vector<int> &directions,
		   const std::vector<int> &distances,
		   Parameters p)
    : shm_bytes_thr(npoints),
      comm_bytes_thr(npoints),
      comm_enter_thr(npoints),
      comm_leave_thr(npoints),
      comm_time_thr(npoints)
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
    this->same_node.resize(npoints);

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

      int fd = _grid->_fdimensions[dimension];
      int rd = _grid->_rdimensions[dimension];
      this->_permute_type[point]=_grid->PermuteType(dimension);

      this->_checkerboard = checkerboard;

      //////////////////////////
      // the permute type
      //////////////////////////
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

	//	int rank           = grid->_processor;
	//	int recv_from_rank;
	//	int xmit_to_rank;

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
  int Gather(const Lattice<vobj> &rhs,int dimension,int shift,int cbmask,compressor & compress,int &face_idx)
  {
    typedef typename cobj::vector_type vector_type;
    typedef typename cobj::scalar_type scalar_type;

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

    int shm_receive_only = 1;
    for(int x=0;x<rd;x++){

      int sx        = (x+sshift)%rd;
      int comm_proc = ((x+sshift)/rd)%pd;

      if (comm_proc) {

	int words = buffer_size;
	if (cbmask != 0x3) words=words>>1;

	int bytes =  words * compress.CommDatumSize();

	int so  = sx*rhs.Grid()->_ostride[dimension]; // base offset for start of plane
	if ( !face_table_computed ) {
	  face_table.resize(face_idx+1);
	  Gather_plane_table_compute ((GridBase *)_grid,dimension,sx,cbmask,u_comm_offset,face_table[face_idx]);
	}

	//      	int rank           = _grid->_processor;
	int recv_from_rank;
	int xmit_to_rank;
	_grid->ShiftedRanks(dimension,comm_proc,xmit_to_rank,recv_from_rank);

	assert (xmit_to_rank   != _grid->ThisRank());
	assert (recv_from_rank != _grid->ThisRank());

	/////////////////////////////////////////////////////////
	// try the direct copy if possible
	/////////////////////////////////////////////////////////
	cobj *send_buf;
	cobj *recv_buf;
	if ( compress.DecompressionStep() ) {
	  recv_buf=u_simd_recv_buf[0];
	} else {
	  recv_buf=this->u_recv_buf_p;
	}

	send_buf = (cobj *)_grid->ShmBufferTranslate(xmit_to_rank,recv_buf);
	if ( send_buf==NULL ) {
	  send_buf = this->u_send_buf_p;
	}

	// Find out if we get the direct copy.
	void *success = (void *) _grid->ShmBufferTranslate(recv_from_rank,this->u_send_buf_p);
	if (success==NULL) {
	  // we found a packet that comes from MPI and contributes to this leg of stencil
	  shm_receive_only = 0;
	}

	gathertime-=usecond();
	assert(send_buf!=NULL);
	Gather_plane_simple_table(face_table[face_idx],rhs,send_buf,compress,u_comm_offset,so);  face_idx++;
	gathertime+=usecond();

	if ( compress.DecompressionStep() ) {

	  if ( shm_receive_only ) { // Early decompress before MPI is finished is possible
	    AddDecompress(&this->u_recv_buf_p[u_comm_offset],
			  &recv_buf[u_comm_offset],
			  words,DecompressionsSHM);
	  } else { // Decompress after MPI is finished
	    AddDecompress(&this->u_recv_buf_p[u_comm_offset],
			  &recv_buf[u_comm_offset],
			  words,Decompressions);
	  }

	  AddPacket((void *)&send_buf[u_comm_offset],
		    (void *)&recv_buf[u_comm_offset],
		    xmit_to_rank,
		    recv_from_rank,
		    bytes);

	} else {
	  AddPacket((void *)&send_buf[u_comm_offset],
		    (void *)&this->u_recv_buf_p[u_comm_offset],
		    xmit_to_rank,
		    recv_from_rank,
		    bytes);
	}
	u_comm_offset+=words;
      }
    }
    return shm_receive_only;
  }

  template<class compressor>
  int  GatherSimd(const Lattice<vobj> &rhs,int dimension,int shift,int cbmask,compressor &compress,int & face_idx)
  {
    const int Nsimd = _grid->Nsimd();

    const int maxl =2;// max layout in a direction
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
    //    std::cout << "SimdNew permute type "<<permute_type<<std::endl;

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
    assert(bytes*simd_layout == reduced_buffer_size*datum_bytes);

    Vector<cobj *> rpointers(maxl);
    Vector<cobj *> spointers(maxl);

    ///////////////////////////////////////////
    // Work out what to send where
    ///////////////////////////////////////////

    int cb    = (cbmask==0x2)? Odd : Even;
    int sshift= _grid->CheckerBoardShiftForCB(rhs.Checkerboard(),dimension,shift,cb);

    // loop over outer coord planes orthog to dim
    int shm_receive_only = 1;
    for(int x=0;x<rd;x++){

      int any_offnode = ( ((x+sshift)%fd) >= rd );

      if ( any_offnode ) {

	for(int i=0;i<maxl;i++){
	  spointers[i] = (cobj *) &u_simd_send_buf[i][u_comm_offset];
	}

	int sx   = (x+sshift)%rd;

	if ( !face_table_computed ) {
	  face_table.resize(face_idx+1);
	  Gather_plane_table_compute ((GridBase *)_grid,dimension,sx,cbmask,u_comm_offset,face_table[face_idx]);
	}
	gathermtime-=usecond();

	Gather_plane_exchange_table(face_table[face_idx],rhs,spointers,dimension,sx,cbmask,compress,permute_type);  face_idx++;

	gathermtime+=usecond();
	//spointers[0] -- low
	//spointers[1] -- high

	for(int i=0;i<maxl;i++){

	  int my_coor  = rd*i + x;            // self explanatory
	  int nbr_coor = my_coor+sshift;      // self explanatory

	  int nbr_proc = ((nbr_coor)/ld) % pd;// relative shift in processors
	  int nbr_lcoor= (nbr_coor%ld);       // local plane coor on neighbour node
	  int nbr_ic   = (nbr_lcoor)/rd;      // inner coord of peer simd lane "i"
	  int nbr_ox   = (nbr_lcoor%rd);      // outer coord of peer "x"

	  int nbr_plane = nbr_ic;
	  assert (sx == nbr_ox);

	  auto rp = &u_simd_recv_buf[i        ][u_comm_offset];
	  auto sp = &u_simd_send_buf[nbr_plane][u_comm_offset];

	  if(nbr_proc){

	    int recv_from_rank;
	    int xmit_to_rank;

	    _grid->ShiftedRanks(dimension,nbr_proc,xmit_to_rank,recv_from_rank);

	    // shm == receive pointer         if offnode
	    // shm == Translate[send pointer] if on node -- my view of his send pointer
	    cobj *shm = (cobj *) _grid->ShmBufferTranslate(recv_from_rank,sp);
	    if (shm==NULL) {
	      shm = rp;
	      // we found a packet that comes from MPI and contributes to this shift.
	      // is_same_node is only used in the WilsonStencil, and gets set for this point in the stencil.
	      // Kernel will add the exterior_terms except if is_same_node.
	      shm_receive_only = 0;
	      // leg of stencil
	    }
	    // if Direct, StencilSendToRecvFrom will suppress copy to a peer on node
	    // assuming above pointer flip
	    rpointers[i] = shm;

	    AddPacket((void *)sp,(void *)rp,xmit_to_rank,recv_from_rank,bytes);


	  } else {

	    rpointers[i] = sp;

	  }
	}

	if ( shm_receive_only ) {
	  AddMerge(&this->u_recv_buf_p[u_comm_offset],rpointers,reduced_buffer_size,permute_type,MergersSHM);
	} else {
	  AddMerge(&this->u_recv_buf_p[u_comm_offset],rpointers,reduced_buffer_size,permute_type,Mergers);
	}

	u_comm_offset     +=buffer_size;
      }
    }
    return shm_receive_only;
  }

  void ZeroCounters(void) {
    gathertime = 0.;
    commtime = 0.;
    mpi3synctime=0.;
    mpi3synctime_g=0.;
    shmmergetime=0.;
    for(int i=0;i<this->_npoints;i++){
      comm_time_thr[i]=0;
      comm_bytes_thr[i]=0;
      comm_enter_thr[i]=0;
      comm_leave_thr[i]=0;
      shm_bytes_thr[i]=0;
    }
    halogtime = 0.;
    mergetime = 0.;
    decompresstime = 0.;
    gathermtime = 0.;
    splicetime = 0.;
    nosplicetime = 0.;
    comms_bytes = 0.;
    shm_bytes = 0.;
    calls = 0.;
  };

  void Report(void) {
#define AVERAGE(A)
#define PRINTIT(A) AVERAGE(A); std::cout << GridLogMessage << " Stencil " << #A << " "<< A/calls<<std::endl;
    RealD NP = _grid->_Nprocessors;
    RealD NN = _grid->NodeCount();
    double t = 0;
    // if comm_time_thr is set they were all done in parallel so take the max
    // but add up the bytes
    int threaded = 0 ;
    for (int i = 0; i < 8; ++i) {
      if ( comm_time_thr[i]>0.0 ) {
	threaded = 1;
	comms_bytes += comm_bytes_thr[i];
	shm_bytes   += shm_bytes_thr[i];
	if (t < comm_time_thr[i]) t = comm_time_thr[i];
      }
    }
    if (threaded) commtime += t;

    _grid->GlobalSum(commtime);    commtime/=NP;
    if ( calls > 0. ) {
      std::cout << GridLogMessage << " Stencil calls "<<calls<<std::endl;
      PRINTIT(halogtime);
      PRINTIT(gathertime);
      PRINTIT(gathermtime);
      PRINTIT(mergetime);
      PRINTIT(decompresstime);
      if(comms_bytes>1.0){
	PRINTIT(comms_bytes);
	PRINTIT(commtime);
	std::cout << GridLogMessage << " Stencil " << comms_bytes/commtime/1000. << " GB/s per rank"<<std::endl;
	std::cout << GridLogMessage << " Stencil " << comms_bytes/commtime/1000.*NP/NN << " GB/s per node"<<std::endl;
      }
      if(shm_bytes>1.0){
	PRINTIT(shm_bytes); // X bytes + R bytes
	                    // Double this to include spin projection overhead with 2:1 ratio in wilson
	auto gatheralltime = gathertime+gathermtime;
	std::cout << GridLogMessage << " Stencil SHM " << (shm_bytes)/gatheralltime/1000. << " GB/s per rank"<<std::endl;
	std::cout << GridLogMessage << " Stencil SHM " << (shm_bytes)/gatheralltime/1000.*NP/NN << " GB/s per node"<<std::endl;

	auto all_bytes = comms_bytes+shm_bytes;
	std::cout << GridLogMessage << " Stencil SHM all " << (all_bytes)/gatheralltime/1000. << " GB/s per rank"<<std::endl;
	std::cout << GridLogMessage << " Stencil SHM all " << (all_bytes)/gatheralltime/1000.*NP/NN << " GB/s per node"<<std::endl;

	auto membytes = (shm_bytes + comms_bytes/2) // read/write
	              + (shm_bytes+comms_bytes)/2 * sizeof(vobj)/sizeof(cobj);
	std::cout << GridLogMessage << " Stencil SHM mem " << (membytes)/gatheralltime/1000. << " GB/s per rank"<<std::endl;
	std::cout << GridLogMessage << " Stencil SHM mem " << (membytes)/gatheralltime/1000.*NP/NN << " GB/s per node"<<std::endl;
      }
      /*
      PRINTIT(mpi3synctime);
      PRINTIT(mpi3synctime_g);
      PRINTIT(shmmergetime);
      PRINTIT(splicetime);
      PRINTIT(nosplicetime);
      */
    }
#undef PRINTIT
#undef AVERAGE
  };

};
NAMESPACE_END(Grid);

#endif
