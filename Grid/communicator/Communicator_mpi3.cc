/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./lib/communicator/Communicator_mpi.cc

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
#include <Grid/GridCore.h>
#include <Grid/communicator/SharedMemory.h>

NAMESPACE_BEGIN(Grid);

Grid_MPI_Comm       CartesianCommunicator::communicator_world;

////////////////////////////////////////////
// First initialise of comms system
////////////////////////////////////////////
void CartesianCommunicator::Init(int *argc, char ***argv)
{

  int flag;
  int provided;

  MPI_Initialized(&flag); // needed to coexist with other libs apparently
  if ( !flag ) {

#ifndef GRID_COMMS_THREADS
    nCommThreads=1;
    // wrong results here too
    // For now: comms-overlap leads to wrong results in Benchmark_wilson even on single node MPI runs
    // other comms schemes are ok
    MPI_Init_thread(argc,argv,MPI_THREAD_SERIALIZED,&provided);
#else
    MPI_Init_thread(argc,argv,MPI_THREAD_MULTIPLE,&provided);
#endif
    //If only 1 comms thread we require any threading mode other than SINGLE, but for multiple comms threads we need MULTIPLE
    if( (nCommThreads == 1) && (provided == MPI_THREAD_SINGLE) ) {
      assert(0);
    }

    if( (nCommThreads > 1) && (provided != MPI_THREAD_MULTIPLE) ) {
      assert(0);
    }
  }

  // Never clean up as done once.
  MPI_Comm_dup (MPI_COMM_WORLD,&communicator_world);

  Grid_quiesce_nodes();
  GlobalSharedMemory::Init(communicator_world);
  GlobalSharedMemory::SharedMemoryAllocate(
		   GlobalSharedMemory::MAX_MPI_SHM_BYTES,
		   GlobalSharedMemory::Hugepages);
  Grid_unquiesce_nodes();
}

///////////////////////////////////////////////////////////////////////////
// Use cartesian communicators now even in MPI3
///////////////////////////////////////////////////////////////////////////
void CartesianCommunicator::ShiftedRanks(int dim,int shift,int &source,int &dest)
{
  int ierr=MPI_Cart_shift(communicator,dim,shift,&source,&dest);
  assert(ierr==0);
}
int CartesianCommunicator::RankFromProcessorCoor(Coordinate &coor)
{
  int rank;
  int ierr=MPI_Cart_rank  (communicator, &coor[0], &rank);
  assert(ierr==0);
  return rank;
}
void  CartesianCommunicator::ProcessorCoorFromRank(int rank, Coordinate &coor)
{
  coor.resize(_ndimension);
  int ierr=MPI_Cart_coords  (communicator, rank, _ndimension,&coor[0]);
  assert(ierr==0);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
// Initialises from communicator_world
////////////////////////////////////////////////////////////////////////////////////////////////////////
CartesianCommunicator::CartesianCommunicator(const Coordinate &processors)
{
  MPI_Comm optimal_comm;
  ////////////////////////////////////////////////////
  // Remap using the shared memory optimising routine
  // The remap creates a comm which must be freed
  ////////////////////////////////////////////////////
  GlobalSharedMemory::OptimalCommunicator    (processors,optimal_comm,_shm_processors);
  InitFromMPICommunicator(processors,optimal_comm);
  SetCommunicator(optimal_comm);
  ///////////////////////////////////////////////////
  // Free the temp communicator
  ///////////////////////////////////////////////////
  MPI_Comm_free(&optimal_comm);
}

//////////////////////////////////
// Try to subdivide communicator
//////////////////////////////////
CartesianCommunicator::CartesianCommunicator(const Coordinate &processors,const CartesianCommunicator &parent,int &srank)
{
  _ndimension = processors.size();  assert(_ndimension>=1);
  int parent_ndimension = parent._ndimension; assert(_ndimension >= parent._ndimension);
  Coordinate parent_processor_coor(_ndimension,0);
  Coordinate parent_processors    (_ndimension,1);
  Coordinate shm_processors       (_ndimension,1);
  // Can make 5d grid from 4d etc...
  int pad = _ndimension-parent_ndimension;
  for(int d=0;d<parent_ndimension;d++){
    parent_processor_coor[pad+d]=parent._processor_coor[d];
    parent_processors    [pad+d]=parent._processors[d];
    shm_processors       [pad+d]=parent._shm_processors[d];
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  // split the communicator
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //  int Nparent = parent._processors ;
  int Nparent;
  MPI_Comm_size(parent.communicator,&Nparent);

  int childsize=1;
  for(int d=0;d<processors.size();d++) {
    childsize *= processors[d];
  }
  int Nchild = Nparent/childsize;
  assert (childsize * Nchild == Nparent);

  Coordinate ccoor(_ndimension); // coor within subcommunicator
  Coordinate scoor(_ndimension); // coor of split within parent
  Coordinate ssize(_ndimension); // coor of split within parent

  for(int d=0;d<_ndimension;d++){
    ccoor[d] = parent_processor_coor[d] % processors[d];
    scoor[d] = parent_processor_coor[d] / processors[d];
    ssize[d] = parent_processors[d]     / processors[d];
    if ( processors[d] < shm_processors[d] ) shm_processors[d] = processors[d]; // subnode splitting.
  }

  // rank within subcomm ; srank is rank of subcomm within blocks of subcomms
  int crank;
  // Mpi uses the reverse Lexico convention to us; so reversed routines called
  Lexicographic::IndexFromCoorReversed(ccoor,crank,processors); // processors is the split grid dimensions
  Lexicographic::IndexFromCoorReversed(scoor,srank,ssize);      // ssize is the number of split grids

  MPI_Comm comm_split;
  if ( Nchild > 1 ) {

    ////////////////////////////////////////////////////////////////
    // Split the communicator
    ////////////////////////////////////////////////////////////////
    int ierr= MPI_Comm_split(parent.communicator,srank,crank,&comm_split);
    assert(ierr==0);

  } else {
    srank = 0;
    int ierr = MPI_Comm_dup (parent.communicator,&comm_split);
    assert(ierr==0);
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  // Set up from the new split communicator
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  InitFromMPICommunicator(processors,comm_split);

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  // Take the right SHM buffers
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  SetCommunicator(comm_split);

  ///////////////////////////////////////////////
  // Free the temp communicator
  ///////////////////////////////////////////////
  MPI_Comm_free(&comm_split);

  if(0){
    std::cout << " ndim " <<_ndimension<<" " << parent._ndimension << std::endl;
    for(int d=0;d<processors.size();d++){
      std::cout << d<< " " << _processor_coor[d] <<" " <<  ccoor[d]<<std::endl;
    }
  }
  for(int d=0;d<processors.size();d++){
    assert(_processor_coor[d] == ccoor[d] );
  }
}

void CartesianCommunicator::InitFromMPICommunicator(const Coordinate &processors, MPI_Comm communicator_base)
{
  ////////////////////////////////////////////////////
  // Creates communicator, and the communicator_halo
  ////////////////////////////////////////////////////
  _ndimension = processors.size();
  _processor_coor.resize(_ndimension);

  /////////////////////////////////
  // Count the requested nodes
  /////////////////////////////////
  _Nprocessors=1;
  _processors = processors;
  for(int i=0;i<_ndimension;i++){
    _Nprocessors*=_processors[i];
  }

  Coordinate periodic(_ndimension,1);
  MPI_Cart_create(communicator_base, _ndimension,&_processors[0],&periodic[0],0,&communicator);
  MPI_Comm_rank(communicator,&_processor);
  MPI_Cart_coords(communicator,_processor,_ndimension,&_processor_coor[0]);

  if ( 0 && (communicator_base != communicator_world) ) {
    std::cout << "InitFromMPICommunicator Cartesian communicator created with a non-world communicator"<<std::endl;
    std::cout << " new communicator rank "<<_processor<< " coor ["<<_ndimension<<"] ";
    for(int d=0;d<_processors.size();d++){
      std::cout << _processor_coor[d]<<" ";
    }
    std::cout << std::endl;
  }

  int Size;
  MPI_Comm_size(communicator,&Size);

  communicator_halo.resize (2*_ndimension);
  for(int i=0;i<_ndimension*2;i++){
    MPI_Comm_dup(communicator,&communicator_halo[i]);
  }
  assert(Size==_Nprocessors);
}

CartesianCommunicator::~CartesianCommunicator()
{
  int MPI_is_finalised;
  MPI_Finalized(&MPI_is_finalised);
  if (communicator && !MPI_is_finalised) {
    MPI_Comm_free(&communicator);
    for(int i=0;i<communicator_halo.size();i++){
      MPI_Comm_free(&communicator_halo[i]);
    }
  }
}
void CartesianCommunicator::GlobalSum(uint32_t &u){
  int ierr=MPI_Allreduce(MPI_IN_PLACE,&u,1,MPI_UINT32_T,MPI_SUM,communicator);
  assert(ierr==0);
}
void CartesianCommunicator::GlobalSum(uint64_t &u){
  int ierr=MPI_Allreduce(MPI_IN_PLACE,&u,1,MPI_UINT64_T,MPI_SUM,communicator);
  assert(ierr==0);
}
void CartesianCommunicator::GlobalSumVector(uint64_t* u,int N){
  int ierr=MPI_Allreduce(MPI_IN_PLACE,u,N,MPI_UINT64_T,MPI_SUM,communicator);
  assert(ierr==0);
}
void CartesianCommunicator::GlobalXOR(uint32_t &u){
  int ierr=MPI_Allreduce(MPI_IN_PLACE,&u,1,MPI_UINT32_T,MPI_BXOR,communicator);
  assert(ierr==0);
}
void CartesianCommunicator::GlobalXOR(uint64_t &u){
  int ierr=MPI_Allreduce(MPI_IN_PLACE,&u,1,MPI_UINT64_T,MPI_BXOR,communicator);
  assert(ierr==0);
}
void CartesianCommunicator::GlobalMax(float &f)
{
  int ierr=MPI_Allreduce(MPI_IN_PLACE,&f,1,MPI_FLOAT,MPI_MAX,communicator);
  assert(ierr==0);
}
void CartesianCommunicator::GlobalMax(double &d)
{
  int ierr = MPI_Allreduce(MPI_IN_PLACE,&d,1,MPI_DOUBLE,MPI_MAX,communicator);
  assert(ierr==0);
}
void CartesianCommunicator::GlobalSum(float &f){
  int ierr=MPI_Allreduce(MPI_IN_PLACE,&f,1,MPI_FLOAT,MPI_SUM,communicator);
  assert(ierr==0);
}
void CartesianCommunicator::GlobalSumVector(float *f,int N)
{
  int ierr=MPI_Allreduce(MPI_IN_PLACE,f,N,MPI_FLOAT,MPI_SUM,communicator);
  assert(ierr==0);
}
void CartesianCommunicator::GlobalSum(double &d)
{
  int ierr = MPI_Allreduce(MPI_IN_PLACE,&d,1,MPI_DOUBLE,MPI_SUM,communicator);
  assert(ierr==0);
}
void CartesianCommunicator::GlobalSumVector(double *d,int N)
{
  int ierr = MPI_Allreduce(MPI_IN_PLACE,d,N,MPI_DOUBLE,MPI_SUM,communicator);
  assert(ierr==0);
}

void CartesianCommunicator::SendToRecvFromBegin(std::vector<CommsRequest_t> &list,
						void *xmit,
						int dest,
						void *recv,
						int from,
						int bytes,int dir)
{
  MPI_Request xrq;
  MPI_Request rrq;

  assert(dest != _processor);
  assert(from != _processor);

  int tag;

  tag= dir+from*32;
  int ierr=MPI_Irecv(recv, bytes, MPI_CHAR,from,tag,communicator,&rrq);
  assert(ierr==0);
  list.push_back(rrq);
  
  tag= dir+_processor*32;
  ierr =MPI_Isend(xmit, bytes, MPI_CHAR,dest,tag,communicator,&xrq);
  assert(ierr==0);
  list.push_back(xrq);
}
void CartesianCommunicator::CommsComplete(std::vector<CommsRequest_t> &list)
{
  int nreq=list.size();

  if (nreq==0) return;

  std::vector<MPI_Status> status(nreq);
  int ierr = MPI_Waitall(nreq,&list[0],&status[0]);
  assert(ierr==0);
  list.resize(0);
}

// Basic Halo comms primitive
void CartesianCommunicator::SendToRecvFrom(void *xmit,
					   int dest,
					   void *recv,
					   int from,
					   int bytes)
{
  std::vector<CommsRequest_t> reqs(0);
  unsigned long  xcrc = crc32(0L, Z_NULL, 0);
  unsigned long  rcrc = crc32(0L, Z_NULL, 0);

  int myrank = _processor;
  int ierr;

  // Enforce no UVM in comms, device or host OK
  assert(acceleratorIsCommunicable(xmit));
  assert(acceleratorIsCommunicable(recv));

  // Give the CPU to MPI immediately; can use threads to overlap optionally
  //  printf("proc %d SendToRecvFrom %d bytes Sendrecv \n",_processor,bytes);
  ierr=MPI_Sendrecv(xmit,bytes,MPI_CHAR,dest,myrank,
		    recv,bytes,MPI_CHAR,from, from,
		    communicator,MPI_STATUS_IGNORE);
  assert(ierr==0);

  //  xcrc = crc32(xcrc,(unsigned char *)xmit,bytes);
  //  rcrc = crc32(rcrc,(unsigned char *)recv,bytes);
  //  printf("proc %d SendToRecvFrom %d bytes xcrc %lx rcrc %lx\n",_processor,bytes,xcrc,rcrc); fflush
}
// Basic Halo comms primitive
double CartesianCommunicator::StencilSendToRecvFrom( void *xmit,
						     int dest, int dox,
						     void *recv,
						     int from, int dor,
						     int bytes,int dir)
{
  std::vector<CommsRequest_t> list;
  double offbytes = StencilSendToRecvFromBegin(list,xmit,dest,dox,recv,from,dor,bytes,bytes,dir);
  StencilSendToRecvFromComplete(list,dir);
  return offbytes;
}

#undef NVLINK_GET // Define to use get instead of put DMA
double CartesianCommunicator::StencilSendToRecvFromBegin(std::vector<CommsRequest_t> &list,
							 void *xmit,
							 int dest,int dox,
							 void *recv,
							 int from,int dor,
							 int xbytes,int rbytes,int dir)
{
  int ncomm  =communicator_halo.size();
  int commdir=dir%ncomm;

  MPI_Request xrq;
  MPI_Request rrq;

  int ierr;
  int gdest = ShmRanks[dest];
  int gfrom = ShmRanks[from];
  int gme   = ShmRanks[_processor];

  assert(dest != _processor);
  assert(from != _processor);
  assert(gme  == ShmRank);
  double off_node_bytes=0.0;
  int tag;

  if ( dor ) {
    if ( (gfrom ==MPI_UNDEFINED) || Stencil_force_mpi ) {
      tag= dir+from*32;
      ierr=MPI_Irecv(recv, rbytes, MPI_CHAR,from,tag,communicator_halo[commdir],&rrq);
      assert(ierr==0);
      list.push_back(rrq);
      off_node_bytes+=rbytes;
    }
#ifdef NVLINK_GET
      void *shm = (void *) this->ShmBufferTranslate(from,xmit);
      assert(shm!=NULL);
      acceleratorCopyDeviceToDeviceAsynch(shm,recv,rbytes);
#endif
  }
  
  if (dox) {
    //  rcrc = crc32(rcrc,(unsigned char *)recv,bytes);
    if ( (gdest == MPI_UNDEFINED) || Stencil_force_mpi ) {
      tag= dir+_processor*32;
      ierr =MPI_Isend(xmit, xbytes, MPI_CHAR,dest,tag,communicator_halo[commdir],&xrq);
      assert(ierr==0);
      list.push_back(xrq);
      off_node_bytes+=xbytes;
    } else {
#ifndef NVLINK_GET
      void *shm = (void *) this->ShmBufferTranslate(dest,recv);
      assert(shm!=NULL);
      acceleratorCopyDeviceToDeviceAsynch(xmit,shm,xbytes);
#endif
      
    }
  }

  return off_node_bytes;
}
void CartesianCommunicator::StencilSendToRecvFromComplete(std::vector<CommsRequest_t> &list,int dir)
{
  int nreq=list.size();

  acceleratorCopySynchronise();

  if (nreq==0) return;

  std::vector<MPI_Status> status(nreq);
  int ierr = MPI_Waitall(nreq,&list[0],&status[0]);
  assert(ierr==0);
  list.resize(0);
}
void CartesianCommunicator::StencilBarrier(void)
{
  MPI_Barrier  (ShmComm);
}
//void CartesianCommunicator::SendToRecvFromComplete(std::vector<CommsRequest_t> &list)
//{
//}
void CartesianCommunicator::Barrier(void)
{
  int ierr = MPI_Barrier(communicator);
  assert(ierr==0);
}
void CartesianCommunicator::Broadcast(int root,void* data, int bytes)
{
  int ierr=MPI_Bcast(data,
		     bytes,
		     MPI_BYTE,
		     root,
		     communicator);
  assert(ierr==0);
}
int CartesianCommunicator::RankWorld(void){
  int r;
  MPI_Comm_rank(communicator_world,&r);
  return r;
}
void CartesianCommunicator::BarrierWorld(void){
  int ierr = MPI_Barrier(communicator_world);
  assert(ierr==0);
}
void CartesianCommunicator::BroadcastWorld(int root,void* data, int bytes)
{
  int ierr= MPI_Bcast(data,
		      bytes,
		      MPI_BYTE,
		      root,
		      communicator_world);
  assert(ierr==0);
}

void CartesianCommunicator::AllToAll(int dim,void  *in,void *out,uint64_t words,uint64_t bytes)
{
  Coordinate row(_ndimension,1);
  assert(dim>=0 && dim<_ndimension);

  //  Split the communicator
  row[dim] = _processors[dim];

  int me;
  CartesianCommunicator Comm(row,*this,me);
  Comm.AllToAll(in,out,words,bytes);
}
void CartesianCommunicator::AllToAll(void  *in,void *out,uint64_t words,uint64_t bytes)
{
  // MPI is a pain and uses "int" arguments
  // 64*64*64*128*16 == 500Million elements of data.
  // When 24*4 bytes multiples get 50x 10^9 >>> 2x10^9 Y2K bug.
  // (Turns up on 32^3 x 64 Gparity too)
  MPI_Datatype object;
  int iwords;
  int ibytes;
  iwords = words;
  ibytes = bytes;
  assert(words == iwords); // safe to cast to int ?
  assert(bytes == ibytes); // safe to cast to int ?
  MPI_Type_contiguous(ibytes,MPI_BYTE,&object);
  MPI_Type_commit(&object);
  MPI_Alltoall(in,iwords,object,out,iwords,object,communicator);
  MPI_Type_free(&object);
}

NAMESPACE_END(Grid);
