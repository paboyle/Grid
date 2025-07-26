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

#ifdef GRID_CHECKSUM_COMMS
extern void * Grid_backtrace_buffer[_NBACKTRACE];
#endif

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

#ifdef GRID_CHECKSUM_COMMS
  checksumIndex = 0;
#endif
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
#ifdef USE_GRID_REDUCTION
void CartesianCommunicator::GlobalSum(float &f){
  FlightRecorder::StepLog("GlobalSumP2P");
  CartesianCommunicator::GlobalSumP2P(f);
}
void CartesianCommunicator::GlobalSum(double &d)
{
  FlightRecorder::StepLog("GlobalSumP2P");
  CartesianCommunicator::GlobalSumP2P(d);
}
#else
void CartesianCommunicator::GlobalSum(float &f){
  FlightRecorder::StepLog("AllReduce");
  int ierr=MPI_Allreduce(MPI_IN_PLACE,&f,1,MPI_FLOAT,MPI_SUM,communicator);
  assert(ierr==0);
}
void CartesianCommunicator::GlobalSum(double &d)
{
  FlightRecorder::StepLog("AllReduce");
  int ierr = MPI_Allreduce(MPI_IN_PLACE,&d,1,MPI_DOUBLE,MPI_SUM,communicator);
  assert(ierr==0);
}
#endif
void CartesianCommunicator::GlobalSum(uint32_t &u){
  FlightRecorder::StepLog("AllReduce");
  int ierr=MPI_Allreduce(MPI_IN_PLACE,&u,1,MPI_UINT32_T,MPI_SUM,communicator);
  assert(ierr==0);
}
void CartesianCommunicator::GlobalSum(uint64_t &u){
  FlightRecorder::StepLog("AllReduce");
  int ierr=MPI_Allreduce(MPI_IN_PLACE,&u,1,MPI_UINT64_T,MPI_SUM,communicator);
  assert(ierr==0);
}
void CartesianCommunicator::GlobalSumVector(uint64_t* u,int N){
  FlightRecorder::StepLog("AllReduceVector");
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
void CartesianCommunicator::GlobalSumVector(float *f,int N)
{
  int ierr=MPI_Allreduce(MPI_IN_PLACE,f,N,MPI_FLOAT,MPI_SUM,communicator);
  assert(ierr==0);
}
void CartesianCommunicator::GlobalSumVector(double *d,int N)
{
  int ierr = MPI_Allreduce(MPI_IN_PLACE,d,N,MPI_DOUBLE,MPI_SUM,communicator);
  assert(ierr==0);
}

void CartesianCommunicator::SendToRecvFromBegin(std::vector<MpiCommsRequest_t> &list,
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
void CartesianCommunicator::CommsComplete(std::vector<MpiCommsRequest_t> &list)
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
  std::vector<MpiCommsRequest_t> reqs(0);

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

}
// Basic Halo comms primitive
double CartesianCommunicator::StencilSendToRecvFrom( void *xmit,
						     int dest, int dox,
						     void *recv,
						     int from, int dor,
						     int bytes,int dir)
{
  std::vector<CommsRequest_t> list;
  double offbytes = StencilSendToRecvFromPrepare(list,xmit,dest,dox,recv,from,dor,bytes,bytes,dir);
  offbytes       += StencilSendToRecvFromBegin(list,xmit,dest,dox,recv,from,dor,bytes,bytes,dir);
  StencilSendToRecvFromComplete(list,dir);
  return offbytes;
}


#ifdef ACCELERATOR_AWARE_MPI
void CartesianCommunicator::StencilSendToRecvFromPollIRecv(std::vector<CommsRequest_t> &list) {};
void CartesianCommunicator::StencilSendToRecvFromPollDtoH(std::vector<CommsRequest_t> &list) {};
double CartesianCommunicator::StencilSendToRecvFromPrepare(std::vector<CommsRequest_t> &list,
							   void *xmit,
							   int dest,int dox,
							   void *recv,
							   int from,int dor,
							   int xbytes,int rbytes,int dir)
{
  return 0.0; // Do nothing -- no preparation required
}
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
    else { 
      void *shm = (void *) this->ShmBufferTranslate(from,xmit);
      assert(shm!=NULL);
      acceleratorCopyDeviceToDeviceAsynch(shm,recv,rbytes);
    }
#endif
  }
  // This is a NVLINK PUT  
  if (dox) {
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
  /*finishes Get/Put*/
  acceleratorCopySynchronise();

  if (nreq==0) return;
  std::vector<MPI_Status> status(nreq);
  int ierr = MPI_Waitall(nreq,&list[0],&status[0]);
  assert(ierr==0);
  list.resize(0);
  this->StencilBarrier(); 
}

#else /* NOT     ... ACCELERATOR_AWARE_MPI */
///////////////////////////////////////////
// Pipeline mode through host memory
///////////////////////////////////////////
  /*
   * In prepare (phase 1):
   * PHASE 1: (prepare)
   * - post MPI receive buffers asynch
   * - post device - host send buffer transfer asynch
   * PHASE 2: (Begin)
   * - complete all copies
   * - post MPI send asynch
   * - post device - device transfers
   * PHASE 3: (Complete)
   * - MPI_waitall
   * - host-device transfers
   *
   *********************************
   * NB could split this further:
   *--------------------------------
   * PHASE 1: (Prepare)
   * - post MPI receive buffers asynch
   * - post device - host send buffer transfer asynch
   * PHASE 2: (BeginInterNode)
   * - complete all copies 
   * - post MPI send asynch
   * PHASE 3: (BeginIntraNode)
   * - post device - device transfers
   * PHASE 4: (Complete)
   * - MPI_waitall
   * - host-device transfers asynch
   * - (complete all copies) 
   */
double CartesianCommunicator::StencilSendToRecvFromPrepare(std::vector<CommsRequest_t> &list,
							   void *xmit,
							   int dest,int dox,
							   void *recv,
							   int from,int dor,
							   int xbytes,int rbytes,int dir)
{
/*
 * Bring sequence from Stencil.h down to lower level.
 * Assume using XeLink is ok
 */  
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

  void * host_recv = NULL;
  void * host_xmit = NULL;

  /*
   * PHASE 1: (Prepare)
   * - post MPI receive buffers asynch
   * - post device - host send buffer transfer asynch
   */

#ifdef GRID_CHECKSUM_COMMS
  rbytes += 8;
  xbytes += 8;
#endif

  if ( dor ) {
    if ( (gfrom ==MPI_UNDEFINED) || Stencil_force_mpi ) {
      tag= dir+from*32;

      host_recv = this->HostBufferMalloc(rbytes);
      ierr=MPI_Irecv(host_recv, rbytes, MPI_CHAR,from,tag,communicator_halo[commdir],&rrq);

      assert(ierr==0);
      CommsRequest_t srq;
      srq.PacketType = InterNodeRecv;
      srq.bytes      = rbytes;
      srq.req        = rrq;
      srq.host_buf   = host_recv;
      srq.device_buf = recv;
      srq.tag        = tag;
      list.push_back(srq);
      off_node_bytes+=rbytes;
    }
  }
  
  if (dox) {
    if ( (gdest == MPI_UNDEFINED) || Stencil_force_mpi ) {

      tag= dir+_processor*32;

      host_xmit = this->HostBufferMalloc(xbytes);

      CommsRequest_t srq;

#ifdef GRID_CHECKSUM_COMMS
      uint64_t xbytes_data = xbytes - 8;
      srq.ev = acceleratorCopyFromDeviceAsynch(xmit, host_xmit,xbytes_data); // Make this Asynch
      assert(xbytes % 8 == 0);
      *(uint64_t*)(((char*)host_xmit) + xbytes_data) = checksum_gpu((uint64_t*)xmit, xbytes_data / 8) ^ (checksumIndex + 1 + 1000 * tag); // flip one bit so that a zero buffer is not consistent
#else
      srq.ev = acceleratorCopyFromDeviceAsynch(xmit, host_xmit,xbytes); // Make this Asynch
#endif
      
      //      ierr =MPI_Isend(host_xmit, xbytes, MPI_CHAR,dest,tag,communicator_halo[commdir],&xrq);
      //      assert(ierr==0);
      //      off_node_bytes+=xbytes;

      srq.PacketType = InterNodeXmit;
      srq.bytes      = xbytes;
      //      srq.req        = xrq;
      srq.host_buf   = host_xmit;
      srq.device_buf = xmit;
      srq.tag        = tag;
      srq.dest       = dest;
      srq.commdir    = commdir;
      list.push_back(srq);
    }
  }

  return off_node_bytes;
}
/*
 * In the interest of better pipelining, poll for completion on each DtoH and 
 * start MPI_ISend in the meantime
 */
void CartesianCommunicator::StencilSendToRecvFromPollIRecv(std::vector<CommsRequest_t> &list)
{
  int pending = 0;
  do {

    pending = 0;

    for(int idx = 0; idx<list.size();idx++){

      if ( list[idx].PacketType==InterNodeRecv ) {

	int flag = 0;
	MPI_Status status;
	int ierr = MPI_Test(&list[idx].req,&flag,&status);
	assert(ierr==0);

	if ( flag ) {
	  //	  std::cout << " PollIrecv "<<idx<<" flag "<<flag<<std::endl;
#ifdef GRID_CHECKSUM_COMMS
	  acceleratorCopyToDeviceAsynch(list[idx].host_buf,list[idx].device_buf,list[idx].bytes - 8);
#else
	  acceleratorCopyToDeviceAsynch(list[idx].host_buf,list[idx].device_buf,list[idx].bytes);
#endif
	  list[idx].PacketType=InterNodeReceiveHtoD;
	} else {
	  pending ++;
	}
      }
    }
    //    std::cout << " PollIrecv "<<pending<<" pending requests"<<std::endl;
  } while ( pending );
  
}
void CartesianCommunicator::StencilSendToRecvFromPollDtoH(std::vector<CommsRequest_t> &list)
{
  int pending = 0;
  do {

    pending = 0;

    for(int idx = 0; idx<list.size();idx++){

      if ( list[idx].PacketType==InterNodeXmit ) {

	if ( acceleratorEventIsComplete(list[idx].ev) ) {

	  void *host_xmit = list[idx].host_buf;
	  uint32_t xbytes = list[idx].bytes;
	  int dest        = list[idx].dest;
	  int tag         = list[idx].tag;
	  int commdir     = list[idx].commdir;
	  ///////////////////
	  // Send packet
	  ///////////////////

	  //	  std::cout << " DtoH is complete for index "<<idx<<" calling MPI_Isend "<<std::endl;
	  
	  MPI_Request xrq;
	  int ierr =MPI_Isend(host_xmit, xbytes, MPI_CHAR,dest,tag,communicator_halo[commdir],&xrq);
	  assert(ierr==0);

	  list[idx].req        = xrq; // Update the MPI request in the list

	  list[idx].PacketType=InterNodeXmitISend;

	} else {
	  // not done, so return to polling loop
	  pending++;
	}
      }
    }
  } while (pending);
}  

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

  void * host_xmit = NULL;

  ////////////////////////////////
  // Receives already posted
  // Copies already started
  ////////////////////////////////
  /*  
   * PHASE 2: (Begin)
   * - complete all copies
   * - post MPI send asynch
   */
#ifdef NVLINK_GET
  if ( dor ) {

    if ( ! ( (gfrom ==MPI_UNDEFINED) || Stencil_force_mpi ) ) {
      // Intranode
      void *shm = (void *) this->ShmBufferTranslate(from,xmit);
      assert(shm!=NULL);

      CommsRequest_t srq;

      srq.ev = acceleratorCopyDeviceToDeviceAsynch(shm,recv,rbytes);

      srq.PacketType = IntraNodeRecv;
      srq.bytes      = xbytes;
      //      srq.req        = xrq;
      srq.host_buf   = NULL;
      srq.device_buf = xmit;
      srq.tag        = -1;
      srq.dest       = dest;
      srq.commdir    = dir;
      list.push_back(srq);
    }
  }  
#else
  if (dox) {

    if ( !( (gdest == MPI_UNDEFINED) || Stencil_force_mpi ) ) {
      // Intranode
      void *shm = (void *) this->ShmBufferTranslate(dest,recv);
      assert(shm!=NULL);

      CommsRequest_t srq;
      
      srq.ev = acceleratorCopyDeviceToDeviceAsynch(xmit,shm,xbytes);

      srq.PacketType = IntraNodeXmit;
      srq.bytes      = xbytes;
      //      srq.req        = xrq;
      srq.host_buf   = NULL;
      srq.device_buf = xmit;
      srq.tag        = -1;
      srq.dest       = dest;
      srq.commdir    = dir;
      list.push_back(srq);
      
    }
  }
#endif
  return off_node_bytes;
}
void CartesianCommunicator::StencilSendToRecvFromComplete(std::vector<CommsRequest_t> &list,int dir)
{
  acceleratorCopySynchronise(); // Complete all pending copy transfers D2D

  std::vector<MPI_Status> status;
  std::vector<MPI_Request> MpiRequests;
    
  for(int r=0;r<list.size();r++){
    // Must check each Send buf is clear to reuse
    if ( list[r].PacketType == InterNodeXmitISend ) MpiRequests.push_back(list[r].req);
    //    if ( list[r].PacketType == InterNodeRecv ) MpiRequests.push_back(list[r].req); // Already "Test" passed
  }

  int nreq=MpiRequests.size();

  if (nreq>0) {
    status.resize(MpiRequests.size());
    int ierr = MPI_Waitall(MpiRequests.size(),&MpiRequests[0],&status[0]); // Sends are guaranteed in order. No harm in not completing.
    assert(ierr==0);
  }
  
  //  for(int r=0;r<nreq;r++){
  //    if ( list[r].PacketType==InterNodeRecv ) {
  //      acceleratorCopyToDeviceAsynch(list[r].host_buf,list[r].device_buf,list[r].bytes);
  //    }
  //  }

#ifdef GRID_CHECKSUM_COMMS
  for(int r=0;r<list.size();r++){
    if ( list[r].PacketType == InterNodeReceiveHtoD ) {
      uint64_t rbytes_data = list[r].bytes - 8;
      uint64_t expected_cs = *(uint64_t*)(((char*)list[r].host_buf) + rbytes_data);
      uint64_t computed_cs = checksum_gpu((uint64_t*)list[r].device_buf, rbytes_data / 8) ^ (checksumIndex + 1 + 1000 * list[r].tag); //
      if (expected_cs != computed_cs) {
	// TODO: error message, backtrace, quit

	fprintf(stderr, "GRID_CHECKSUM_COMMS error:\n");
	fprintf(stderr, " processor = %d\n", (int)_processor);
	for(int d=0;d<_processors.size();d++)
	  fprintf(stderr, " processor_coord[%d] = %d\n", d, _processor_coor[d]);
	fprintf(stderr, " hostname: %s\n", GridHostname());
	fprintf(stderr, " expected_cs: %ld\n", expected_cs);
	fprintf(stderr, " computed_cs: %ld\n", computed_cs);
	fprintf(stderr, " dest: %d\n", list[r].dest);
	fprintf(stderr, " tag: %d\n", list[r].tag);
	fprintf(stderr, " commdir: %d\n", list[r].commdir);
	fprintf(stderr, " bytes: %ld\n", (uint64_t)list[r].bytes);

	fflush(stderr);

	// backtrace
	int symbols = backtrace(Grid_backtrace_buffer,_NBACKTRACE);
	backtrace_symbols_fd(Grid_backtrace_buffer,symbols, 2);

	exit(1);
      }
    }
  }

  incrementChecksumIndex();
#endif

  list.resize(0);               // Delete the list
  this->HostBufferFreeAll();    // Clean up the buffer allocs
#ifndef NVLINK_GET
  this->StencilBarrier(); // if PUT must check our nbrs have filled our receive buffers.
#endif   
}
#endif
////////////////////////////////////////////
// END PIPELINE MODE / NO CUDA AWARE MPI
////////////////////////////////////////////

void CartesianCommunicator::StencilBarrier(void)
{
  FlightRecorder::StepLog("NodeBarrier");
  MPI_Barrier  (ShmComm);
}
//void CartesianCommunicator::SendToRecvFromComplete(std::vector<CommsRequest_t> &list)
//{
//}
void CartesianCommunicator::Barrier(void)
{
  FlightRecorder::StepLog("GridBarrier");
  int ierr = MPI_Barrier(communicator);
  assert(ierr==0);
}
void CartesianCommunicator::Broadcast(int root,void* data, int bytes)
{
  FlightRecorder::StepLog("Broadcast");
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
  FlightRecorder::StepLog("BroadcastWorld");
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
  FlightRecorder::StepLog("AllToAll");
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
