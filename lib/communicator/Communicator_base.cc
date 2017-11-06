    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/communicator/Communicator_none.cc

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
#include <fcntl.h>
#include <unistd.h>
#include <limits.h>
#include <sys/mman.h>

namespace Grid {

///////////////////////////////////////////////////////////////
// Info that is setup once and indept of cartesian layout
///////////////////////////////////////////////////////////////
void *              CartesianCommunicator::ShmCommBuf;
uint64_t            CartesianCommunicator::MAX_MPI_SHM_BYTES   = 1024LL*1024LL*1024LL; 
CartesianCommunicator::CommunicatorPolicy_t  
CartesianCommunicator::CommunicatorPolicy= CartesianCommunicator::CommunicatorPolicyConcurrent;
int CartesianCommunicator::nCommThreads = -1;
int CartesianCommunicator::Hugepages = 0;

/////////////////////////////////
// Alloc, free shmem region
/////////////////////////////////
void *CartesianCommunicator::ShmBufferMalloc(size_t bytes){
  //  bytes = (bytes+sizeof(vRealD))&(~(sizeof(vRealD)-1));// align up bytes
  void *ptr = (void *)heap_top;
  heap_top  += bytes;
  heap_bytes+= bytes;
  if (heap_bytes >= MAX_MPI_SHM_BYTES) {
    std::cout<< " ShmBufferMalloc exceeded shared heap size -- try increasing with --shm <MB> flag" <<std::endl;
    std::cout<< " Parameter specified in units of MB (megabytes) " <<std::endl;
    std::cout<< " Current value is " << (MAX_MPI_SHM_BYTES/(1024*1024)) <<std::endl;
    assert(heap_bytes<MAX_MPI_SHM_BYTES);
  }
  return ptr;
}
void CartesianCommunicator::ShmBufferFreeAll(void) { 
  heap_top  =(size_t)ShmBufferSelf();
  heap_bytes=0;
}

/////////////////////////////////
// Grid information queries
/////////////////////////////////
int                      CartesianCommunicator::Dimensions(void)        { return _ndimension; };
int                      CartesianCommunicator::IsBoss(void)            { return _processor==0; };
int                      CartesianCommunicator::BossRank(void)          { return 0; };
int                      CartesianCommunicator::ThisRank(void)          { return _processor; };
const std::vector<int> & CartesianCommunicator::ThisProcessorCoor(void) { return _processor_coor; };
const std::vector<int> & CartesianCommunicator::ProcessorGrid(void)     { return _processors; };
int                      CartesianCommunicator::ProcessorCount(void)    { return _Nprocessors; };

////////////////////////////////////////////////////////////////////////////////
// very VERY rarely (Log, serial RNG) we need world without a grid
////////////////////////////////////////////////////////////////////////////////

void CartesianCommunicator::GlobalSum(ComplexF &c)
{
  GlobalSumVector((float *)&c,2);
}
void CartesianCommunicator::GlobalSumVector(ComplexF *c,int N)
{
  GlobalSumVector((float *)c,2*N);
}
void CartesianCommunicator::GlobalSum(ComplexD &c)
{
  GlobalSumVector((double *)&c,2);
}
void CartesianCommunicator::GlobalSumVector(ComplexD *c,int N)
{
  GlobalSumVector((double *)c,2*N);
}


#if defined( GRID_COMMS_MPI) || defined (GRID_COMMS_MPIT) || defined (GRID_COMMS_MPI3)

CartesianCommunicator::CartesianCommunicator(const std::vector<int> &processors,const CartesianCommunicator &parent,int &srank) 
{
  _ndimension = processors.size();
  assert(_ndimension = parent._ndimension);
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  // split the communicator
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  int Nparent;
  MPI_Comm_size(parent.communicator,&Nparent);

  int childsize=1;
  for(int d=0;d<processors.size();d++) {
    childsize *= processors[d];
  }
  int Nchild = Nparent/childsize;
  assert (childsize * Nchild == Nparent);

  std::vector<int> ccoor(_ndimension); // coor within subcommunicator
  std::vector<int> scoor(_ndimension); // coor of split within parent
  std::vector<int> ssize(_ndimension); // coor of split within parent

  for(int d=0;d<_ndimension;d++){
    ccoor[d] = parent._processor_coor[d] % processors[d];
    scoor[d] = parent._processor_coor[d] / processors[d];
    ssize[d] = parent._processors[d]     / processors[d];
  }
  int crank;  // rank within subcomm ; srank is rank of subcomm within blocks of subcomms
  // Mpi uses the reverse Lexico convention to us
  Lexicographic::IndexFromCoorReversed(ccoor,crank,processors);
  Lexicographic::IndexFromCoorReversed(scoor,srank,ssize);

  MPI_Comm comm_split;
  if ( Nchild > 1 ) { 

    /*
    std::cout << GridLogMessage<<"Child communicator of "<< std::hex << parent.communicator << std::dec<<std::endl;
    std::cout << GridLogMessage<<" parent grid["<< parent._ndimension<<"]    ";
    for(int d=0;d<parent._processors.size();d++)  std::cout << parent._processors[d] << " ";
    std::cout<<std::endl;

    std::cout << GridLogMessage<<" child grid["<< _ndimension <<"]    ";
    for(int d=0;d<processors.size();d++)  std::cout << processors[d] << " ";
    std::cout<<std::endl;

    std::cout << GridLogMessage<<" old rank "<< parent._processor<<" coor ["<< _ndimension <<"]    ";
    for(int d=0;d<processors.size();d++)  std::cout << parent._processor_coor[d] << " ";
    std::cout<<std::endl;

    std::cout << GridLogMessage<<" new rank "<< crank<<" coor ["<< _ndimension <<"]    ";
    for(int d=0;d<processors.size();d++)  std::cout << ccoor[d] << " ";
    std::cout<<std::endl;

    std::cout << GridLogMessage<<" new coor ["<< _ndimension <<"]    ";
    for(int d=0;d<processors.size();d++)  std::cout << parent._processor_coor[d] << " ";
    std::cout<<std::endl;
    */

    int ierr= MPI_Comm_split(parent.communicator,srank,crank,&comm_split);
    assert(ierr==0);
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    // Declare victory
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    /*
    std::cout << GridLogMessage<<"Divided communicator "<< parent._Nprocessors<<" into "
	      << Nchild <<" communicators with " << childsize << " ranks"<<std::endl;
    */
  } else {
    comm_split=parent.communicator;
    srank = 0;
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  // Set up from the new split communicator
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  InitFromMPICommunicator(processors,comm_split);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////
// Take an MPI_Comm and self assemble
//////////////////////////////////////////////////////////////////////////////////////////////////////
void CartesianCommunicator::InitFromMPICommunicator(const std::vector<int> &processors, MPI_Comm communicator_base)
{
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

  std::vector<int> periodic(_ndimension,1);
  MPI_Cart_create(communicator_base, _ndimension,&_processors[0],&periodic[0],0,&communicator);
  MPI_Comm_rank(communicator,&_processor);
  MPI_Cart_coords(communicator,_processor,_ndimension,&_processor_coor[0]);

  if ( communicator_base != communicator_world ) {
    std::cout << "Cartesian communicator created with a non-world communicator"<<std::endl;
    
    std::cout << " new communicator rank "<<_processor<< " coor ["<<_ndimension<<"] ";
    for(int d=0;d<_processors.size();d++){
      std::cout << _processor_coor[d]<<" ";
    }
    std::cout << std::endl;
  }

  int Size;
  MPI_Comm_size(communicator,&Size);

#ifdef GRID_COMMS_MPIT
  communicator_halo.resize (2*_ndimension);
  for(int i=0;i<_ndimension*2;i++){
    MPI_Comm_dup(communicator,&communicator_halo[i]);
  }
#endif
  
  assert(Size==_Nprocessors);
}

CartesianCommunicator::CartesianCommunicator(const std::vector<int> &processors) 
{
  InitFromMPICommunicator(processors,communicator_world);
}

#endif

#if !defined( GRID_COMMS_MPI3) 

int                      CartesianCommunicator::NodeCount(void)    { return ProcessorCount();};
int                      CartesianCommunicator::RankCount(void)    { return ProcessorCount();};
#endif
#if !defined( GRID_COMMS_MPI3) && !defined (GRID_COMMS_MPIT)
double CartesianCommunicator::StencilSendToRecvFrom( void *xmit,
						     int xmit_to_rank,
						     void *recv,
						     int recv_from_rank,
						     int bytes, int dir)
{
  std::vector<CommsRequest_t> list;
  // Discard the "dir"
  SendToRecvFromBegin   (list,xmit,xmit_to_rank,recv,recv_from_rank,bytes);
  SendToRecvFromComplete(list);
  return 2.0*bytes;
}
double CartesianCommunicator::StencilSendToRecvFromBegin(std::vector<CommsRequest_t> &list,
							 void *xmit,
							 int xmit_to_rank,
							 void *recv,
							 int recv_from_rank,
							 int bytes, int dir)
{
  // Discard the "dir"
  SendToRecvFromBegin(list,xmit,xmit_to_rank,recv,recv_from_rank,bytes);
  return 2.0*bytes;
}
void CartesianCommunicator::StencilSendToRecvFromComplete(std::vector<CommsRequest_t> &waitall,int dir)
{
  SendToRecvFromComplete(waitall);
}
#endif

#if !defined( GRID_COMMS_MPI3) 

void CartesianCommunicator::StencilBarrier(void){};

commVector<uint8_t> CartesianCommunicator::ShmBufStorageVector;

void *CartesianCommunicator::ShmBufferSelf(void) { return ShmCommBuf; }

void *CartesianCommunicator::ShmBuffer(int rank) {
  return NULL;
}
void *CartesianCommunicator::ShmBufferTranslate(int rank,void * local_p) { 
  return NULL;
}
void CartesianCommunicator::ShmInitGeneric(void){
#if 1
  int mmap_flag =0;
#ifdef MAP_ANONYMOUS
  mmap_flag = mmap_flag| MAP_SHARED | MAP_ANONYMOUS;
#endif
#ifdef MAP_ANON
  mmap_flag = mmap_flag| MAP_SHARED | MAP_ANON;
#endif
#ifdef MAP_HUGETLB
  if ( Hugepages ) mmap_flag |= MAP_HUGETLB;
#endif
  ShmCommBuf =(void *) mmap(NULL, MAX_MPI_SHM_BYTES, PROT_READ | PROT_WRITE, mmap_flag, -1, 0); 
  if (ShmCommBuf == (void *)MAP_FAILED) {
    perror("mmap failed ");
    exit(EXIT_FAILURE);  
  }
#ifdef MADV_HUGEPAGE
  if (!Hugepages ) madvise(ShmCommBuf,MAX_MPI_SHM_BYTES,MADV_HUGEPAGE);
#endif
#else 
  ShmBufStorageVector.resize(MAX_MPI_SHM_BYTES);
  ShmCommBuf=(void *)&ShmBufStorageVector[0];
#endif
  bzero(ShmCommBuf,MAX_MPI_SHM_BYTES);
}

#endif
  
}

