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
#include "Grid.h"
#include <mpi.h>

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Workarounds:
/// i) bloody mac os doesn't implement unnamed semaphores since it is "optional" posix.
///    darwin dispatch semaphores don't seem to be multiprocess.
///
/// ii) openmpi under --mca shmem posix works with two squadrons per node; 
///     openmpi under default mca settings (I think --mca shmem mmap) on MacOS makes two squadrons map the SAME
///     memory as each other, despite their living on different communicators. This appears to be a bug in OpenMPI.
///
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include <semaphore.h>
typedef sem_t *Grid_semaphore;

#define SEM_INIT(S)      S = sem_open(sem_name,0,0600,0); assert ( S != SEM_FAILED );
#define SEM_INIT_EXCL(S) sem_unlink(sem_name); S = sem_open(sem_name,O_CREAT|O_EXCL,0600,0); assert ( S != SEM_FAILED );
#define SEM_POST(S) assert ( sem_post(S) == 0 ); 
#define SEM_WAIT(S) assert ( sem_wait(S) == 0 );

#include <sys/mman.h>


namespace Grid {

enum { COMMAND_ISEND, COMMAND_IRECV, COMMAND_WAITALL };

struct Descriptor {
  uint64_t buf;
  size_t bytes;
  int rank;
  int tag;
  int command;
  MPI_Request request;
};

const int pool = 48;

class SlaveState {
public:
  volatile int head;
  volatile int start;
  volatile int tail;
  volatile Descriptor Descrs[pool];
};

class Slave {
public:
  Grid_semaphore  sem_head;
  Grid_semaphore  sem_tail;
  SlaveState *state;
  MPI_Comm squadron;
  uint64_t     base;
  int universe_rank;
  int vertical_rank;
  char sem_name [NAME_MAX];
  ////////////////////////////////////////////////////////////
  // Descriptor circular pointers
  ////////////////////////////////////////////////////////////
  Slave() {};

  void Init(SlaveState * _state,MPI_Comm _squadron,int _universe_rank,int _vertical_rank);

  void SemInit(void) {
    sprintf(sem_name,"/Grid_mpi3_sem_head_%d",universe_rank);
    printf("SEM_NAME: %s \n",sem_name);
    SEM_INIT(sem_head);
    sprintf(sem_name,"/Grid_mpi3_sem_tail_%d",universe_rank);
    printf("SEM_NAME: %s \n",sem_name);
    SEM_INIT(sem_tail);
  }  
  void SemInitExcl(void) {
    sprintf(sem_name,"/Grid_mpi3_sem_head_%d",universe_rank);
    printf("SEM_INIT_EXCL: %s \n",sem_name);
    SEM_INIT_EXCL(sem_head);
    sprintf(sem_name,"/Grid_mpi3_sem_tail_%d",universe_rank);
    printf("SEM_INIT_EXCL: %s \n",sem_name);
    SEM_INIT_EXCL(sem_tail);
  }  
  void WakeUpDMA(void) { 
    SEM_POST(sem_head);
  };
  void WakeUpCompute(void) { 
    SEM_POST(sem_tail);
  };
  void WaitForCommand(void) { 
    SEM_WAIT(sem_head);
  };
  void WaitForComplete(void) { 
    SEM_WAIT(sem_tail);
  };
  void EventLoop (void) {
    std::cout<< " Entering event loop "<<std::endl;
    while(1){
      WaitForCommand();
      //      std::cout << "Getting command "<<std::endl;
      Event();
    }
  }

  int Event (void) ;

  uint64_t QueueCommand(int command,void *buf, int bytes, int hashtag, MPI_Comm comm,int u_rank) ;

  void WaitAll() {
    //    std::cout << "Queueing WAIT command  "<<std::endl;
    QueueCommand(COMMAND_WAITALL,0,0,0,squadron,0);
    //    std::cout << "Waking up DMA "<<std::endl;
    WakeUpDMA();
    //    std::cout << "Waiting from semaphore "<<std::endl;
    WaitForComplete();
    //    std::cout << "Checking FIFO is empty "<<std::endl;
    assert ( state->tail == state->head );
  }
};

////////////////////////////////////////////////////////////////////////
// One instance of a data mover.
// Master and Slave must agree on location in shared memory
////////////////////////////////////////////////////////////////////////

class MPIoffloadEngine { 
public:

  static std::vector<Slave> Slaves;

  static int ShmSetup;
  
  static int UniverseRank;
  static int UniverseSize;
  
  static MPI_Comm communicator_universe;
  static MPI_Comm communicator_cached;

  static MPI_Comm HorizontalComm;
  static int HorizontalRank;
  static int HorizontalSize;
  
  static MPI_Comm VerticalComm;
  static MPI_Win  VerticalWindow; 
  static int VerticalSize;
  static int VerticalRank;
  
  static std::vector<void *> VerticalShmBufs;
  static std::vector<std::vector<int> > UniverseRanks;
  static std::vector<int> UserCommunicatorToWorldRanks; 
  
  static MPI_Group WorldGroup, CachedGroup;
  
  static void CommunicatorInit (MPI_Comm &communicator_world,
				MPI_Comm &ShmComm,
				void * &ShmCommBuf);

  static void MapCommRankToWorldRank(int &hashtag, int & comm_world_peer,int tag, MPI_Comm comm,int commrank);

  /////////////////////////////////////////////////////////
  // routines for master proc must handle any communicator
  /////////////////////////////////////////////////////////

  static void QueueSend(int slave,void *buf, int bytes, int tag, MPI_Comm comm,int rank) {
     //    std::cout<< " Queueing send  "<< bytes<< " slave "<< slave << " to comm "<<rank  <<std::endl;
    Slaves[slave].QueueCommand(COMMAND_ISEND,buf,bytes,tag,comm,rank);
    //    std::cout << "Queued send command to rank "<< rank<< " via "<<slave <<std::endl;
    Slaves[slave].WakeUpDMA();
    //    std::cout << "Waking up DMA "<< slave<<std::endl;
  };

  static void QueueRecv(int slave, void *buf, int bytes, int tag, MPI_Comm comm,int rank) {
    //    std::cout<< " Queueing recv "<< bytes<< " slave "<< slave << " from comm "<<rank  <<std::endl;
    Slaves[slave].QueueCommand(COMMAND_IRECV,buf,bytes,tag,comm,rank);
    //    std::cout << "Queued recv command from rank "<< rank<< " via "<<slave <<std::endl;
    Slaves[slave].WakeUpDMA();
    //    std::cout << "Waking up DMA "<< slave<<std::endl;
  };

  static void WaitAll() {
    for(int s=1;s<VerticalSize;s++) {
      //      std::cout << "Waiting for slave "<< s<<std::endl;
      Slaves[s].WaitAll();
    }
    //    std::cout << " Wait all Complete "<<std::endl;
  };

  static void GetWork(int nwork, int me, int & mywork, int & myoff,int units){
    int basework = nwork/units;
    int backfill = units-(nwork%units);
    if ( me >= units ) { 
      mywork = myoff = 0;
    } else { 
      mywork = (nwork+me)/units;
      myoff  = basework * me;
      if ( me > backfill ) 
	myoff+= (me-backfill);
    }
    return;
  };

  static void QueueMultiplexedSend(void *buf, int bytes, int tag, MPI_Comm comm,int rank) {
    uint8_t * cbuf = (uint8_t *) buf;
    int mywork, myoff, procs;
    procs = VerticalSize-1;
    for(int s=0;s<procs;s++) {
      GetWork(bytes,s,mywork,myoff,procs);
      QueueSend(s+1,&cbuf[myoff],mywork,tag,comm,rank);
    }
  };

  static void QueueMultiplexedRecv(void *buf, int bytes, int tag, MPI_Comm comm,int rank) {
    uint8_t * cbuf = (uint8_t *) buf;
    int mywork, myoff, procs;
    procs = VerticalSize-1;
    for(int s=0;s<procs;s++) {
      GetWork(bytes,s,mywork,myoff,procs);
      QueueRecv(s+1,&cbuf[myoff],mywork,tag,comm,rank);
    }
  };

};


///////////////////////////////////////////////////////////////////////////////////////////////////
// Info that is setup once and indept of cartesian layout
///////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<Slave> MPIoffloadEngine::Slaves;
    
int MPIoffloadEngine::UniverseRank;
int MPIoffloadEngine::UniverseSize;

MPI_Comm  MPIoffloadEngine::communicator_universe;
MPI_Comm  MPIoffloadEngine::communicator_cached;
MPI_Group MPIoffloadEngine::WorldGroup;
MPI_Group MPIoffloadEngine::CachedGroup;

MPI_Comm MPIoffloadEngine::HorizontalComm;
int      MPIoffloadEngine::HorizontalRank;
int      MPIoffloadEngine::HorizontalSize;

MPI_Comm MPIoffloadEngine::VerticalComm;
int      MPIoffloadEngine::VerticalSize;
int      MPIoffloadEngine::VerticalRank;
MPI_Win  MPIoffloadEngine::VerticalWindow; 
std::vector<void *>            MPIoffloadEngine::VerticalShmBufs;
std::vector<std::vector<int> > MPIoffloadEngine::UniverseRanks;
std::vector<int>               MPIoffloadEngine::UserCommunicatorToWorldRanks; 

int MPIoffloadEngine::ShmSetup = 0;

void MPIoffloadEngine::CommunicatorInit (MPI_Comm &communicator_world,
					 MPI_Comm &ShmComm,
					 void * &ShmCommBuf)
{      
  int flag;
  assert(ShmSetup==0);  
  
  //////////////////////////////////////////////////////////////////////
  // Universe is all nodes prior to squadron grouping
  //////////////////////////////////////////////////////////////////////
  MPI_Comm_dup (MPI_COMM_WORLD,&communicator_universe);
  MPI_Comm_rank(communicator_universe,&UniverseRank);
  MPI_Comm_size(communicator_universe,&UniverseSize);
  
  /////////////////////////////////////////////////////////////////////
  // Split into groups that can share memory (Verticals)
  /////////////////////////////////////////////////////////////////////
#define MPI_SHARED_MEM_DEBUG
#ifdef  MPI_SHARED_MEM_DEBUG
  MPI_Comm_split(communicator_universe,(UniverseRank/4),UniverseRank,&VerticalComm);
#else 
  MPI_Comm_split_type(communicator_universe, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL,&VerticalComm);
#endif
  MPI_Comm_rank(VerticalComm     ,&VerticalRank);
  MPI_Comm_size(VerticalComm     ,&VerticalSize);
  
  //////////////////////////////////////////////////////////////////////
  // Split into horizontal groups by rank in squadron
  //////////////////////////////////////////////////////////////////////
  MPI_Comm_split(communicator_universe,VerticalRank,UniverseRank,&HorizontalComm);
  MPI_Comm_rank(HorizontalComm,&HorizontalRank);
  MPI_Comm_size(HorizontalComm,&HorizontalSize);
  assert(HorizontalSize*VerticalSize==UniverseSize);
  
  ////////////////////////////////////////////////////////////////////////////////
  // What is my place in the world
  ////////////////////////////////////////////////////////////////////////////////
  int WorldRank=0;
  if(VerticalRank==0) WorldRank = HorizontalRank;
  int ierr=MPI_Allreduce(MPI_IN_PLACE,&WorldRank,1,MPI_INT,MPI_SUM,VerticalComm);
  assert(ierr==0);
  
  ////////////////////////////////////////////////////////////////////////////////
  // Where is the world in the universe?
  ////////////////////////////////////////////////////////////////////////////////
  UniverseRanks = std::vector<std::vector<int> >(HorizontalSize,std::vector<int>(VerticalSize,0));
  UniverseRanks[WorldRank][VerticalRank] = UniverseRank;
  for(int w=0;w<HorizontalSize;w++){
    ierr=MPI_Allreduce(MPI_IN_PLACE,&UniverseRanks[w][0],VerticalSize,MPI_INT,MPI_SUM,communicator_universe);
    assert(ierr==0);
  }
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  // allocate the shared window for our group, pass back Shm info to CartesianCommunicator
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  VerticalShmBufs.resize(VerticalSize);

#undef MPI_SHARED_MEM
#ifdef MPI_SHARED_MEM
  ierr = MPI_Win_allocate_shared(CartesianCommunicator::MAX_MPI_SHM_BYTES,1,MPI_INFO_NULL,VerticalComm,&ShmCommBuf,&VerticalWindow);
  ierr|= MPI_Win_lock_all (MPI_MODE_NOCHECK, VerticalWindow);
  assert(ierr==0);
  //  std::cout<<"SHM "<<ShmCommBuf<<std::endl;

  for(int r=0;r<VerticalSize;r++){
    MPI_Aint sz;
    int dsp_unit;
    MPI_Win_shared_query (VerticalWindow, r, &sz, &dsp_unit, &VerticalShmBufs[r]);
    //    std::cout<<"SHM "<<r<<" " <<VerticalShmBufs[r]<<std::endl;
  }
#else 
  char shm_name [NAME_MAX];
  MPI_Barrier(VerticalComm);

  if ( VerticalRank == 0 ) {
    for(int r=0;r<VerticalSize;r++){

      size_t size = CartesianCommunicator::MAX_MPI_SHM_BYTES;
      if ( r>0 ) size = sizeof(SlaveState);

      sprintf(shm_name,"/Grid_mpi3_shm_%d_%d",WorldRank,r);
      
      shm_unlink(shm_name);

      int fd=shm_open(shm_name,O_RDWR|O_CREAT,0600);
      if ( fd < 0 ) {
	perror("failed shm_open");
	assert(0);
      }

      ftruncate(fd, size);

      VerticalShmBufs[r] = mmap(NULL,size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);

      if ( VerticalShmBufs[r] == MAP_FAILED ) { 
	perror("failed mmap");
	assert(0);
      }

      uint64_t * check = (uint64_t *) VerticalShmBufs[r];
      check[0] = WorldRank;
      check[1] = r;

      //      std::cout<<"SHM "<<r<<" " <<VerticalShmBufs[r]<<std::endl;
    }
  }

  MPI_Barrier(VerticalComm);

  if ( VerticalRank != 0 ) { 
  for(int r=0;r<VerticalSize;r++){

    size_t size = CartesianCommunicator::MAX_MPI_SHM_BYTES ;
    if ( r>0 ) size = sizeof(SlaveState);
    
    sprintf(shm_name,"/Grid_mpi3_shm_%d_%d",WorldRank,r);
    
    int fd=shm_open(shm_name,O_RDWR|O_CREAT,0600);
    if ( fd<0 ) {
      perror("failed shm_open");
      assert(0);
    }
    VerticalShmBufs[r] = mmap(NULL,size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);

    uint64_t * check = (uint64_t *) VerticalShmBufs[r];
    assert(check[0]== WorldRank);
    assert(check[1]== r);
    std::cerr<<"SHM "<<r<<" " <<VerticalShmBufs[r]<<std::endl;
  }
  }
#endif
  MPI_Barrier(VerticalComm);

  //////////////////////////////////////////////////////////////////////
  // Map rank of leader on node in their in new world, to the
  // rank in this vertical plane's horizontal communicator
  //////////////////////////////////////////////////////////////////////
  communicator_world = HorizontalComm;
  ShmComm            = VerticalComm;
  ShmCommBuf         = VerticalShmBufs[0];
  MPI_Comm_group (communicator_world, &WorldGroup); 
  
  ///////////////////////////////////////////////////////////
  // Start the slave data movers
  ///////////////////////////////////////////////////////////
  if ( VerticalRank != 0 ) {
    Slave indentured;
    indentured.Init( (SlaveState *) VerticalShmBufs[VerticalRank], VerticalComm, UniverseRank,VerticalRank);
    indentured.SemInitExcl();// init semaphore in shared memory
    MPI_Barrier(VerticalComm);
    MPI_Barrier(VerticalComm);
    indentured.EventLoop();
    assert(0);
  } else {
    Slaves.resize(VerticalSize);
    for(int i=1;i<VerticalSize;i++){
      Slaves[i].Init((SlaveState *)VerticalShmBufs[i],VerticalComm, UniverseRanks[HorizontalRank][i],i);
    }
    MPI_Barrier(VerticalComm);
    for(int i=1;i<VerticalSize;i++){
      Slaves[i].SemInit();// init semaphore in shared memory
    }
    MPI_Barrier(VerticalComm);
  }
  
  ///////////////////////////////////////////////////////////
  // Verbose for now
  ///////////////////////////////////////////////////////////
  
  ShmSetup=1;
  
  if (UniverseRank == 0){
      
    std::cout<<GridLogMessage << "Grid MPI-3 configuration: detected ";
    std::cout<<UniverseSize   << " Ranks " ;
    std::cout<<HorizontalSize << " Nodes " ;
    std::cout<<VerticalSize   << " with ranks-per-node "<<std::endl;
    
    std::cout<<GridLogMessage << "Grid MPI-3 configuration: using one lead process per node " << std::endl;
    std::cout<<GridLogMessage << "Grid MPI-3 configuration: reduced communicator has size " << HorizontalSize << std::endl;
    
    for(int g=0;g<HorizontalSize;g++){
      std::cout<<GridLogMessage<<" Node "<<g<<" led by MPI rank "<< UniverseRanks[g][0]<<std::endl;
    }
    
    for(int g=0;g<HorizontalSize;g++){
      std::cout<<GridLogMessage<<" { ";
      for(int s=0;s<VerticalSize;s++){
	std::cout<< UniverseRanks[g][s];
	if ( s<VerticalSize-1 ) {
	  std::cout<<",";
	}
      }
      std::cout<<" } "<<std::endl;
    }
  }
};

  ///////////////////////////////////////////////////////////////////////////////////////////////
  // Map the communicator into communicator_world, and find the neighbour.
  // Cache the mappings; cache size is 1.
  ///////////////////////////////////////////////////////////////////////////////////////////////
void MPIoffloadEngine::MapCommRankToWorldRank(int &hashtag, int & comm_world_peer,int tag, MPI_Comm comm,int rank) {

  if ( comm == HorizontalComm ) {
    comm_world_peer = rank;
    //    std::cout << " MapCommRankToWorldRank  horiz " <<rank<<"->"<<comm_world_peer<<std::endl;
  } else if ( comm == communicator_cached ) {
    comm_world_peer = UserCommunicatorToWorldRanks[rank];
    //    std::cout << " MapCommRankToWorldRank  cached " <<rank<<"->"<<comm_world_peer<<std::endl;
  } else { 
    
    int size;

    MPI_Comm_size(comm,&size);

    UserCommunicatorToWorldRanks.resize(size);

    std::vector<int> cached_ranks(size); 

    for(int r=0;r<size;r++) {
      cached_ranks[r]=r;
    }

    communicator_cached=comm;
    
    MPI_Comm_group(communicator_cached, &CachedGroup);
    
    MPI_Group_translate_ranks(CachedGroup,size,&cached_ranks[0],WorldGroup, &UserCommunicatorToWorldRanks[0]); 
    
    comm_world_peer = UserCommunicatorToWorldRanks[rank];
    //    std::cout << " MapCommRankToWorldRank  cache miss " <<rank<<"->"<<comm_world_peer<<std::endl;
    
    assert(comm_world_peer != MPI_UNDEFINED);
  }

  assert( (tag & (~0xFFFFL)) ==0); 
  
  uint64_t icomm = (uint64_t)comm;
  int comm_hash = ((icomm>>0 )&0xFFFF)^((icomm>>16)&0xFFFF)
                ^ ((icomm>>32)&0xFFFF)^((icomm>>48)&0xFFFF);
  
  //  hashtag = (comm_hash<<15) | tag;      
  hashtag = tag;      

};

void Slave::Init(SlaveState * _state,MPI_Comm _squadron,int _universe_rank,int _vertical_rank)
{
  squadron=_squadron;
  universe_rank=_universe_rank;
  vertical_rank=_vertical_rank;
  state   =_state;
  std::cout << "state "<<_state<<" comm "<<_squadron<<" universe_rank"<<universe_rank <<std::endl;
  state->head = state->tail = state->start = 0;
  base = (uint64_t)MPIoffloadEngine::VerticalShmBufs[0];
  int rank; MPI_Comm_rank(_squadron,&rank);
}
#define PERI_PLUS(A) ( (A+1)%pool )
int Slave::Event (void) {

  static int tail_last;
  static int head_last;
  static int start_last;
  int ierr;

  ////////////////////////////////////////////////////
  // Try to advance the start pointers
  ////////////////////////////////////////////////////
  int s=state->start;
  if ( s != state->head ) {
    switch ( state->Descrs[s].command ) {
    case COMMAND_ISEND:
      /*
            std::cout<< " Send "<<s << " ptr "<< state<<" "<< state->Descrs[s].buf<< "["<<state->Descrs[s].bytes<<"]"
      	       << " to " << state->Descrs[s].rank<< " tag" << state->Descrs[s].tag
       << " Comm " << MPIoffloadEngine::communicator_universe<< " me " <<universe_rank<< std::endl;
      */
      ierr = MPI_Isend((void *)(state->Descrs[s].buf+base), 
		       state->Descrs[s].bytes, 
		       MPI_CHAR,
		       state->Descrs[s].rank,
		       state->Descrs[s].tag,
		       MPIoffloadEngine::communicator_universe,
		       (MPI_Request *)&state->Descrs[s].request);
      assert(ierr==0);
      state->start = PERI_PLUS(s);
      return 1;
      break;

    case COMMAND_IRECV:
      /*
      std::cout<< " Recv "<<s << " ptr "<< state<<" "<< state->Descrs[s].buf<< "["<<state->Descrs[s].bytes<<"]"
	       << " from " << state->Descrs[s].rank<< " tag" << state->Descrs[s].tag
	       << " Comm " << MPIoffloadEngine::communicator_universe<< " me "<< universe_rank<< std::endl;
      */
      ierr=MPI_Irecv((void *)(state->Descrs[s].buf+base), 
		     state->Descrs[s].bytes, 
		     MPI_CHAR,
		     state->Descrs[s].rank,
		     state->Descrs[s].tag,
		     MPIoffloadEngine::communicator_universe,
		     (MPI_Request *)&state->Descrs[s].request);

      //      std::cout<< " Request is "<<state->Descrs[s].request<<std::endl;
      //      std::cout<< " Request0 is "<<state->Descrs[0].request<<std::endl;
      assert(ierr==0);
      state->start = PERI_PLUS(s);
      return 1;
      break;

    case COMMAND_WAITALL:

      for(int t=state->tail;t!=s; t=PERI_PLUS(t) ){
	MPI_Wait((MPI_Request *)&state->Descrs[t].request,MPI_STATUS_IGNORE);
      };
      s=PERI_PLUS(s);
      state->start = s;
      state->tail  = s;

      WakeUpCompute();

      return 1;
      break;

    default:
      assert(0);
      break;
    }
  }
  return 0;
}
  //////////////////////////////////////////////////////////////////////////////
  // External interaction with the queue
  //////////////////////////////////////////////////////////////////////////////
  
uint64_t Slave::QueueCommand(int command,void *buf, int bytes, int tag, MPI_Comm comm,int commrank) 
{
  /////////////////////////////////////////
  // Spin; if FIFO is full until not full
  /////////////////////////////////////////
  int head =state->head;
  int next = PERI_PLUS(head);
    
  // Set up descriptor
  int worldrank;
  int hashtag;
  MPI_Comm    communicator;
  MPI_Request request;
  
  MPIoffloadEngine::MapCommRankToWorldRank(hashtag,worldrank,tag,comm,commrank);

  uint64_t relative= (uint64_t)buf - base;
  state->Descrs[head].buf    = relative;
  state->Descrs[head].bytes  = bytes;
  state->Descrs[head].rank   = MPIoffloadEngine::UniverseRanks[worldrank][vertical_rank];
  state->Descrs[head].tag    = hashtag;
  state->Descrs[head].command= command;

  /*  
  if ( command == COMMAND_ISEND ) { 
  std::cout << "QueueSend from "<< universe_rank <<" to commrank " << commrank 
            << " to worldrank " << worldrank <<std::endl;
  std::cout << " via VerticalRank "<< vertical_rank <<" to universerank " << MPIoffloadEngine::UniverseRanks[worldrank][vertical_rank]<<std::endl;
  std::cout << " QueueCommand "<<buf<<"["<<bytes<<"]" << std::endl;
  } 
  if ( command == COMMAND_IRECV ) { 
  std::cout << "QueueRecv on "<< universe_rank <<" from commrank " << commrank 
            << " from worldrank " << worldrank <<std::endl;
  std::cout << " via VerticalRank "<< vertical_rank <<" from universerank " << MPIoffloadEngine::UniverseRanks[worldrank][vertical_rank]<<std::endl;
  std::cout << " QueueSend "<<buf<<"["<<bytes<<"]" << std::endl;
  } 
  */
  // Block until FIFO has space
  while( state->tail==next );

  // Msync on weak order architectures
  // Advance pointer
  state->head = next;

  return 0;
}
  

///////////////////////////////////////////////////////////////////////////////////////////////////
// Info that is setup once and indept of cartesian layout
///////////////////////////////////////////////////////////////////////////////////////////////////

MPI_Comm CartesianCommunicator::communicator_world;

void CartesianCommunicator::Init(int *argc, char ***argv) 
{
  int flag;
  MPI_Initialized(&flag); // needed to coexist with other libs apparently
  if ( !flag ) {
    MPI_Init(argc,argv);
  }
  communicator_world = MPI_COMM_WORLD;
  MPI_Comm ShmComm;
  MPIoffloadEngine::CommunicatorInit (communicator_world,ShmComm,ShmCommBuf);
}
void CartesianCommunicator::ShiftedRanks(int dim,int shift,int &source,int &dest)
{
  int ierr=MPI_Cart_shift(communicator,dim,shift,&source,&dest);
  assert(ierr==0);
}
int CartesianCommunicator::RankFromProcessorCoor(std::vector<int> &coor)
{
  int rank;
  int ierr=MPI_Cart_rank  (communicator, &coor[0], &rank);
  assert(ierr==0);
  return rank;
}
void  CartesianCommunicator::ProcessorCoorFromRank(int rank, std::vector<int> &coor)
{
  coor.resize(_ndimension);
  int ierr=MPI_Cart_coords  (communicator, rank, _ndimension,&coor[0]);
  assert(ierr==0);
}

CartesianCommunicator::CartesianCommunicator(const std::vector<int> &processors)
{ 
  _ndimension = processors.size();
  std::vector<int> periodic(_ndimension,1);

  _Nprocessors=1;
  _processors = processors;

  for(int i=0;i<_ndimension;i++){
    _Nprocessors*=_processors[i];
  }

  int Size; 
  MPI_Comm_size(communicator_world,&Size);
  assert(Size==_Nprocessors);

  _processor_coor.resize(_ndimension);
  MPI_Cart_create(communicator_world, _ndimension,&_processors[0],&periodic[0],1,&communicator);
  MPI_Comm_rank  (communicator,&_processor);
  MPI_Cart_coords(communicator,_processor,_ndimension,&_processor_coor[0]);
};

void CartesianCommunicator::GlobalSum(uint32_t &u){
  int ierr=MPI_Allreduce(MPI_IN_PLACE,&u,1,MPI_UINT32_T,MPI_SUM,communicator);
  assert(ierr==0);
}
void CartesianCommunicator::GlobalSum(uint64_t &u){
  int ierr=MPI_Allreduce(MPI_IN_PLACE,&u,1,MPI_UINT64_T,MPI_SUM,communicator);
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

// Basic Halo comms primitive
void CartesianCommunicator::SendToRecvFrom(void *xmit,
					   int dest,
					   void *recv,
					   int from,
					   int bytes)
{
  std::vector<CommsRequest_t> reqs(0);
  SendToRecvFromBegin(reqs,xmit,dest,recv,from,bytes);
  SendToRecvFromComplete(reqs);
}

void CartesianCommunicator::SendRecvPacket(void *xmit,
					   void *recv,
					   int sender,
					   int receiver,
					   int bytes)
{
  MPI_Status stat;
  assert(sender != receiver);
  int tag = sender;
  if ( _processor == sender ) {
    MPI_Send(xmit, bytes, MPI_CHAR,receiver,tag,communicator);
  }
  if ( _processor == receiver ) { 
    MPI_Recv(recv, bytes, MPI_CHAR,sender,tag,communicator,&stat);
  }
}

// Basic Halo comms primitive
void CartesianCommunicator::SendToRecvFromBegin(std::vector<CommsRequest_t> &list,
						void *xmit,
						int dest,
						void *recv,
						int from,
						int bytes)
{
  MPI_Request xrq;
  MPI_Request rrq;
  int rank = _processor;
  int ierr;
  ierr =MPI_Isend(xmit, bytes, MPI_CHAR,dest,_processor,communicator,&xrq);
  ierr|=MPI_Irecv(recv, bytes, MPI_CHAR,from,from,communicator,&rrq);
  
  assert(ierr==0);

  list.push_back(xrq);
  list.push_back(rrq);
}

void CartesianCommunicator::StencilSendToRecvFromBegin(std::vector<CommsRequest_t> &list,
						       void *xmit,
						       int dest,
						       void *recv,
						       int from,
						       int bytes)
{
  uint64_t xmit_i = (uint64_t) xmit;
  uint64_t recv_i = (uint64_t) recv;
  uint64_t shm    = (uint64_t) ShmCommBuf;
  // assert xmit and recv lie in shared memory region
  assert( (xmit_i >= shm) && (xmit_i+bytes <= shm+MAX_MPI_SHM_BYTES) );
  assert( (recv_i >= shm) && (recv_i+bytes <= shm+MAX_MPI_SHM_BYTES) );
  assert(from!=_processor);
  assert(dest!=_processor);
  MPIoffloadEngine::QueueMultiplexedSend(xmit,bytes,_processor,communicator,dest);
  MPIoffloadEngine::QueueMultiplexedRecv(recv,bytes,from,communicator,from);
}


void CartesianCommunicator::StencilSendToRecvFromComplete(std::vector<CommsRequest_t> &list)
{
  MPIoffloadEngine::WaitAll();
}

void CartesianCommunicator::StencilBarrier(void)
{
}

void CartesianCommunicator::SendToRecvFromComplete(std::vector<CommsRequest_t> &list)
{
  int nreq=list.size();
  std::vector<MPI_Status> status(nreq);
  int ierr = MPI_Waitall(nreq,&list[0],&status[0]);
  assert(ierr==0);
}

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

void CartesianCommunicator::BroadcastWorld(int root,void* data, int bytes)
{
  int ierr= MPI_Bcast(data,
		      bytes,
		      MPI_BYTE,
		      root,
		      communicator_world);
  assert(ierr==0);
}

void *CartesianCommunicator::ShmBufferSelf(void) { return ShmCommBuf; }

void *CartesianCommunicator::ShmBuffer(int rank) {
  return NULL;
}
void *CartesianCommunicator::ShmBufferTranslate(int rank,void * local_p) { 
  return NULL;
}


};

