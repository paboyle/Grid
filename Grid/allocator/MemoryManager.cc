#include <Grid/GridCore.h>

NAMESPACE_BEGIN(Grid);

/*Allocation types, saying which pointer cache should be used*/
#define Cpu      (0)
#define CpuSmall (1)
#define Acc      (2)
#define AccSmall (3)

//////////////////////////////////////////////////////////////////////
// Data tables for recently freed pooiniter caches
//////////////////////////////////////////////////////////////////////
MemoryManager::AllocationCacheEntry MemoryManager::Entries[MemoryManager::NallocType][MemoryManager::NallocCacheMax];
int MemoryManager::Victim[MemoryManager::NallocType];
int MemoryManager::Ncache[MemoryManager::NallocType];

//////////////////////////////////////////////////////////////////////
// Actual allocation and deallocation utils
//////////////////////////////////////////////////////////////////////
void *MemoryManager::AcceleratorAllocate(size_t bytes)
{
  void *ptr = (void *) Lookup(bytes,Acc);

  if ( ptr == (void *) NULL ) {
    ptr = (void *) acceleratorAllocDevice(bytes);
    //    std::cout <<"AcceleratorAllocate: allocated Accelerator pointer "<<std::hex<<ptr<<std::endl;
  }

  return ptr;
}
void  MemoryManager::AcceleratorFree    (void *ptr,size_t bytes)
{
  void *__freeme = Insert(ptr,bytes,Acc);

  if ( __freeme ) acceleratorFreeDevice(__freeme);
}
void *MemoryManager::CpuAllocate(size_t bytes)
{
  void *ptr = (void *) Lookup(bytes,Cpu);

  if ( ptr == (void *) NULL ) {
    ptr = (void *) acceleratorAllocShared(bytes);
    //    std::cout <<"CpuAllocate: allocated Cpu pointer "<<std::hex<<ptr<<std::endl;
  }

  return ptr;
}
void  MemoryManager::CpuFree    (void *_ptr,size_t bytes)
{
  NotifyDeletion(_ptr);

  // If present remove entry and free accelerator too.
  // Can we ever hit a free event with a view still in scope?
  void *__freeme = Insert(_ptr,bytes,Cpu);
  if ( __freeme ) acceleratorFreeShared(__freeme);
}
//////////////////////////////////////////
// call only once
//////////////////////////////////////////
void MemoryManager::Init(void)
{
  Ncache[Cpu] = 8;
  Ncache[Acc] = 8;
  Ncache[CpuSmall] = 32;
  Ncache[AccSmall] = 32;

  char * str;
  int Nc;
  int NcS;
  
  str= getenv("GRID_ALLOC_NCACHE_LARGE");
  if ( str ) {
    Nc = atoi(str);
    if ( (Nc>=0) && (Nc < NallocCacheMax)) {
      Ncache[Cpu]=Nc;
      Ncache[Acc]=Nc;
    }
  }

  str= getenv("GRID_ALLOC_NCACHE_SMALL");
  if ( str ) {
    Nc = atoi(str);
    if ( (Nc>=0) && (Nc < NallocCacheMax)) {
      Ncache[CpuSmall]=Nc;
      Ncache[AccSmall]=Nc;
    }
  }
  std::cout << GridLogMessage<< "MemoryManager::Init() setting up"<<std::endl;
#ifdef ALLOCATION_CACHE
  std::cout << GridLogMessage<< "MemoryManager::Init() cache pool for recent allocations: SMALL "<<Ncache[CpuSmall]<<" LARGE "<<Ncache[Cpu]<<std::endl;
#endif

#ifdef GRID_UVM
  std::cout << GridLogMessage<< "MemoryManager::Init() Unified memory space"<<std::endl;
#ifdef GRID_CUDA
  std::cout << GridLogMessage<< "MemoryManager::Init() Using cudaMallocManaged"<<std::endl;
#endif
#ifdef GRID_HIP
  std::cout << GridLogMessage<< "MemoryManager::Init() Using hipMallocManaged"<<std::endl;
#endif
#ifdef GRID_SYCL
  std::cout << GridLogMessage<< "MemoryManager::Init() Using SYCL malloc_shared"<<std::endl;
#endif
#else
  std::cout << GridLogMessage<< "MemoryManager::Init() Non unified: Caching accelerator data in dedicated memory"<<std::endl;
#ifdef GRID_CUDA
  std::cout << GridLogMessage<< "MemoryManager::Init() Using cudaMalloc"<<std::endl;
#endif
#ifdef GRID_HIP
  std::cout << GridLogMessage<< "MemoryManager::Init() Using hipMalloc"<<std::endl;
#endif
#ifdef GRID_SYCL
  std::cout << GridLogMessage<< "MemoryManager::Init() Using SYCL malloc_device"<<std::endl;
#endif
#endif
}

void *MemoryManager::Insert(void *ptr,size_t bytes,int type) 
{
#ifdef ALLOCATION_CACHE
  bool small = (bytes < GRID_ALLOC_SMALL_LIMIT);
  int cache = type + small;
  return Insert(ptr,bytes,Entries[cache],Ncache[cache],Victim[cache]);  
#else
  return ptr;
#endif
}

void *MemoryManager::Insert(void *ptr,size_t bytes,AllocationCacheEntry *entries,int ncache,int &victim) 
{
  assert(ncache>0);
#ifdef GRID_OMP
  assert(omp_in_parallel()==0);
#endif 

  void * ret = NULL;
  int v = -1;

  for(int e=0;e<ncache;e++) {
    if ( entries[e].valid==0 ) {
      v=e; 
      break;
    }
  }

  if ( v==-1 ) {
    v=victim;
    victim = (victim+1)%ncache;
  }

  if ( entries[v].valid ) {
    ret = entries[v].address;
    entries[v].valid = 0;
    entries[v].address = NULL;
    entries[v].bytes = 0;
  }

  entries[v].address=ptr;
  entries[v].bytes  =bytes;
  entries[v].valid  =1;

  return ret;
}

void *MemoryManager::Lookup(size_t bytes,int type)
{
#ifdef ALLOCATION_CACHE
  bool small = (bytes < GRID_ALLOC_SMALL_LIMIT);
  int cache = type+small;
  return Lookup(bytes,Entries[cache],Ncache[cache]);
#else
  return NULL;
#endif
}

void *MemoryManager::Lookup(size_t bytes,AllocationCacheEntry *entries,int ncache) 
{
  assert(ncache>0);
#ifdef GRID_OMP
  assert(omp_in_parallel()==0);
#endif 
  for(int e=0;e<ncache;e++){
    if ( entries[e].valid && ( entries[e].bytes == bytes ) ) {
      entries[e].valid = 0;
      return entries[e].address;
    }
  }
  return NULL;
}


NAMESPACE_END(Grid);

