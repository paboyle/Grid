#include <Grid/GridCore.h>

NAMESPACE_BEGIN(Grid);

/*Allocation types, saying which pointer cache should be used*/
#define Cpu      (0)
#define CpuSmall (1)
#define Acc      (2)
#define AccSmall (3)
#define Shared   (4)
#define SharedSmall (5)
#undef GRID_MM_VERBOSE 
uint64_t total_shared;
uint64_t total_device;
uint64_t total_host;;
void MemoryManager::PrintBytes(void)
{
  std::cout << " MemoryManager : ------------------------------------ "<<std::endl;
  std::cout << " MemoryManager : PrintBytes "<<std::endl;
  std::cout << " MemoryManager : ------------------------------------ "<<std::endl;
  std::cout << " MemoryManager : "<<(total_shared>>20)<<" shared      Mbytes "<<std::endl;
  std::cout << " MemoryManager : "<<(total_device>>20)<<" accelerator Mbytes "<<std::endl;
  std::cout << " MemoryManager : "<<(total_host>>20)  <<" cpu         Mbytes "<<std::endl;
  uint64_t cacheBytes;
  cacheBytes = CacheBytes[Cpu];
  std::cout << " MemoryManager : "<<(cacheBytes>>20) <<" cpu cache Mbytes "<<std::endl;
  cacheBytes = CacheBytes[Acc];
  std::cout << " MemoryManager : "<<(cacheBytes>>20) <<" acc cache Mbytes "<<std::endl;
  cacheBytes = CacheBytes[Shared];
  std::cout << " MemoryManager : "<<(cacheBytes>>20) <<" shared cache Mbytes "<<std::endl;
  
#ifdef GRID_CUDA
  cuda_mem();
#endif
  
}

//////////////////////////////////////////////////////////////////////
// Data tables for recently freed pooiniter caches
//////////////////////////////////////////////////////////////////////
MemoryManager::AllocationCacheEntry MemoryManager::Entries[MemoryManager::NallocType][MemoryManager::NallocCacheMax];
int MemoryManager::Victim[MemoryManager::NallocType];
int MemoryManager::Ncache[MemoryManager::NallocType] = { 2, 8, 8, 16, 8, 16 };
uint64_t MemoryManager::CacheBytes[MemoryManager::NallocType];
//////////////////////////////////////////////////////////////////////
// Actual allocation and deallocation utils
//////////////////////////////////////////////////////////////////////
void *MemoryManager::AcceleratorAllocate(size_t bytes)
{
  total_device+=bytes;
  void *ptr = (void *) Lookup(bytes,Acc);
  if ( ptr == (void *) NULL ) {
    ptr = (void *) acceleratorAllocDevice(bytes);
  }
#ifdef GRID_MM_VERBOSE
  std::cout <<"AcceleratorAllocate "<<std::endl;
  PrintBytes();
#endif
  return ptr;
}
void  MemoryManager::AcceleratorFree    (void *ptr,size_t bytes)
{
  total_device-=bytes;
  void *__freeme = Insert(ptr,bytes,Acc);
  if ( __freeme ) {
    acceleratorFreeDevice(__freeme);
  }
#ifdef GRID_MM_VERBOSE
  std::cout <<"AcceleratorFree "<<std::endl;
  PrintBytes();
#endif
}
void *MemoryManager::SharedAllocate(size_t bytes)
{
  total_shared+=bytes;
  void *ptr = (void *) Lookup(bytes,Shared);
  if ( ptr == (void *) NULL ) {
    ptr = (void *) acceleratorAllocShared(bytes);
  }
#ifdef GRID_MM_VERBOSE
  std::cout <<"SharedAllocate "<<std::endl;
  PrintBytes();
#endif
  return ptr;
}
void  MemoryManager::SharedFree    (void *ptr,size_t bytes)
{
  total_shared-=bytes;
  void *__freeme = Insert(ptr,bytes,Shared);
  if ( __freeme ) {
    acceleratorFreeShared(__freeme);
  }
#ifdef GRID_MM_VERBOSE
  std::cout <<"SharedFree "<<std::endl;
  PrintBytes();
#endif
}
#ifdef GRID_UVM
void *MemoryManager::CpuAllocate(size_t bytes)
{
  total_host+=bytes;
  void *ptr = (void *) Lookup(bytes,Cpu);
  if ( ptr == (void *) NULL ) {
    ptr = (void *) acceleratorAllocShared(bytes);
  }
#ifdef GRID_MM_VERBOSE
  std::cout <<"CpuAllocate "<<std::endl;
  PrintBytes();
#endif
  return ptr;
}
void  MemoryManager::CpuFree    (void *_ptr,size_t bytes)
{
  total_host-=bytes;
  NotifyDeletion(_ptr);
  void *__freeme = Insert(_ptr,bytes,Cpu);
  if ( __freeme ) { 
    acceleratorFreeShared(__freeme);
  }
#ifdef GRID_MM_VERBOSE
  std::cout <<"CpuFree "<<std::endl;
  PrintBytes();
#endif
}
#else
void *MemoryManager::CpuAllocate(size_t bytes)
{
  total_host+=bytes;
  void *ptr = (void *) Lookup(bytes,Cpu);
  if ( ptr == (void *) NULL ) {
    ptr = (void *) acceleratorAllocCpu(bytes);
  }
#ifdef GRID_MM_VERBOSE
  std::cout <<"CpuAllocate "<<std::endl;
  PrintBytes();
#endif
  return ptr;
}
void  MemoryManager::CpuFree    (void *_ptr,size_t bytes)
{
  total_host-=bytes;
  NotifyDeletion(_ptr);
  void *__freeme = Insert(_ptr,bytes,Cpu);
  if ( __freeme ) { 
    acceleratorFreeCpu(__freeme);
  }
#ifdef GRID_MM_VERBOSE
  std::cout <<"CpuFree "<<std::endl;
  PrintBytes();
#endif
}
#endif

//////////////////////////////////////////
// call only once
//////////////////////////////////////////
void MemoryManager::Init(void)
{

  char * str;
  int Nc;
  
  str= getenv("GRID_ALLOC_NCACHE_LARGE");
  if ( str ) {
    Nc = atoi(str);
    if ( (Nc>=0) && (Nc < NallocCacheMax)) {
      Ncache[Cpu]=Nc;
      Ncache[Acc]=Nc;
      Ncache[Shared]=Nc;
    }
  }

  str= getenv("GRID_ALLOC_NCACHE_SMALL");
  if ( str ) {
    Nc = atoi(str);
    if ( (Nc>=0) && (Nc < NallocCacheMax)) {
      Ncache[CpuSmall]=Nc;
      Ncache[AccSmall]=Nc;
      Ncache[SharedSmall]=Nc;
    }
  }

}

void MemoryManager::InitMessage(void) {

#ifndef GRID_UVM
  std::cout << GridLogMessage << "MemoryManager Cache "<< MemoryManager::DeviceMaxBytes <<" bytes "<<std::endl;
#endif
  
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
  return Insert(ptr,bytes,Entries[cache],Ncache[cache],Victim[cache],CacheBytes[cache]);  
#else
  return ptr;
#endif
}

void *MemoryManager::Insert(void *ptr,size_t bytes,AllocationCacheEntry *entries,int ncache,int &victim, uint64_t &cacheBytes) 
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
    cacheBytes -= entries[v].bytes;
    entries[v].valid = 0;
    entries[v].address = NULL;
    entries[v].bytes = 0;
  }

  entries[v].address=ptr;
  entries[v].bytes  =bytes;
  entries[v].valid  =1;
  cacheBytes += bytes;

  return ret;
}

void *MemoryManager::Lookup(size_t bytes,int type)
{
#ifdef ALLOCATION_CACHE
  bool small = (bytes < GRID_ALLOC_SMALL_LIMIT);
  int cache = type+small;
  return Lookup(bytes,Entries[cache],Ncache[cache],CacheBytes[cache]);
#else
  return NULL;
#endif
}

void *MemoryManager::Lookup(size_t bytes,AllocationCacheEntry *entries,int ncache,uint64_t & cacheBytes) 
{
  assert(ncache>0);
#ifdef GRID_OMP
  assert(omp_in_parallel()==0);
#endif 
  for(int e=0;e<ncache;e++){
    if ( entries[e].valid && ( entries[e].bytes == bytes ) ) {
      entries[e].valid = 0;
      cacheBytes -= entries[e].bytes;
      return entries[e].address;
    }
  }
  return NULL;
}


NAMESPACE_END(Grid);

