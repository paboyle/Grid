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
AllocationCache::AllocationCacheEntry AllocationCache::Entries[AllocationCache::NallocType][AllocationCache::NallocCacheMax];
int AllocationCache::Victim[AllocationCache::NallocType];
int AllocationCache::Ncache[AllocationCache::NallocType];

//////////////////////////////////////////////////////////////////////
// Actual allocation and deallocation utils
//////////////////////////////////////////////////////////////////////
void *AllocationCache::AcceleratorAllocate(size_t bytes)
{
  void *ptr = (void *) Lookup(bytes,Acc);

  if ( ptr == (void *) NULL ) {
    ptr = (void *) acceleratorAllocDevice(bytes);
    //    std::cout <<"AcceleratorAllocate: allocated Accelerator pointer "<<std::hex<<ptr<<std::endl;
  }

  return ptr;
}
void  AllocationCache::AcceleratorFree    (void *ptr,size_t bytes)
{
  void *__freeme = Insert(ptr,bytes,Acc);

  if ( __freeme ) acceleratorFreeDevice(__freeme);
}
void *AllocationCache::CpuAllocate(size_t bytes)
{
  void *ptr = (void *) Lookup(bytes,Cpu);

  if ( ptr == (void *) NULL ) {
    ptr = (void *) acceleratorAllocShared(bytes);
    //    std::cout <<"CpuAllocate: allocated Cpu pointer "<<std::hex<<ptr<<std::endl;
  }

  return ptr;
}
void  AllocationCache::CpuFree    (void *ptr,size_t bytes)
{
  // Look up in ViewCache
  int e=CpuViewLookup(ptr);
  if(e>=0){ Discard(e); }

  // If present remove entry and free accelerator too.
  // Can we ever hit a free event with a view still in scope?
  void *__freeme = Insert(ptr,bytes,Cpu);
  //  std::cout <<"CpuFree cached pointer "<<std::hex<<ptr<<std::endl;
  //  std::cout <<"CpuFree deallocating pointer "<<std::hex<<__freeme<<std::endl;
  if ( __freeme ) acceleratorFreeShared(__freeme);
}
//////////////////////////////////////////
// call only once
//////////////////////////////////////////
void AllocationCache::Init(void)
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
  std::cout << "MemoryManager::Init() SMALL "<<Ncache[CpuSmall]<<" LARGE "<<Ncache[Cpu]<<std::endl;
}

void *AllocationCache::Insert(void *ptr,size_t bytes,int type) 
{
#ifdef ALLOCATION_CACHE
  bool small = (bytes < GRID_ALLOC_SMALL_LIMIT);
  int cache = type + small;
  return Insert(ptr,bytes,Entries[cache],Ncache[cache],Victim[cache]);  
#else
  return ptr;
#endif
}
void *AllocationCache::Insert(void *ptr,size_t bytes,AllocationCacheEntry *entries,int ncache,int &victim) 
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

void *AllocationCache::Lookup(size_t bytes,int type)
{
#ifdef ALLOCATION_CACHE
  bool small = (bytes < GRID_ALLOC_SMALL_LIMIT);
  int cache = type+small;
  return Lookup(bytes,Entries[cache],Ncache[cache]);
#else
  return NULL;
#endif
}
void *AllocationCache::Lookup(size_t bytes,AllocationCacheEntry *entries,int ncache) 
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

