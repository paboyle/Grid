#include <Grid/GridCore.h>
#include <fcntl.h>

NAMESPACE_BEGIN(Grid);

MemoryStats *MemoryProfiler::stats = nullptr;
bool         MemoryProfiler::debug = false;

int PointerCache::NcacheSmall = PointerCache::NcacheSmallMax;
#ifdef GRID_CUDA
int PointerCache::Ncache      = 32;
#else 
int PointerCache::Ncache      = 8;
#endif
int PointerCache::Victim;
int PointerCache::VictimSmall;
PointerCache::PointerCacheEntry PointerCache::Entries[PointerCache::NcacheMax];
PointerCache::PointerCacheEntry PointerCache::EntriesSmall[PointerCache::NcacheSmallMax];

void PointerCache::Init(void)
{
  char * str;

  str= getenv("GRID_ALLOC_NCACHE_LARGE");
  if ( str ) Ncache = atoi(str);
  if ( (Ncache<0) || (Ncache > NcacheMax)) Ncache = NcacheMax;

  str= getenv("GRID_ALLOC_NCACHE_SMALL");
  if ( str ) NcacheSmall = atoi(str);
  if ( (NcacheSmall<0) || (NcacheSmall > NcacheSmallMax)) NcacheSmall = NcacheSmallMax;

  //  printf("Aligned alloocator cache: large %d/%d small %d/%d\n",Ncache,NcacheMax,NcacheSmall,NcacheSmallMax);
}
void *PointerCache::Insert(void *ptr,size_t bytes) 
{
  if (bytes < GRID_ALLOC_SMALL_LIMIT ) 
    return Insert(ptr,bytes,EntriesSmall,NcacheSmall,VictimSmall);
  return Insert(ptr,bytes,Entries,Ncache,Victim);  
}
void *PointerCache::Insert(void *ptr,size_t bytes,PointerCacheEntry *entries,int ncache,int &victim) 
{
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

void *PointerCache::Lookup(size_t bytes)
{
  if (bytes < GRID_ALLOC_SMALL_LIMIT ) 
    return Lookup(bytes,EntriesSmall,NcacheSmall);
  return Lookup(bytes,Entries,Ncache);
}
void *PointerCache::Lookup(size_t bytes,PointerCacheEntry *entries,int ncache) 
{
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


void check_huge_pages(void *Buf,uint64_t BYTES)
{
#ifdef __linux__
  int fd = open("/proc/self/pagemap", O_RDONLY);
  assert(fd >= 0);
  const int page_size = 4096;
  uint64_t virt_pfn = (uint64_t)Buf / page_size;
  off_t offset = sizeof(uint64_t) * virt_pfn;
  uint64_t npages = (BYTES + page_size-1) / page_size;
  uint64_t pagedata[npages];
  uint64_t ret = lseek(fd, offset, SEEK_SET);
  assert(ret == offset);
  ret = ::read(fd, pagedata, sizeof(uint64_t)*npages);
  assert(ret == sizeof(uint64_t) * npages);
  int nhugepages = npages / 512;
  int n4ktotal, nnothuge;
  n4ktotal = 0;
  nnothuge = 0;
  for (int i = 0; i < nhugepages; ++i) {
    uint64_t baseaddr = (pagedata[i*512] & 0x7fffffffffffffULL) * page_size;
    for (int j = 0; j < 512; ++j) {
      uint64_t pageaddr = (pagedata[i*512+j] & 0x7fffffffffffffULL) * page_size;
      ++n4ktotal;
      if (pageaddr != baseaddr + j * page_size)
	++nnothuge;
    }
  }
  int rank = CartesianCommunicator::RankWorld();
  printf("rank %d Allocated %d 4k pages, %d not in huge pages\n", rank, n4ktotal, nnothuge);
#endif
}

std::string sizeString(const size_t bytes)
{
  constexpr unsigned int bufSize = 256;
  const char             *suffixes[7] = {"", "K", "M", "G", "T", "P", "E"};
  char                   buf[256];
  size_t                 s     = 0;
  double                 count = bytes;
  
  while (count >= 1024 && s < 7)
    {
      s++;
      count /= 1024;
    }
  if (count - floor(count) == 0.0)
    {
      snprintf(buf, bufSize, "%d %sB", (int)count, suffixes[s]);
    }
  else
    {
      snprintf(buf, bufSize, "%.1f %sB", count, suffixes[s]);
    }
  
  return std::string(buf);
}

NAMESPACE_END(Grid);

