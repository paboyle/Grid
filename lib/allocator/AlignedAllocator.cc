#include <Grid/GridCore.h>
#include <fcntl.h>

namespace Grid {

MemoryStats *MemoryProfiler::stats = nullptr;
bool         MemoryProfiler::debug = false;

int PointerCache::victim;

PointerCache::PointerCacheEntry PointerCache::Entries[PointerCache::Ncache];

void *PointerCache::Insert(void *ptr,size_t bytes) {

  if (bytes < 4096 ) return ptr;

#ifdef GRID_OMP
  assert(omp_in_parallel()==0);
#endif 

  void * ret = NULL;
  int v = -1;

  for(int e=0;e<Ncache;e++) {
    if ( Entries[e].valid==0 ) {
      v=e; 
      break;
    }
  }

  if ( v==-1 ) {
    v=victim;
    victim = (victim+1)%Ncache;
  }

  if ( Entries[v].valid ) {
    ret = Entries[v].address;
    Entries[v].valid = 0;
    Entries[v].address = NULL;
    Entries[v].bytes = 0;
  }

  Entries[v].address=ptr;
  Entries[v].bytes  =bytes;
  Entries[v].valid  =1;

  return ret;
}

void *PointerCache::Lookup(size_t bytes) {

 if (bytes < 4096 ) return NULL;

#ifdef _OPENMP
  assert(omp_in_parallel()==0);
#endif 

  for(int e=0;e<Ncache;e++){
    if ( Entries[e].valid && ( Entries[e].bytes == bytes ) ) {
      Entries[e].valid = 0;
      return Entries[e].address;
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

}
