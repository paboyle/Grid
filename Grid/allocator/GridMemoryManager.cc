#include <Grid/GridCore.h>

NAMESPACE_BEGIN(Grid);

#define _GRID_MEM_PAGE_SIZE 4096
void* _grid_mem_base = 0;
size_t _grid_mem_pages;
struct _grid_mem_range {
  size_t page_start, page_end;
};
std::vector<_grid_mem_range> _grid_mem_avail;
std::map<void*,_grid_mem_range> _grid_mem_alloc;

void gridMemoryInit() {
#ifdef GRID_NVCC
  size_t free,total;
  cudaMemGetInfo(&free,&total);
  
  char* ev = getenv("GRID_DEVICE_BYTES_FOR_CACHE");
  if (ev) {
    long bytes;
    assert(sscanf(ev,"%ld",&bytes)==1);
    free -= bytes;
  }

  _grid_mem_pages = free / _GRID_MEM_PAGE_SIZE;
  size_t sz = _grid_mem_pages * _GRID_MEM_PAGE_SIZE;

  assert(cudaSuccess == cudaMallocManaged(&_grid_mem_base,sz));
  
  int target;
  cudaGetDevice(&target);
  cudaMemAdvise(_grid_mem_base,sz,cudaMemAdviseSetPreferredLocation,target);

  assert(cudaSuccess == cudaMemset(_grid_mem_base,0,sz)); // touch on device
  std::cout << GridLogMessage << "gridMemoryInit: " << sz << " bytes" << std::endl;

  _grid_mem_avail.push_back( { 0, _grid_mem_pages } );
#endif
}

void gridMallocManaged(void** pp, size_t sz) {
#ifdef GRID_NVCC
  if (_grid_mem_avail.empty())
    gridMemoryInit();

  size_t pages = (sz + _GRID_MEM_PAGE_SIZE - 1) / _GRID_MEM_PAGE_SIZE;
  // find free block
  size_t m;
  for (m=0;m<_grid_mem_avail.size();m++) {
    auto & b = _grid_mem_avail[m];
    if (b.page_end - b.page_start >= pages)
      break;
  }
  if (m == _grid_mem_avail.size()) {
    std::cout << GridLogMessage << "Out of memory" << std::endl;
    assert(0);
  }
  *pp = (char*)_grid_mem_base + _GRID_MEM_PAGE_SIZE*_grid_mem_avail[m].page_start;
  _grid_mem_alloc[*pp] = { _grid_mem_avail[m].page_start, _grid_mem_avail[m].page_start + pages };
  _grid_mem_avail[m].page_start += pages;
#else
  *pp = malloc(sz);
#endif
}

void gridFree(void* p) {
#ifdef GRID_NVCC
  if (_grid_mem_avail.empty())
    gridMemoryInit();

  auto & alloc = _grid_mem_alloc[p];
  if (alloc.page_start == alloc.page_end) {
    free(p);
    //cudaFreeHost(p);
  } else {
    // can we enlarge existing one?
    for (size_t m=0;m<_grid_mem_avail.size();m++) {
      auto & b = _grid_mem_avail[m];
      if (b.page_start == alloc.page_end) {
	b.page_start = alloc.page_start;
	return;
      }
      if (b.page_end == alloc.page_start) {
	b.page_end = alloc.page_end;
	return;
      }
    }
    // fragment memory
    _grid_mem_avail.push_back( alloc );  
  }
  _grid_mem_alloc.erase(p);
#else
  free(p);
#endif
}

void gridAcceleratorPrefetch(void* p, size_t sz) {
#ifdef GRID_NVCC
  auto & alloc = _grid_mem_alloc[p];
  if (alloc.page_start == alloc.page_end) // pinned to host
    return;

  int target;
  cudaGetDevice(&target);
  cudaMemPrefetchAsync(p,sz,target);
#endif
}

void gridMemGetInfo(size_t* pfree, size_t* ptotal) {
#ifdef GRID_NVCC
  if (_grid_mem_avail.empty())
    gridMemoryInit();

  *ptotal = _grid_mem_pages * _GRID_MEM_PAGE_SIZE;
  *pfree = 0;
  for (auto & a : _grid_mem_avail)
    *pfree += (a.page_end - a.page_start) * _GRID_MEM_PAGE_SIZE;
#else
  *pfree = 0;
  *ptotal = 0;
#endif
}

void gridMoveToHost(void** pp) {
#ifdef GRID_NVCC
  if (_grid_mem_avail.empty())
    gridMemoryInit();

  auto & alloc = _grid_mem_alloc[*pp];
  if (alloc.page_start == alloc.page_end) // already on host
    return;

  size_t sz = (alloc.page_end - alloc.page_start) * _GRID_MEM_PAGE_SIZE;
  void*pn;
  //assert(cudaSuccess == cudaMallocHost(&pn,sz));
  pn = malloc(sz);
  memcpy(pn,*pp,sz);
  gridFree(*pp);
  *pp = pn;
  _grid_mem_alloc[pn] = { 0,0 };
#endif
}

NAMESPACE_END(Grid);
