#include <Grid/GridCore.h>
#ifdef GRID_UVM

#warning "Grid is assuming unified virtual memory address space"
NAMESPACE_BEGIN(Grid);
/////////////////////////////////////////////////////////////////////////////////
// View management is 1:1 address space mapping
/////////////////////////////////////////////////////////////////////////////////

void  AllocationCache::AcceleratorViewClose(void* AccPtr){};
void *AllocationCache::AcceleratorViewOpen(void* CpuPtr,size_t bytes,ViewMode mode,ViewAdvise hint){ return CpuPtr; }
void  AllocationCache::CpuViewClose(void* Ptr){};
void *AllocationCache::CpuViewOpen(void* CpuPtr,size_t bytes,ViewMode mode,ViewAdvise hint){ return CpuPtr; }
int   AllocationCache::CpuViewLookup(void *CpuPtr){ return 0;}
/////////////////////////////////////
// Dummy stubs
/////////////////////////////////////
void  AllocationCache::CpuDiscard(int e)      { return;}
void  AllocationCache::Discard(int e)      { return;}
void  AllocationCache::Evict(int e)      { return; }
void  AllocationCache::Flush(int e)      { assert(0);}
void  AllocationCache::Clone(int e)      { assert(0);}
int   AllocationCache::ViewVictim(void)  { assert(0); return 0;}
void  AllocationCache::ViewClose(void* AccPtr,ViewMode mode){};
void *AllocationCache::ViewOpen (void* CpuPtr,size_t bytes,ViewMode mode,ViewAdvise hint){return CpuPtr;};

NAMESPACE_END(Grid);
#endif
