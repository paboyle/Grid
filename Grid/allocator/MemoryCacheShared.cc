#include <Grid/GridCore.h>
#ifdef GRID_UNIFIED

#warning "Grid is assuming unified virtual memory address space"
NAMESPACE_BEGIN(Grid);
/////////////////////////////////////////////////////////////////////////////////
// View management is 1:1 address space mapping
/////////////////////////////////////////////////////////////////////////////////

void *AllocationCache::CpuViewOpen(void* CpuPtr,size_t bytes,int mode,int transient) { return CpuPtr; }
void *AllocationCache::AccViewOpen(void* CpuPtr,size_t bytes,int mode,int transient) { return CpuPtr; }
void  AllocationCache::AccViewClose(void* AccPtr){}
void  AllocationCache::CpuViewClose(void* CpuPtr){}

/////////////////////////////////////
// Dummy stubs
/////////////////////////////////////
int  AllocationCache::ViewVictim(void)  { assert(0); return 0;}
void AllocationCache::Evict(int e)      { assert(0);}
void AllocationCache::Flush(int e)      { assert(0);}
void AllocationCache::Clone(int e)      { assert(0);}

int   AllocationCache::CpuViewLookup(void *CpuPtr){assert(0); return 0;}
int   AllocationCache::AccViewLookup(void *AccPtr){assert(0); return 0;}

NAMESPACE_END(Grid);
#endif
