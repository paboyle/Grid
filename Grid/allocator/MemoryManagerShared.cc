#include <Grid/GridCore.h>
#ifdef GRID_UVM

#warning "Grid is assuming unified virtual memory address space"
NAMESPACE_BEGIN(Grid);
/////////////////////////////////////////////////////////////////////////////////
// View management is 1:1 address space mapping
/////////////////////////////////////////////////////////////////////////////////

void  MemoryManager::ViewClose(void* AccPtr,ViewMode mode){};
void *MemoryManager::ViewOpen(void* CpuPtr,size_t bytes,ViewMode mode,ViewAdvise hint){ return CpuPtr; };
void  MemoryManager::Print(void){};
void  MemoryManager::NotifyDeletion(void *ptr){};

NAMESPACE_END(Grid);
#endif
