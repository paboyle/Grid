#include <Grid/GridCore.h>
#ifdef GRID_UVM

#warning "Grid is assuming unified virtual memory address space"
NAMESPACE_BEGIN(Grid);
/////////////////////////////////////////////////////////////////////////////////
// View management is 1:1 address space mapping
/////////////////////////////////////////////////////////////////////////////////
uint64_t  MemoryManager::DeviceBytes;
uint64_t  MemoryManager::DeviceLRUBytes;
uint64_t  MemoryManager::DeviceMaxBytes = 1024*1024*128;
uint64_t  MemoryManager::HostToDeviceBytes;
uint64_t  MemoryManager::DeviceToHostBytes;
uint64_t  MemoryManager::HostToDeviceXfer;
uint64_t  MemoryManager::DeviceToHostXfer;

void  MemoryManager::ViewClose(void* AccPtr,ViewMode mode){};
void *MemoryManager::ViewOpen(void* CpuPtr,size_t bytes,ViewMode mode,ViewAdvise hint){ return CpuPtr; };
int   MemoryManager::isOpen   (void* CpuPtr) { return 0;}
void  MemoryManager::Print(void){};
void  MemoryManager::NotifyDeletion(void *ptr){};

NAMESPACE_END(Grid);
#endif
