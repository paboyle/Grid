#include <Grid/GridCore.h>
#ifdef GRID_UVM

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
uint64_t  MemoryManager::DeviceEvictions;
uint64_t  MemoryManager::DeviceDestroy;

void  MemoryManager::Audit(std::string s){};
void  MemoryManager::ViewClose(void* AccPtr,ViewMode mode){};
void *MemoryManager::ViewOpen(void* CpuPtr,size_t bytes,ViewMode mode,ViewAdvise hint){ return CpuPtr; };
int   MemoryManager::isOpen   (void* CpuPtr) { return 0;}
void  MemoryManager::PrintState(void* CpuPtr)
{
std::cout << GridLogMessage << "Host<->Device memory movement not currently managed by Grid." << std::endl;
};
void  MemoryManager::Print(void){};
void  MemoryManager::PrintAll(void){};
void  MemoryManager::NotifyDeletion(void *ptr){};

NAMESPACE_END(Grid);
#endif
