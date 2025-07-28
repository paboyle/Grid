#include <Grid/GridCore.h>

NAMESPACE_BEGIN(Grid);

#define DEVICE_MEMORY_ALLOCATOR_PAGE_SIZE (64*1024)


struct DeviceMemoryAllocator {

  bool initialized;
  char* base;
  size_t size;
  size_t offset;

  DeviceMemoryAllocator() {
    initialized = false;
    base = 0;
    size = 0;
    offset = 0;
  }

  ~DeviceMemoryAllocator() {
    if (initialized) {
      acceleratorFreeDevice(base);
      initialized = false;
    }
  }

  std::vector<size_t> pages;
  std::map<size_t, std::vector<size_t> > size_map;

  void Init(size_t _size) {
    assert(!initialized);
    
    size_t n_pages = (_size + DEVICE_MEMORY_ALLOCATOR_PAGE_SIZE - 1) / DEVICE_MEMORY_ALLOCATOR_PAGE_SIZE;
    size = n_pages * DEVICE_MEMORY_ALLOCATOR_PAGE_SIZE;
    std::cout << GridLogMessage << "Init with " <<
      size << " bytes" << std::endl;

    base = (char*)acceleratorAllocDeviceInternal(size);
    assert(base);
    
    std::cout << GridLogMessage << "Initialize memory to zero" << std::endl;

    {
      uint64_t* ba = (uint64_t*)base;
      size_t n = size / sizeof(uint64_t);
      accelerator_for(i, n, 1, {
	  ba[i] = (uint64_t)-1;
	});
    }  
    
    std::cout << GridLogMessage << "Done" << std::endl;
    
    offset = 0;
    
    pages.resize(n_pages, 0);
    
    std::cout << GridLogMessage << "Pages initialized" << std::endl;
  
    initialized = true;
  }
};

static DeviceMemoryAllocator dma;

void *acceleratorAllocDevice(size_t bytes) {
  if (!dma.initialized)
    dma.Init(MemoryManager::DeviceMaxBytes);
  
  if (!bytes)
    bytes++;
  
  size_t n_pages = (bytes + DEVICE_MEMORY_ALLOCATOR_PAGE_SIZE - 1) / DEVICE_MEMORY_ALLOCATOR_PAGE_SIZE;
  
  // first check if
  auto & sm = dma.size_map[n_pages];
  if (sm.size()) {
    size_t index = sm.back();
    sm.pop_back();
    assert(dma.pages[index] == 0);
    dma.pages[index] = n_pages;
    
    return dma.base + index * DEVICE_MEMORY_ALLOCATOR_PAGE_SIZE;
  }

  size_t end = (dma.offset + n_pages) * DEVICE_MEMORY_ALLOCATOR_PAGE_SIZE;
  assert(end <= dma.size);
  dma.pages[dma.offset] = n_pages;

  void* ptr = dma.base + dma.offset * DEVICE_MEMORY_ALLOCATOR_PAGE_SIZE;
  dma.offset += n_pages;

  std::cout << GridLogMessage << (dma.size - end) / DEVICE_MEMORY_ALLOCATOR_PAGE_SIZE << " pages left to allocate ("
	    << (dma.size - end) * 100 / dma.size << "% free)" << std::endl;

  return ptr;

}

void acceleratorFreeDevice(void *ptr) {
  if (!dma.initialized)
    return;
  
  size_t index = ((size_t)((char*)ptr - dma.base)) / DEVICE_MEMORY_ALLOCATOR_PAGE_SIZE;
  size_t n_pages = dma.pages[index];
  //std::cout << GridLogMessage << "Freeing ptr " << ptr << " has " << n_pages << " pages" << std::endl;
  dma.pages[index] = 0;
  auto & sm = dma.size_map[n_pages];
  sm.push_back(index);
  
}

NAMESPACE_END(Grid);
