#include <Grid/GridCore.h>

NAMESPACE_BEGIN(Grid);
uint32_t accelerator_threads;
uint32_t acceleratorThreads(void)       {return accelerator_threads;};
void     acceleratorThreads(uint32_t t) {accelerator_threads = t;};
#ifdef GRID_SYCL
cl::sycl::queue *theGridAccelerator;
#endif
NAMESPACE_END(Grid);
