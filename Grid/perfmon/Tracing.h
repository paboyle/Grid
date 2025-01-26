#pragma once

NAMESPACE_BEGIN(Grid);

#ifdef GRID_TRACING_NVTX
#include <nvToolsExt.h>
class GridTracer {
public:
  GridTracer(const char* name) {
    nvtxRangePushA(name);
  }
  ~GridTracer() {
    nvtxRangePop();
  }
};
inline void tracePush(const char *name) { nvtxRangePushA(name); }
inline void tracePop(const char *name) { nvtxRangePop(); }
inline int  traceStart(const char *name) {  }
inline void traceStop(int ID) {  }
#endif

#ifdef GRID_TRACING_ROCTX
#include <roctracer/roctx.h>
class GridTracer {
 public:
  GridTracer(const char* name) {
    roctxRangePushA(name);
    std::cout << "roctxRangePush "<<name<<std::endl;
  }
  ~GridTracer() {
    roctxRangePop();
    std::cout << "roctxRangePop "<<std::endl;
  }
};
inline void tracePush(const char *name) { roctxRangePushA(name); }
inline void tracePop(const char *name) { roctxRangePop(); }
inline int  traceStart(const char *name) { return roctxRangeStart(name); }
inline void traceStop(int ID) { roctxRangeStop(ID); }
#endif

#ifdef GRID_TRACING_TIMER
class GridTracer {
 public:
  const char *name;
  double elapsed;
  GridTracer(const char* _name) {
    name = _name;
    elapsed=-usecond();
  }
  ~GridTracer() {
    elapsed+=usecond();
    std::cout << GridLogTracing << name << " took " <<elapsed<< " us" <<std::endl;
  }
};
inline void tracePush(const char *name) {  }
inline void tracePop(const char *name) {  }
inline int  traceStart(const char *name) { return 0; }
inline void traceStop(int ID) {  }
#endif

#ifdef GRID_TRACING_NONE
#define GRID_TRACE(name) 
inline void tracePush(const char *name) {  }
inline void tracePop(const char *name) {  }
inline int  traceStart(const char *name) { return 0;  }
inline void traceStop(int ID) {  }
#else
#define GRID_TRACE(name) GridTracer uniq_name_using_macros##__COUNTER__(name);
#endif
NAMESPACE_END(Grid);
