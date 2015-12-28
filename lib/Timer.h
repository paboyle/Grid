#ifndef GRID_TIME_H
#define GRID_TIME_H

#include <sys/time.h>
#include <ctime>
#include <chrono>

namespace Grid {


  // Dress the output; use std::chrono

// C++11 time facilities better?
double usecond(void);

typedef  std::chrono::system_clock          GridClock;
typedef  std::chrono::time_point<GridClock> GridTimePoint;
typedef  std::chrono::milliseconds          GridTime;

inline std::ostream& operator<< (std::ostream & stream, const std::chrono::milliseconds & time)
{
  stream << time.count()<<" ms";
  return stream;
}
 
class GridStopWatch {
private:
  bool running;
  GridTimePoint start;
  GridTime accumulator;
public:
  GridStopWatch () { 
    Reset();
  }
  void     Start(void) { 
    assert(running == false);
    start = GridClock::now(); 
    running = true;
  }
  void     Stop(void)  { 
    assert(running == true);
    accumulator+= std::chrono::duration_cast<GridTime>(GridClock::now()-start); 
    running = false; 
  };
  void     Reset(void){
    running = false;
    start = GridClock::now();
    accumulator = std::chrono::duration_cast<GridTime>(start-start); 
  }
  GridTime Elapsed(void) {
    assert(running == false);
    return accumulator;
  }
};

}
#endif
