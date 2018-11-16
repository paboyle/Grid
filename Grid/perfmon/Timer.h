    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/Timer.h

    Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

    See the full license in the file "LICENSE" in the top level distribution directory
    *************************************************************************************/
    /*  END LEGAL */
#ifndef GRID_TIME_H
#define GRID_TIME_H

#include <sys/time.h>
#include <ctime>
#include <chrono>

namespace Grid {


  // Dress the output; use std::chrono

// C++11 time facilities better?
inline double usecond(void) {
  struct timeval tv;
#ifdef TIMERS_ON
  gettimeofday(&tv,NULL);
#endif
  return 1.0*tv.tv_usec + 1.0e6*tv.tv_sec;
}

typedef  std::chrono::system_clock          GridClock;
typedef  std::chrono::time_point<GridClock> GridTimePoint;

typedef  std::chrono::seconds               GridSecs;
typedef  std::chrono::milliseconds          GridMillisecs;
typedef  std::chrono::microseconds          GridUsecs;
typedef  std::chrono::microseconds          GridTime;

inline std::ostream& operator<< (std::ostream & stream, const GridSecs & time)
{
  stream << time.count()<<" s";
  return stream;
}
inline std::ostream& operator<< (std::ostream & stream, const GridMillisecs & now)
{
  GridSecs second(1);
  auto     secs       = now/second ; 
  auto     subseconds = now%second ; 
  stream << secs<<"."<<std::setw(3)<<std::setfill('0')<<subseconds.count()<<" s";
  return stream;
}
inline std::ostream& operator<< (std::ostream & stream, const GridUsecs & now)
{
  GridSecs second(1);
  auto     seconds    = now/second ; 
  auto     subseconds = now%second ; 
  stream << seconds<<"."<<std::setw(6)<<std::setfill('0')<<subseconds.count()<<" s";
  return stream;
}


class GridStopWatch {
private:
  bool running;
  GridTimePoint start;
  GridUsecs accumulator;
public:
  GridStopWatch () { 
    Reset();
  }
  void     Start(void) { 
    assert(running == false);
#ifdef TIMERS_ON
    start = GridClock::now(); 
#endif
    running = true;
  }
  void     Stop(void)  { 
    assert(running == true);
#ifdef TIMERS_ON
    accumulator+= std::chrono::duration_cast<GridUsecs>(GridClock::now()-start); 
#endif
    running = false; 
  };
  void     Reset(void){
    running = false;
#ifdef TIMERS_ON
    start = GridClock::now();
#endif
    accumulator = std::chrono::duration_cast<GridUsecs>(start-start); 
  }
  GridTime Elapsed(void) {
    assert(running == false);
    return std::chrono::duration_cast<GridTime>( accumulator );
  }
  uint64_t useconds(void){
    assert(running == false);
    return (uint64_t) accumulator.count();
  }
  bool isRunning(void){
    return running;
  }
};

}
#endif
