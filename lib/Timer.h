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
