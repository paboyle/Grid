    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/Log.h

    Copyright (C) 2015

Author: Antonin Portelli <antonin.portelli@me.com>
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
#ifndef GRID_LOG_H
#define GRID_LOG_H

#ifdef HAVE_EXECINFO_H
#include <execinfo.h>
#endif

namespace Grid {

// Dress the output; use std::chrono for time stamping via the StopWatch class
int Rank(void); // used for early stage debug before library init


class Logger {
protected:
    int active;
    std::string name, topName, COLOUR;
public:
    static GridStopWatch StopWatch;
    static std::ostream devnull;

    static std::string BLACK;
    static std::string RED  ;
    static std::string GREEN;
    static std::string YELLOW;
    static std::string BLUE  ;
    static std::string PURPLE;
    static std::string CYAN  ;
    static std::string WHITE ;
    static std::string NORMAL;
    
 Logger(std::string topNm, int on, std::string nm,std::string col)
   : active(on), name(nm), topName(topNm), COLOUR(col) {};
    
    void Active(int on) {active = on;};
    int  isActive(void) {return active;};
    
    friend std::ostream& operator<< (std::ostream& stream, const Logger& log){
        if ( log.active ) {
            StopWatch.Stop();
            GridTime now = StopWatch.Elapsed();
            StopWatch.Start();
            stream << BLACK<< log.topName << BLACK<< " : ";
            stream << log.COLOUR <<std::setw(10) << std::left << log.name << BLACK << " : ";
	    stream << YELLOW<< now <<BLACK << " : " ;
	    stream << log.COLOUR;
            return stream;
        } else { 
            return devnull;
        }
    }
    
};
    
class GridLogger: public Logger {
public:
 GridLogger(int on, std::string nm, std::string col = Logger::BLACK): Logger("Grid", on, nm, col){};
};

void GridLogConfigure(std::vector<std::string> &logstreams);

extern GridLogger GridLogError;
extern GridLogger GridLogWarning;
extern GridLogger GridLogMessage;
extern GridLogger GridLogDebug  ;
extern GridLogger GridLogPerformance;
extern GridLogger GridLogIterative  ;
extern GridLogger GridLogIntegrator  ;


#define _NBACKTRACE (256)
extern void * Grid_backtrace_buffer[_NBACKTRACE];

#define BACKTRACEFILE() {\
    char string[20];					\
    std::sprintf(string,"backtrace.%d",Rank());				\
    std::FILE * fp = std::fopen(string,"w");				\
    BACKTRACEFP(fp)\
    std::fclose(fp);	    \
}


#ifdef HAVE_EXECINFO_H
#define BACKTRACEFP(fp) { \
  int symbols    = backtrace        (Grid_backtrace_buffer,_NBACKTRACE);\
  char **strings = backtrace_symbols(Grid_backtrace_buffer,symbols);\
  for (int i = 0; i < symbols; i++){\
    std::fprintf (fp,"BackTrace Strings: %d %s\n",i, strings[i]); std::fflush(fp); \
  }\
}
#else 
#define BACKTRACEFP(fp) { \
    std::fprintf (fp,"BT %d %lx\n",0, __builtin_return_address(0)); std::fflush(fp); \
    std::fprintf (fp,"BT %d %lx\n",1, __builtin_return_address(1)); std::fflush(fp); \
    std::fprintf (fp,"BT %d %lx\n",2, __builtin_return_address(2)); std::fflush(fp); \
    std::fprintf (fp,"BT %d %lx\n",3, __builtin_return_address(3)); std::fflush(fp); \
}
#endif

#define BACKTRACE() BACKTRACEFP(stdout) 

}
#endif
