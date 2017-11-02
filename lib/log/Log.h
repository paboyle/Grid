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

#include <map>

#ifndef GRID_LOG_H
#define GRID_LOG_H

#ifdef HAVE_EXECINFO_H
#include <execinfo.h>
#endif

namespace Grid {

//////////////////////////////////////////////////////////////////////////////////////////////////
// Dress the output; use std::chrono for time stamping via the StopWatch class
//////////////////////////////////////////////////////////////////////////////////////////////////


class Colours{
protected:
  bool is_active;
public:
  std::map<std::string, std::string> colour;

  Colours(bool activate=false){
    Active(activate);
  };

  void Active(bool activate){
    is_active=activate;
    if (is_active){
     colour["BLACK"]  ="\033[30m";
     colour["RED"]    ="\033[31m";
     colour["GREEN"]  ="\033[32m";
     colour["YELLOW"] ="\033[33m";
     colour["BLUE"]   ="\033[34m";
     colour["PURPLE"] ="\033[35m";
     colour["CYAN"]   ="\033[36m";
     colour["WHITE"]  ="\033[37m";
     colour["NORMAL"] ="\033[0;39m";
    } else {
      colour["BLACK"] ="";
      colour["RED"]   ="";
      colour["GREEN"] ="";
      colour["YELLOW"]="";
      colour["BLUE"]  ="";
      colour["PURPLE"]="";
      colour["CYAN"]  ="";
      colour["WHITE"] ="";
      colour["NORMAL"]="";
    }
  };
};


class Logger {
protected:
  Colours &Painter;
  int active;
  int timing_mode;
  static int timestamp;
  std::string name, topName;
  std::string COLOUR;

public:
  static GridStopWatch GlobalStopWatch;
  GridStopWatch         LocalStopWatch;
  GridStopWatch *StopWatch;
  static std::ostream devnull;

  std::string background() {return Painter.colour["NORMAL"];}
  std::string evidence() {return Painter.colour["YELLOW"];}
  std::string colour() {return Painter.colour[COLOUR];}

  Logger(std::string topNm, int on, std::string nm, Colours& col_class, std::string col)  : active(on),
    name(nm),
    topName(topNm),
    Painter(col_class),
    timing_mode(0),
    COLOUR(col) 
    {
      StopWatch = & GlobalStopWatch;
    };
  
  void Active(int on) {active = on;};
  int  isActive(void) {return active;};
  static void Timestamp(int on) {timestamp = on;};
  void Reset(void) { 
    StopWatch->Reset(); 
    StopWatch->Start(); 
  }
  void TimingMode(int on) { 
    timing_mode = on; 
    if(on) { 
      StopWatch = &LocalStopWatch;
      Reset(); 
    }
  }

  friend std::ostream& operator<< (std::ostream& stream, Logger& log){

    if ( log.active ) {
      stream << log.background()<<  std::left << log.topName << log.background()<< " : ";
      stream << log.colour() <<  std::left << log.name << log.background() << " : ";
      if ( log.timestamp ) {
	log.StopWatch->Stop();
	GridTime now = log.StopWatch->Elapsed();
	if ( log.timing_mode==1 ) log.StopWatch->Reset();
	log.StopWatch->Start();
	stream << log.evidence()<< std::setw(6)<<now << log.background() << " : " ;
      }
      stream << log.colour();
      return stream;
    } else { 
      return devnull;
    }
  }

};

class GridLogger: public Logger {
public:
  GridLogger(int on, std::string nm, Colours&col_class, std::string col_key = "NORMAL"):
  Logger("Grid", on, nm, col_class, col_key){};
};

void GridLogConfigure(std::vector<std::string> &logstreams);

extern GridLogger GridLogIRL;
extern GridLogger GridLogSolver;
extern GridLogger GridLogError;
extern GridLogger GridLogWarning;
extern GridLogger GridLogMessage;
extern GridLogger GridLogDebug  ;
extern GridLogger GridLogPerformance;
extern GridLogger GridLogIterative  ;
extern GridLogger GridLogIntegrator  ;
extern Colours    GridLogColours;

 std::string demangle(const char* name) ;

#define _NBACKTRACE (256)
extern void * Grid_backtrace_buffer[_NBACKTRACE];

#define BACKTRACEFILE() {\
char string[20];					\
std::sprintf(string,"backtrace.%d",CartesianCommunicator::RankWorld()); \
std::FILE * fp = std::fopen(string,"w");				\
BACKTRACEFP(fp)\
std::fclose(fp);	    \
}


#ifdef HAVE_EXECINFO_H
#define BACKTRACEFP(fp) { \
int symbols    = backtrace        (Grid_backtrace_buffer,_NBACKTRACE);\
char **strings = backtrace_symbols(Grid_backtrace_buffer,symbols);\
for (int i = 0; i < symbols; i++){\
  std::fprintf (fp,"BackTrace Strings: %d %s\n",i, demangle(strings[i]).c_str()); std::fflush(fp); \
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
