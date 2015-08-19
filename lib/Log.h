#ifndef GRID_LOG_H
#define GRID_LOG_H
namespace Grid {

// Dress the output; use std::chrono for time stamping via the StopWatch class

std::ostream& operator<< (std::ostream& stream, const GridTime& time);

class GridLogger { 
  int active;
  std::string name;
public:

  static GridStopWatch StopWatch;
  static std::ostream devnull;
  
  GridLogger(int on, std::string nm): active(on), name(nm) { 
  };
  
  void Active(int on) {active = on;};

  friend std::ostream& operator<< (std::ostream& stream, const GridLogger& log){
    if ( log.active ) {
      StopWatch.Stop();
      GridTime now = StopWatch.Elapsed();
      StopWatch.Start();
      stream << "Grid : "<<log.name << " : " << now << " : ";
      return stream;
    } else { 
      return devnull;
    }
  }

};

void GridLogConfigure(std::vector<std::string> &logstreams);

extern GridLogger GridLogError;
extern GridLogger GridLogWarning;
extern GridLogger GridLogMessage;
extern GridLogger GridLogDebug  ;
extern GridLogger GridLogPerformance;
extern GridLogger GridLogIterative  ;

}
#endif
