#ifndef GRID_LOG_H
#define GRID_LOG_H
namespace Grid {

// Dress the output; use std::chrono for time stamping via the StopWatch class


class Logger {
protected:
    int active;
    std::string name, topName;
public:
    static GridStopWatch StopWatch;
    static std::ostream devnull;
    
    Logger(std::string topNm, int on, std::string nm)
    : active(on), name(nm), topName(topNm) {};
    
    void Active(int on) {active = on;};
    int  isActive(void) {return active;};
    
    friend std::ostream& operator<< (std::ostream& stream, const Logger& log){
        if ( log.active ) {
            StopWatch.Stop();
            GridTime now = StopWatch.Elapsed();
            StopWatch.Start();
            stream << std::setw(8) << std::left << log.topName << " : ";
            stream << std::setw(12) << std::left << log.name << " : ";
            stream << now << " : ";
            return stream;
        } else { 
            return devnull;
        }
    }
    
};
    
class GridLogger: public Logger {
public:
  GridLogger(int on, std::string nm): Logger("Grid", on, nm){};
};

void GridLogConfigure(std::vector<std::string> &logstreams);

extern GridLogger GridLogError;
extern GridLogger GridLogWarning;
extern GridLogger GridLogMessage;
extern GridLogger GridLogDebug  ;
extern GridLogger GridLogPerformance;
extern GridLogger GridLogIterative  ;
extern GridLogger GridLogIntegrator  ;

}
#endif
