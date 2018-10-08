#ifndef Hadrons_TimerArray_hpp_
#define Hadrons_TimerArray_hpp_

#include <Hadrons/Global.hpp>

BEGIN_HADRONS_NAMESPACE

class TimerArray
{
public:
    TimerArray(void) = default;
    virtual ~TimerArray(void) = default;
    void                            startTimer(const std::string &name);
    GridTime                        getTimer(const std::string &name);
    double                          getDTimer(const std::string &name);
    void                            startCurrentTimer(const std::string &name);
    void                            stopTimer(const std::string &name);
    void                            stopCurrentTimer(void);
    void                            stopAllTimers(void);
    void                            resetTimers(void);
    std::map<std::string, GridTime> getTimings(void);
private:
    std::string                          currentTimer_;
    std::map<std::string, GridStopWatch> timer_; 
};

END_HADRONS_NAMESPACE

#endif // Hadrons_TimerArray_hpp_
