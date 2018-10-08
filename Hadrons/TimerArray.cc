#include <Hadrons/TimerArray.hpp>

using namespace Grid;
using namespace QCD;
using namespace Hadrons;

void TimerArray::startTimer(const std::string &name)
{
    if (!name.empty())
    {
        timer_[name].Start();
    }
}

GridTime TimerArray::getTimer(const std::string &name)
{
    GridTime t;
    
    if (!name.empty())
    {
        try
        {
            bool running = timer_.at(name).isRunning();

            if (running) stopTimer(name);
            t = timer_.at(name).Elapsed();
            if (running) startTimer(name);
        }
        catch (std::out_of_range &)
        {
            t = GridTime::zero();
        }
    }
    else
    {
        t = GridTime::zero();
    }

    return t;
}

double TimerArray::getDTimer(const std::string &name)
{
    return static_cast<double>(getTimer(name).count());
}

void TimerArray::startCurrentTimer(const std::string &name)
{
    if (!name.empty())
    {
        stopCurrentTimer();
        startTimer(name);
        currentTimer_ = name;
    }
}

void TimerArray::stopTimer(const std::string &name)
{
    if (timer_.at(name).isRunning())
    {
        timer_.at(name).Stop();
    }
}

void TimerArray::stopCurrentTimer(void)
{
    if (!currentTimer_.empty())
    {
        stopTimer(currentTimer_);
        currentTimer_ = "";
    }
}

void TimerArray::stopAllTimers(void)
{
    for (auto &t: timer_)
    {
        stopTimer(t.first);
    }
    currentTimer_ = "";
}

void TimerArray::resetTimers(void)
{
    timer_.clear();
    currentTimer_ = "";
}

std::map<std::string, GridTime> TimerArray::getTimings(void)
{
    std::map<std::string, GridTime> timing;

    for (auto &t: timer_)
    {
        timing[t.first] = t.second.Elapsed();
    }

    return timing;
}
