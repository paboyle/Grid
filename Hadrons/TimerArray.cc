/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/TimerArray.cc

Copyright (C) 2015-2018

Author: Antonin Portelli <antonin.portelli@me.com>

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
