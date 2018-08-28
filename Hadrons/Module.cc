/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: extras/Hadrons/Module.cc

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

#include <Hadrons/Module.hpp>

using namespace Grid;
using namespace QCD;
using namespace Hadrons;

/******************************************************************************
 *                       ModuleBase implementation                            *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
ModuleBase::ModuleBase(const std::string name)
: name_(name)
{}

// access //////////////////////////////////////////////////////////////////////
std::string ModuleBase::getName(void) const
{
    return name_;
}

// get factory registration name if available
std::string ModuleBase::getRegisteredName(void)
{
    HADRONS_ERROR(Definition, "module '" + getName() + "' has no registered type"
                 + " in the factory");
}

// execution ///////////////////////////////////////////////////////////////////
void ModuleBase::operator()(void)
{
    resetTimers();
    startTimer("_total");
    startTimer("_setup");
    setup();
    stopTimer("_setup");
    startTimer("_execute");
    execute();
    stopAllTimers();
}

// timers //////////////////////////////////////////////////////////////////////
void ModuleBase::startTimer(const std::string &name)
{
    if (!name.empty())
    {
        timer_[name].Start();
    }
}

GridTime ModuleBase::getTimer(const std::string &name)
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

double ModuleBase::getDTimer(const std::string &name)
{
    return static_cast<double>(getTimer(name).count());
}

void ModuleBase::startCurrentTimer(const std::string &name)
{
    if (!name.empty())
    {
        stopCurrentTimer();
        startTimer(name);
        currentTimer_ = name;
    }
}

void ModuleBase::stopTimer(const std::string &name)
{
    if (timer_.at(name).isRunning())
    {
        timer_.at(name).Stop();
    }
}

void ModuleBase::stopCurrentTimer(void)
{
    if (!currentTimer_.empty())
    {
        stopTimer(currentTimer_);
        currentTimer_ = "";
    }
}

void ModuleBase::stopAllTimers(void)
{
    for (auto &t: timer_)
    {
        stopTimer(t.first);
    }
    currentTimer_ = "";
}

void ModuleBase::resetTimers(void)
{
    timer_.clear();
    currentTimer_ = "";
}

std::map<std::string, GridTime> ModuleBase::getTimings(void)
{
    std::map<std::string, GridTime> timing;

    for (auto &t: timer_)
    {
        timing[t.first] = t.second.Elapsed();
    }

    return timing;
}

std::string ModuleBase::makeSeedString(void)
{
    std::string seed;

    if (!vm().getRunId().empty())
    {
        seed += vm().getRunId() + "-";
    }
    seed += getName() + "-" + std::to_string(vm().getTrajectory());

    return seed;
}

GridParallelRNG & ModuleBase::rng4d(void)
{
    auto &r = *env().get4dRng();

    if (makeSeedString() != seed_)
    {
        seed_ = makeSeedString();
        LOG(Message) << "Seeding 4D RNG " << &r << " with string '" 
                     << seed_ << "'" << std::endl;
        r.SeedUniqueString(seed_);
    }

    return r;
}
