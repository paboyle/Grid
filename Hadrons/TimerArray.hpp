/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/TimerArray.hpp

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
