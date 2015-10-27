/*
 * Application.hpp, part of Grid
 *
 * Copyright (C) 2015 Antonin Portelli
 *
 * Grid is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Grid is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Grid.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef Hadrons_Application_hpp_
#define Hadrons_Application_hpp_

#include <Hadrons/Global.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         Main program manager                               *
 ******************************************************************************/
class Application
{
public:
    // constructor
    Application(int argc, char *argv[]);
    // destructor
    virtual ~Application(void);
    // execute
    void run(void);
private:
    void parseParameters(void);
private:
    std::string parameterFileName_;
};

END_HADRONS_NAMESPACE

#endif // Hadrons_Application_hpp_
