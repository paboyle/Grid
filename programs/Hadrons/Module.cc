/*
 * Module.cc, part of Grid
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

#include <Hadrons/Module.hpp>

using namespace std;
using namespace Grid;
using namespace Hadrons;

/******************************************************************************
 *                           Module implementation                            *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
Module::Module(const string &name)
: name_(name)
{}

// access //////////////////////////////////////////////////////////////////////
string Module::getName(void) const
{
    return name_;
}
