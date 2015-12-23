/*
 * CMeson.cc, part of Grid
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

#include <Hadrons/CMeson.hpp>

using namespace Grid;
using namespace Hadrons;

/******************************************************************************
 *                          CMeson implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
CMeson::CMeson(const std::string &name)
: Module(name)
{}

// parse parameters ////////////////////////////////////////////////////////////
void CMeson::parseParameters(XmlReader &reader, const std::string &name)
{
    read(reader, name, par_);
}

// dependency relation /////////////////////////////////////////////////////////
std::vector<std::string> CMeson::getInput(void)
{
    std::vector<std::string> input = {par_.q1, par_.q2};
    
    return input;
}

std::vector<std::string> CMeson::getOutput(void)
{
    std::vector<std::string> output = {getName()};
    
    return output;
}

// memory footprint ////////////////////////////////////////////////////////////
double CMeson::nCreatedProp(void)
{
    return 0.;
}

// execution ///////////////////////////////////////////////////////////////////
void CMeson::operator()(Environment &env)
{
    LOG(Message) << "computing meson contraction '" << getName() << "'"
                 << std::endl;
}
