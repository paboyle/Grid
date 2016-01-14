/*
 * Environment.cc, part of Grid
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

#include <Hadrons/Environment.hpp>

using namespace Grid;
using namespace QCD;
using namespace Hadrons;

/******************************************************************************
 *                       Environment implementation                           *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
Environment::Environment(void)
{
    grid4d_.reset(SpaceTimeGrid::makeFourDimGrid(
        GridDefaultLatt(), GridDefaultSimd(Nd, vComplex::Nsimd()),
        GridDefaultMpi()));
}

// dry run /////////////////////////////////////////////////////////////////////
void Environment::dryRun(const bool isDry)
{
    dryRun_ = isDry;
}

bool Environment::isDryRun(void)
{
    return dryRun_;
}

// quark propagators ///////////////////////////////////////////////////////////
void Environment::addProp(const std::string name, const unsigned int Ls)
{
    if (propExists(name))
    {
        HADRON_ERROR("propagator '" + name + "' already exists");
    }
    if (Ls > 1)
    {
        GridCartesian *pt;
        
        try
        {
            pt = grid5d_.at(Ls).get();
        }
        catch(std::out_of_range &)
        {
            grid5d_[Ls].reset(SpaceTimeGrid::makeFiveDimGrid(Ls,
                                                             grid4d_.get()));
            pt = grid5d_[Ls].get();
        }
        if (!isDryRun())
        {
            prop_[name].reset(new LatticePropagator(pt));
        }
        else
        {
            prop_[name].reset(nullptr);
        }
        propSize_[name] = Ls;
    }
    else
    {
        if (!isDryRun())
        {
            prop_[name].reset(new LatticePropagator(grid4d_.get()));
        }
        else
        {
            prop_[name].reset(nullptr);
        }
        propSize_[name] = 1;
    }
}

void Environment::freeProp(const std::string name)
{
    if (propExists(name))
    {
        prop_.erase(name);
        propSize_.erase(name);
    }
    else
    {
        HADRON_ERROR("trying to free unknown propagator '" + name + "'");
    }
}

LatticePropagator * Environment::getProp(const std::string name)
{
    if (propExists(name))
    {
        return prop_[name].get();
    }
    else
    {
        HADRON_ERROR("propagator '" + name + "' unknown");
        
        return nullptr;
    }
}

bool Environment::propExists(const std::string name)
{
    auto it = prop_.find(name);
    
    if (it == prop_.end())
    {
        return false;
    }
    else
    {
        return true;
    }
}

unsigned int Environment::nProp(void)
{
    unsigned int size = 0;
    
    for (auto &s: propSize_)
    {
        size += s.second;
    }
    
    return size;
}

// general free ////////////////////////////////////////////////////////////////
void Environment::free(const std::string name)
{
    if (propExists(name))
    {
        freeProp(name);
    }
}

void Environment::freeAll(void)
{
    prop_.clear();
    propSize_.clear();
}
