/*******************************************************************************
Grid physics library, www.github.com/paboyle/Grid 

Source file: programs/Hadrons/Environment.cc

Copyright (C) 2015

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

See the full license in the file "LICENSE" in the top level distribution 
directory.
*******************************************************************************/

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
    std::vector<int> seed4d({1,2,3,4});
    
    grid4d_.reset(SpaceTimeGrid::makeFourDimGrid(
        GridDefaultLatt(), GridDefaultSimd(Nd, vComplex::Nsimd()),
        GridDefaultMpi()));
    gridRb4d_.reset(SpaceTimeGrid::makeFourDimRedBlackGrid(grid4d_.get()));
    rng4d_.reset(new GridParallelRNG(grid4d_.get()));
    rng4d_->SeedFixedIntegers(seed4d);
    gauge_.reset(new LatticeGaugeField(grid4d_.get()));
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

// grids ///////////////////////////////////////////////////////////////////////
GridCartesian * Environment::get4dGrid(void)
{
    return grid4d_.get();
}

GridRedBlackCartesian * Environment::getRb4dGrid(void)
{
    return gridRb4d_.get();
}

GridCartesian * Environment::get5dGrid(const unsigned int Ls)
{
    try
    {
        return grid5d_.at(Ls).get();
    }
    catch(std::out_of_range &)
    {
        HADRON_ERROR("no 5D grid with Ls= " << Ls);
    }
}

GridRedBlackCartesian * Environment::getRb5dGrid(const unsigned int Ls)
{
    try
    {
        return gridRb5d_.at(Ls).get();
    }
    catch(std::out_of_range &)
    {
        HADRON_ERROR("no red-black 5D grid with Ls= " << Ls);
    }
}

// quark propagators ///////////////////////////////////////////////////////////
void Environment::addProp(const std::string name, const unsigned int Ls)
{
    GridCartesian *p4 = grid4d_.get();

    if (propExists(name))
    {
        HADRON_ERROR("propagator '" + name + "' already exists");
    }
    if (Ls > 1)
    {
        GridCartesian *p;
        
        try
        {
            p = grid5d_.at(Ls).get();
        }
        catch(std::out_of_range &)
        {
            grid5d_[Ls].reset(SpaceTimeGrid::makeFiveDimGrid(Ls, p4));
            gridRb5d_[Ls].reset(SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls, p4));
            p = grid5d_[Ls].get();
        }
        if (!isDryRun())
        {
            prop_[name].reset(new LatticePropagator(p));
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
            prop_[name].reset(new LatticePropagator(p4));
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

// gauge configuration /////////////////////////////////////////////////////////
LatticeGaugeField * Environment::getGauge(void)
{
    return gauge_.get();
}

void Environment::loadUnitGauge(void)
{
    SU3::ColdConfiguration(*rng4d_, *gauge_);
}

void Environment::loadRandomGauge(void)
{
    SU3::HotConfiguration(*rng4d_, *gauge_);
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
