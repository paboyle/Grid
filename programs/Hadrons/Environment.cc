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
    grid4d_.reset(SpaceTimeGrid::makeFourDimGrid(
        GridDefaultLatt(), GridDefaultSimd(Nd, vComplex::Nsimd()),
        GridDefaultMpi()));
    gridRb4d_.reset(SpaceTimeGrid::makeFourDimRedBlackGrid(grid4d_.get()));
    rng4d_.reset(new GridParallelRNG(grid4d_.get()));
    gauge_.reset(new LatticeGaugeField(grid4d_.get()));
    loadUnitGauge();
}

// dry run /////////////////////////////////////////////////////////////////////
void Environment::dryRun(const bool isDry)
{
    dryRun_ = isDry;
}

bool Environment::isDryRun(void) const
{
    return dryRun_;
}

// grids ///////////////////////////////////////////////////////////////////////
GridCartesian * Environment::getGrid(const unsigned int Ls) const
{
    try
    {
        if (Ls == 1)
        {
            return grid4d_.get();
        }
        else
        {
            return grid5d_.at(Ls).get();
        }
    }
    catch(std::out_of_range &)
    {
        HADRON_ERROR("no 5D grid with Ls= " << Ls);
    }
}

GridRedBlackCartesian * Environment::getRbGrid(const unsigned int Ls) const
{
    try
    {
        if (Ls == 1)
        {
            return gridRb4d_.get();
        }
        else
        {
            return gridRb5d_.at(Ls).get();
        }
    }
    catch(std::out_of_range &)
    {
        HADRON_ERROR("no red-black 5D grid with Ls= " << Ls);
    }
}

// fermion actions /////////////////////////////////////////////////////////////
void Environment::addFermionAction(FActionPt action)
{
    fAction_[action->getName()] = std::move(action);
}

FermionAction * Environment::getFermionAction(const std::string name) const
{
    try
    {
        return fAction_.at(name).get();
    }
    catch(std::out_of_range &)
    {
        try
        {
            return fAction_.at(solverAction_.at(name)).get();
        }
        catch (std::out_of_range &)
        {
            HADRON_ERROR("no action with name '" << name << "'");
        }
    }
}

// solvers /////////////////////////////////////////////////////////////////////
void Environment::addSolver(const std::string name, Solver s,
                            const std::string actionName)
{
    solver_[name]       = s;
    solverAction_[name] = actionName;
}

void Environment::callSolver(const std::string name, LatticeFermion &sol,
                             const LatticeFermion &source) const
{
    try
    {
        solver_.at(name)(sol, source);
    }
    catch(std::out_of_range &)
    {
        HADRON_ERROR("no solver with name '" << name << "'");
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

bool Environment::isProp5d(const std::string name) const
{
    if (propExists(name))
    {
        return (getProp(name)->_grid->GlobalDimensions().size() == Nd + 1);
    }
    else
    {
        HADRON_ERROR("propagator '" + name + "' unknown");
        
        return false;
    }
}

unsigned int Environment::getPropLs(const std::string name) const
{
    if (propExists(name))
    {
        if (isProp5d(name))
        {
            return getProp(name)->_grid->GlobalDimensions()[0];
        }
        else
        {
            return 1;
        }
    }
    else
    {
        HADRON_ERROR("propagator '" + name + "' unknown");
        
        return 0;
    }
}

LatticePropagator * Environment::getProp(const std::string name) const
{
    if (propExists(name))
    {
        return prop_.at(name).get();
    }
    else
    {
        HADRON_ERROR("propagator '" + name + "' unknown");
        
        return nullptr;
    }
}

bool Environment::propExists(const std::string name) const
{
    return (prop_.find(name) != prop_.end());
}

unsigned int Environment::nProp(void) const
{
    unsigned int size = 0;
    
    for (auto &s: propSize_)
    {
        size += s.second;
    }
    
    return size;
}

// gauge configuration /////////////////////////////////////////////////////////
LatticeGaugeField * Environment::getGauge(void) const
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

// random number generator /////////////////////////////////////////////////////
void Environment::setSeed(const std::vector<int> &seed)
{
    rng4d_->SeedFixedIntegers(seed);
}

GridParallelRNG * Environment::get4dRng(void) const
{
    return rng4d_.get();
}

// general free ////////////////////////////////////////////////////////////////
void Environment::free(const std::string name)
{
    if (propExists(name))
    {
        LOG(Message) << "freeing '" << name << "'" << std::endl;
        freeProp(name);
    }
}

void Environment::freeAll(void)
{
    prop_.clear();
    propSize_.clear();
}
