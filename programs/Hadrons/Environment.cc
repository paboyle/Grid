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

// trajectory number ///////////////////////////////////////////////////////////
void Environment::setTrajectory(const unsigned int traj)
{
    traj_ = traj;
}

unsigned int Environment::getTrajectory(void) const
{
    return traj_;
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
void Environment::addFermionMatrix(const std::string name, FMat *fMat)
{
    fMat_[name].reset(fMat);
}

Environment::FMat * Environment::getFermionMatrix(const std::string name) const
{
    try
    {
        return fMat_.at(name).get();
    }
    catch(std::out_of_range &)
    {
        try
        {
            return fMat_.at(solverAction_.at(name)).get();
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
void Environment::createProp(const std::string name, const unsigned int Ls)
{
    GridCartesian *g4 = getGrid();

    if (propExists(name))
    {
        HADRON_ERROR("propagator '" + name + "' already exists");
    }
    if (Ls > 1)
    {
        GridCartesian *g;
        
        try
        {
            g = grid5d_.at(Ls).get();
        }
        catch(std::out_of_range &)
        {
            grid5d_[Ls].reset(SpaceTimeGrid::makeFiveDimGrid(Ls, g4));
            gridRb5d_[Ls].reset(SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls, g4));
            g = grid5d_[Ls].get();
        }
        if (!isDryRun())
        {
            prop_[name].reset(new LatticePropagator(g));
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
            prop_[name].reset(new LatticePropagator(g4));
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
void Environment::createGauge(const std::string name)
{
    if (gaugeExists(name))
    {
        HADRON_ERROR("gauge field '" + name + "' already exists");
    }
    gauge_[name].reset(new LatticeGaugeField(getGrid()));
}

void Environment::freeGauge(const std::string name)
{
    if (gaugeExists(name))
    {
        gauge_.erase(name);
    }
    else
    {
        HADRON_ERROR("trying to free unknown gauge field '" + name + "'");
    }
}

LatticeGaugeField * Environment::getGauge(const std::string name) const
{
    if (gaugeExists(name))
    {
        return gauge_.at(name).get();
    }
    else
    {
        HADRON_ERROR("gauge field '" + name + "' unknown");
        
        return nullptr;
    }
}

bool Environment::gaugeExists(const std::string name) const
{
    return (gauge_.find(name) != gauge_.end());
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
        LOG(Message) << "freeing propagator '" << name << "'" << std::endl;
        freeProp(name);
    }
    else if (gaugeExists(name))
    {
        LOG(Message) << "freeing gauge field '" << name << "'" << std::endl;
        freeGauge(name);
    }
}

void Environment::freeAll(void)
{
    prop_.clear();
    propSize_.clear();
    gauge_.clear();
}
