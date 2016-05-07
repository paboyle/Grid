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
void Environment::createGrid(const unsigned int Ls)
{
    auto g = getGrid();
    
    grid5d_[Ls].reset(SpaceTimeGrid::makeFiveDimGrid(Ls, g));
    gridRb5d_[Ls].reset(SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls, g));
}

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
        HADRON_ERROR("no grid with Ls= " << Ls);
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
    if (hasObject(name))
    {
        fMat_[name].reset(fMat);
    }
    else
    {
        HADRON_ERROR("no object named '" << name << "'");
    }
}

Environment::FMat * Environment::getFermionMatrix(const std::string name) const
{
    if (hasFermionMatrix(name))
    {
        return fMat_.at(name).get();
    }
    else
    {
        if (hasSolver(name))
        {
            return fMat_.at(solverAction_.at(name)).get();
        }
        else
        {
            HADRON_ERROR("no action/solver with name '" << name << "'");
        }
    }
}

void Environment::freeFermionMatrix(const std::string name)
{
    if (hasFermionMatrix(name))
    {
        LOG(Message) << "freeing fermion matrix '" << name << "'" << std::endl;
        fMat_.erase(name);
        object_.erase(name);
    }
    else
    {
        HADRON_ERROR("trying to free unknown fermion matrix '" + name + "'");
    }
}

bool Environment::hasFermionMatrix(const std::string name) const
{
    return (hasObject(name) and (fMat_.find(name) != fMat_.end()));
}

// solvers /////////////////////////////////////////////////////////////////////
void Environment::addSolver(const std::string name, Solver s)
{
    if (hasObject(name))
    {
        solver_[name] = s;
    }
    else
    {
        HADRON_ERROR("no object named '" << name << "'");
    }
}

bool Environment::hasSolver(const std::string name) const
{
    return (hasObject(name) and (solver_.find(name) != solver_.end()));
}

void Environment::setSolverAction(const std::string name,
                                  const std::string actionName)
{
    if (hasObject(name))
    {
        solverAction_[name] = actionName;
    }
    else
    {
        HADRON_ERROR("no object named '" << name << "'");
    }
}

std::string Environment::getSolverAction(const std::string name) const
{
    if (hasObject(name))
    {
        try
        {
            return solverAction_.at(name);
        }
        catch (std::out_of_range &)
        {
            HADRON_ERROR("no action registered for solver '" << name << "'")
        }
    }
    else
    {
        HADRON_ERROR("no object with name '" << name << "'");
    }
}

void Environment::callSolver(const std::string name, LatticeFermion &sol,
                             const LatticeFermion &source) const
{
    if (hasSolver(name))
    {
        solver_.at(name)(sol, source);
    }
    else
    {
        HADRON_ERROR("no solver with name '" << name << "'");
    }
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

// lattice store ///////////////////////////////////////////////////////////////
void Environment::freeLattice(const std::string name)
{
    if (hasLattice(name))
    {
        LOG(Message) << "freeing lattice '" << name << "'" << std::endl;
        lattice_.erase(name);
        object_.erase(name);
    }
    else
    {
        HADRON_ERROR("trying to free unknown lattice '" + name + "'");
    }
}

bool Environment::hasLattice(const std::string name) const
{
    return (hasObject(name) and (lattice_.find(name) != lattice_.end()));
}

// general memory management ///////////////////////////////////////////////////
bool Environment::hasObject(const std::string name) const
{
    return (object_.find(name) != object_.end());
}

void Environment::registerObject(const std::string name,
                                 const unsigned int size, const unsigned int Ls)
{
    if (!hasObject(name))
    {
        ObjInfo info{size, Ls};
        
        object_[name] = info;
    }
    else
    {
        HADRON_ERROR("object '" + name + "' already exists");
    }
}

unsigned int Environment::getObjectSize(const std::string name) const
{
    if (hasObject(name))
    {
        return object_.at(name).size;
    }
    else
    {
        HADRON_ERROR("no object named '" + name + "'");
    }
}

long unsigned int Environment::getTotalSize(void) const
{
    long unsigned int size = 0;
    
    for (auto &s: object_)
    {
        size += s.second.size;
    }
    
    return size;
}

unsigned int Environment::getObjectLs(const std::string name) const
{
    if (hasObject(name))
    {
        return object_.at(name).Ls;
    }
    else
    {
        HADRON_ERROR("no object named '" + name + "'");
    }
}

bool Environment::isObject5d(const std::string name) const
{
    return (getObjectLs(name) > 1);
}

void Environment::addOwnership(const std::string owner,
                               const std::string property)
{
    owners_[property].insert(owner);
    properties_[owner].insert(property);
}

bool Environment::hasOwners(const std::string name) const
{
    try
    {
        return (!owners_.at(name).empty());
    }
    catch (std::out_of_range &)
    {
        return false;
    }
}

bool Environment::freeObject(const std::string name)
{
    if (!hasOwners(name))
    {
        for (auto &p: properties_[name])
        {
            owners_[p].erase(name);
        }
        properties_[name].clear();
        if (hasLattice(name))
        {
            freeLattice(name);
        }
        else if (hasFermionMatrix(name))
        {
            freeFermionMatrix(name);
        }
        else if (hasObject(name))
        {
            object_.erase(name);
        }
        
        return true;
    }
    else
    {
        return false;
    }
}

void Environment::freeAll(void)
{
    object_.clear();
    lattice_.clear();
    fMat_.clear();
    solver_.clear();
    solverAction_.clear();
    owners_.clear();
    properties_.clear();
}
