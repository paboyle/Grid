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
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

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
    auto loc = getGrid()->LocalDimensions();
    locVol_ = 1;
    for (unsigned int d = 0; d < loc.size(); ++d)
    {
        locVol_ *= loc[d];
    }
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
    if (grid5d_.find(Ls) == grid5d_.end())
    {
        auto g = getGrid();
        
        grid5d_[Ls].reset(SpaceTimeGrid::makeFiveDimGrid(Ls, g));
        gridRb5d_[Ls].reset(SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls, g));
    }
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

// random number generator /////////////////////////////////////////////////////
void Environment::setSeed(const std::vector<int> &seed)
{
    rng4d_->SeedFixedIntegers(seed);
}

GridParallelRNG * Environment::get4dRng(void) const
{
    return rng4d_.get();
}

// module management ///////////////////////////////////////////////////////////
void Environment::createModule(const std::string name, const std::string type,
                               XmlReader &reader)
{
    auto addObject = [this](const std::string name, const int moduleAddress)
    {
        ObjInfo info;
        
        object_.push_back(info);
        objectName_.push_back(name);
        objectAddress_[name] = object_.size() - 1;
        objectModule_.push_back(moduleAddress);
        owners_.push_back(std::set<unsigned int>());
        properties_.push_back(std::set<unsigned int>());
    };
    
    if (!hasModule(name))
    {
        auto                      &factory = ModuleFactory::getInstance();
        std::vector<unsigned int> inputAddress;

        module_.push_back(factory.create(type, name));
        moduleType_.push_back(type);
        moduleName_.push_back(name);
        moduleAddress_[name] = module_.size() - 1;
        module_.back()->parseParameters(reader, "options");
        auto input  = module_.back()->getInput();
        for (auto &in: input)
        {
            if (!hasObject(in))
            {
                addObject(in , -1);
            }
            inputAddress.push_back(objectAddress_[in]);
        }
        moduleInput_.push_back(inputAddress);
        auto output = module_.back()->getOutput();
        for (auto &out: output)
        {
            if (!hasObject(out))
            {
                addObject(out , module_.size() - 1);
            }
            else
            {
                if (objectModule_[objectAddress_[out]] < 0)
                {
                    objectModule_[objectAddress_[out]] = module_.size() - 1;
                }
                else
                {
                    HADRON_ERROR("object '" + out
                             + "' is already produced by module '"
                             + moduleName_[objectModule_[getObjectAddress(out)]]
                             + "' (while creating module '" + name + "')");
                }
            }
        }
    }
    else
    {
        HADRON_ERROR("module '" + name + "' already exists");
    }
}

ModuleBase * Environment::getModule(const unsigned int address) const
{
    if (hasModule(address))
    {
        return module_[address].get();
    }
    else
    {
        HADRON_ERROR("no module with address " + std::to_string(address));
    }
}

ModuleBase * Environment::getModule(const std::string name) const
{
    return getModule(getModuleAddress(name));
}

unsigned int Environment::getModuleAddress(const std::string name) const
{
    if (hasModule(name))
    {
        return moduleAddress_.at(name);
    }
    else
    {
        HADRON_ERROR("no module with name '" + name + "'");
    }
}

std::string Environment::getModuleName(const unsigned int address) const
{
    if (hasModule(address))
    {
        return moduleName_[address];
    }
    else
    {
        HADRON_ERROR("no module with address " + std::to_string(address));
    }
}

std::string Environment::getModuleType(const unsigned int address) const
{
    if (hasModule(address))
    {
        return moduleType_[address];
    }
    else
    {
        HADRON_ERROR("no module with address " + std::to_string(address));
    }
}

std::string Environment::getModuleType(const std::string name) const
{
    return getModuleType(getModuleAddress(name));
}

bool Environment::hasModule(const unsigned int address) const
{
    return (address < module_.size());
}

bool Environment::hasModule(const std::string name) const
{
    return (moduleAddress_.find(name) != moduleAddress_.end());
}

Graph<unsigned int> Environment::makeModuleGraph(void) const
{
    Graph<unsigned int> moduleGraph;
    
    for (unsigned int i = 0; i < module_.size(); ++i)
    {
        moduleGraph.addVertex(i);
        for (auto &j: moduleInput_[i])
        {
            moduleGraph.addEdge(objectModule_[j], i);
        }
    }
    
    return moduleGraph;
}

#define BIG_SEP "==============="
#define SEP     "---------------"
#define MEM_MSG(size)\
sizeString((size)*locVol_) << " (" << sizeString(size)  << "/site)"

unsigned int Environment::executeProgram(const std::vector<unsigned int> &p)
{
    unsigned int                        memPeak = 0, sizeBefore, sizeAfter;
    std::vector<std::set<unsigned int>> freeProg;
    bool                                continueCollect, nothingFreed;
    
    // build garbage collection schedule
    freeProg.resize(p.size());
    for (unsigned int i = 0; i < object_.size(); ++i)
    {
        auto pred = [i, this](const unsigned int j)
        {
            auto &in = moduleInput_[j];
            auto it  = std::find(in.begin(), in.end(), i);
            
            return (it != in.end()) or (j == objectModule_[i]);
        };
        auto it = std::find_if(p.rbegin(), p.rend(), pred);
        if (it != p.rend())
        {
            freeProg[p.rend() - it - 1].insert(i);
        }
    }
    
    // program execution
    for (unsigned int i = 0; i < p.size(); ++i)
    {
        // execute module
        if (!isDryRun())
        {
            LOG(Message) << SEP << " Measurement step " << i+1 << "/"
                         << p.size() << " (module '" << moduleName_[p[i]]
                         << "') " << SEP << std::endl;
        }
        (*module_[p[i]])();
        sizeBefore = getTotalSize();
        // print used memory after execution
        if (!isDryRun())
        {
            LOG(Message) << "Allocated objects: " << MEM_MSG(sizeBefore)
                         << std::endl;
        }
        if (sizeBefore > memPeak)
        {
            memPeak = sizeBefore;
        }
        // garbage collection for step i
        if (!isDryRun())
        {
            LOG(Message) << "Garbage collection..." << std::endl;
        }
        nothingFreed = true;
        do
        {
            continueCollect = false;
            auto toFree = freeProg[i];
            for (auto &j: toFree)
            {
                // continue garbage collection while there are still
                // objects without owners
                continueCollect = continueCollect or !hasOwners(j);
                if(freeObject(j))
                {
                    // if an object has been freed, remove it from
                    // the garbage collection schedule
                    freeProg[i].erase(j);
                    nothingFreed = false;
                }
            }
        } while (continueCollect);
        // any remaining objects in step i garbage collection schedule
        // is scheduled for step i + 1
        if (i + 1 < p.size())
        {
            for (auto &j: freeProg[i])
            {
                freeProg[i + 1].insert(j);
            }
        }
        // print used memory after garbage collection if necessary
        sizeAfter = getTotalSize();
        if (!isDryRun())
        {
            if (sizeBefore != sizeAfter)
            {
                LOG(Message) << "Allocated objects: " << MEM_MSG(sizeAfter)
                             << std::endl;
            }
            else
            {
                LOG(Message) << "Nothing to free" << std::endl;
            }
        }
    }
    
    return memPeak;
}

unsigned int Environment::executeProgram(const std::vector<std::string> &p)
{
    std::vector<unsigned int> pAddress;
    
    for (auto &n: p)
    {
        pAddress.push_back(getModuleAddress(n));
    }
    
    return executeProgram(pAddress);
}

// lattice store ///////////////////////////////////////////////////////////////
void Environment::freeLattice(const unsigned int address)
{
    if (hasLattice(address))
    {
        if (!isDryRun())
        {
            LOG(Message) << "Freeing lattice '" << moduleName_[address]
                         << "'" << std::endl;
        }
        lattice_.erase(address);
        object_[address] = ObjInfo();
    }
    else
    {
        HADRON_ERROR("trying to free unknown lattice (address "
                     + std::to_string(address) + ")");
    }
}

bool Environment::hasLattice(const unsigned int address) const
{
    return (hasRegisteredObject(address)
            and (lattice_.find(address) != lattice_.end()));
}

bool Environment::hasLattice(const std::string name) const
{
    if (hasObject(name))
    {
        return hasLattice(getObjectAddress(name));
    }
    else
    {
        return false;
    }
}

// fermion actions /////////////////////////////////////////////////////////////
void Environment::addFermionMatrix(const std::string name, FMat *fMat)
{
    if (hasRegisteredObject(name))
    {
        fMat_[getObjectAddress(name)].reset(fMat);
    }
    else
    {
        HADRON_ERROR("no object named '" << name << "'");
    }
}

Environment::FMat * Environment::getFermionMatrix(const std::string name) const
{
    unsigned int i;
    
    if (hasFermionMatrix(name))
    {
        i = getObjectAddress(name);
        
        return fMat_.at(i).get();
    }
    else
    {
        if (hasSolver(name))
        {
            i = getObjectAddress(solverAction_.at(name));
            
            return fMat_.at(i).get();
        }
        else
        {
            HADRON_ERROR("no action/solver with name '" << name << "'");
        }
    }
}

bool Environment::hasFermionMatrix(const unsigned int address) const
{
    return (hasRegisteredObject(address)
            and (fMat_.find(address) != fMat_.end()));
}

bool Environment::hasFermionMatrix(const std::string name) const
{
    if (hasObject(name))
    {
        return hasFermionMatrix(getObjectAddress(name));
    }
    else
    {
        return false;
    }
}

void Environment::freeFermionMatrix(const unsigned int address)
{
    if (hasFermionMatrix(address))
    {
        if (!isDryRun())
        {
            LOG(Message) << "Freeing fermion matrix '" << objectName_[address]
                         << "'" << std::endl;
        }
        fMat_.erase(address);
        object_[address] = ObjInfo();
    }
    else
    {
        HADRON_ERROR("trying to free unknown fermion matrix (address "
                     + std::to_string(address) + ")");
    }
}

void Environment::freeFermionMatrix(const std::string name)
{
    freeFermionMatrix(getObjectAddress(name));
}

// solvers /////////////////////////////////////////////////////////////////////
void Environment::addSolver(const std::string name, Solver s)
{
    auto address = getObjectAddress(name);
    
    if (hasRegisteredObject(address))
    {
        solver_[address] = s;
    }
    else
    {
        HADRON_ERROR("object with name '" + name
                     + "' exsists but is not registered");
    }
}

bool Environment::hasSolver(const unsigned int address) const
{
    return (hasRegisteredObject(address)
            and (solver_.find(address) != solver_.end()));
}

bool Environment::hasSolver(const std::string name) const
{
    if (hasObject(name))
    {
        return hasSolver(getObjectAddress(name));
    }
    else
    {
        return false;
    }
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
        solver_.at(getObjectAddress(name))(sol, source);
    }
    else
    {
        HADRON_ERROR("no solver with name '" << name << "'");
    }
}

// general memory management ///////////////////////////////////////////////////
void Environment::registerObject(const unsigned int address,
                                 const unsigned int size, const unsigned int Ls)
{
    if (!hasRegisteredObject(address))
    {
        if (hasObject(address))
        {
            ObjInfo info;
            
            info.size         = size;
            info.Ls           = Ls;
            info.isRegistered = true;
            object_[address]  = info;
        }
        else
        {
            HADRON_ERROR("no object with address " + std::to_string(address));
        }
    }
    else
    {
        HADRON_ERROR("object with address " + std::to_string(address)
                     + " already registered");
    }
}

void Environment::registerObject(const std::string name,
                                 const unsigned int size, const unsigned int Ls)
{
    registerObject(getObjectAddress(name), size, Ls);
}

unsigned int Environment::getObjectAddress(const std::string name) const
{
    if (hasObject(name))
    {
        return objectAddress_.at(name);
    }
    else
    {
        HADRON_ERROR("no object with name '" + name + "'");
    }
}

std::string Environment::getObjectName(const unsigned int address) const
{
    if (hasObject(address))
    {
        return objectName_[address];
    }
    else
    {
        HADRON_ERROR("no object with address " + std::to_string(address));
    }
}

unsigned int Environment::getObjectSize(const unsigned int address) const
{
    if (hasRegisteredObject(address))
    {
        return object_[address].size;
    }
    else if (hasObject(address))
    {
        HADRON_ERROR("object with address " + std::to_string(address)
                     + " exsists but is not registered");
    }
    else
    {
        HADRON_ERROR("no object with address " + std::to_string(address));
    }
}

unsigned int Environment::getObjectSize(const std::string name) const
{
    return getObjectSize(getObjectAddress(name));
}

unsigned int Environment::getObjectLs(const unsigned int address) const
{
    if (hasRegisteredObject(address))
    {
        return object_[address].Ls;
    }
    else if (hasObject(address))
    {
        HADRON_ERROR("object with address " + std::to_string(address)
                     + " exsists but is not registered");
    }
    else
    {
        HADRON_ERROR("no object with address " + std::to_string(address));
    }
}

unsigned int Environment::getObjectLs(const std::string name) const
{
    return getObjectLs(getObjectAddress(name));
}

bool Environment::hasObject(const unsigned int address) const
{
    return (address < object_.size());
}

bool Environment::hasObject(const std::string name) const
{
    auto it = objectAddress_.find(name);
    
    return ((it != objectAddress_.end()) and hasObject(it->second));
}

bool Environment::hasRegisteredObject(const unsigned int address) const
{
    if (hasObject(address))
    {
        return object_[address].isRegistered;
    }
    else
    {
        return false;
    }
}

bool Environment::hasRegisteredObject(const std::string name) const
{
    if (hasObject(name))
    {
        return hasRegisteredObject(getObjectAddress(name));
    }
    else
    {
        return false;
    }
}

bool Environment::isObject5d(const unsigned int address) const
{
    return (getObjectLs(address) > 1);
}

bool Environment::isObject5d(const std::string name) const
{
    return (getObjectLs(name) > 1);
}

long unsigned int Environment::getTotalSize(void) const
{
    long unsigned int size = 0;
    
    for (auto &o: object_)
    {
        if (o.isRegistered)
        {
            size += o.size;
        }
    }
    
    return size;
}

void Environment::addOwnership(const unsigned int owner,
                               const unsigned int property)
{
    owners_[property].insert(owner);
    properties_[owner].insert(property);
}

void Environment::addOwnership(const std::string owner,
                               const std::string property)
{
    addOwnership(getObjectAddress(owner), getObjectAddress(property));
}

bool Environment::hasOwners(const unsigned int address) const
{
    
    if (hasObject(address))
    {
        return (!owners_[address].empty());
    }
    else
    {
        HADRON_ERROR("no object with address " + std::to_string(address));
    }
}

bool Environment::hasOwners(const std::string name) const
{
    return hasOwners(getObjectAddress(name));
}

bool Environment::freeObject(const unsigned int address)
{
    if (!hasOwners(address))
    {
        for (auto &p: properties_[address])
        {
            owners_[p].erase(address);
        }
        properties_[address].clear();
        if (hasLattice(address))
        {
            freeLattice(address);
        }
        else if (hasFermionMatrix(address))
        {
            freeFermionMatrix(address);
        }
        else if (hasObject(address))
        {
            object_[address] = ObjInfo();
        }
        
        return true;
    }
    else
    {
        return false;
    }
}

bool Environment::freeObject(const std::string name)
{
    return freeObject(getObjectAddress(name));
}

void Environment::freeAll(void)
{
    lattice_.clear();
    fMat_.clear();
    solver_.clear();
    solverAction_.clear();
    owners_.clear();
    properties_.clear();
}

void Environment::printContent(void)
{
    LOG(Message) << "Modules: " << std::endl;
    for (unsigned int i = 0; i < module_.size(); ++i)
    {
        LOG(Message) << std::setw(4) << std::right << i << ": "
                     << moduleName_[i] << " ("
                     << moduleType_[i] << ")" << std::endl;
    }
    LOG(Message) << "Objects: " << std::endl;
    for (unsigned int i = 0; i < object_.size(); ++i)
    {
        LOG(Message) << std::setw(4) << std::right << i << ": "
                     << objectName_[i] << std::endl;
    }
}
