/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: extras/Hadrons/Environment.cc

Copyright (C) 2015
Copyright (C) 2016

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

#include <Grid/Hadrons/Environment.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>

using namespace Grid;
using namespace QCD;
using namespace Hadrons;

/******************************************************************************
 *                       Environment implementation                           *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
Environment::Environment(void)
{
    dim_ = GridDefaultLatt();
    nd_  = dim_.size();
    grid4d_.reset(SpaceTimeGrid::makeFourDimGrid(
        dim_, GridDefaultSimd(nd_, vComplex::Nsimd()),
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

void Environment::memoryProfile(const bool doMemoryProfile)
{
    memoryProfile_ = doMemoryProfile;
}

bool Environment::doMemoryProfile(void) const
{
    return memoryProfile_;
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

unsigned int Environment::getNd(void) const
{
    return nd_;
}

std::vector<int> Environment::getDim(void) const
{
    return dim_;
}

int Environment::getDim(const unsigned int mu) const
{
    return dim_[mu];
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
void Environment::pushModule(Environment::ModPt &pt)
{
    std::string name = pt->getName();
    
    if (!hasModule(name))
    {
        std::vector<unsigned int> inputAddress;
        unsigned int              address;
        ModuleInfo                m;
        
        m.data = std::move(pt);
        m.type = typeIdPt(*m.data.get());
        m.name = name;
        auto input  = m.data->getInput();
        for (auto &in: input)
        {
            if (!hasObject(in))
            {
                addObject(in , -1);
            }
            m.input.push_back(objectAddress_[in]);
        }
        auto output = m.data->getOutput();
        module_.push_back(std::move(m));
        address              = static_cast<unsigned int>(module_.size() - 1);
        moduleAddress_[name] = address;
        for (auto &out: output)
        {
            if (!hasObject(out))
            {
                addObject(out, address);
            }
            else
            {
                if (object_[objectAddress_[out]].module < 0)
                {
                    object_[objectAddress_[out]].module = address;
                }
                else
                {
                    HADRON_ERROR("object '" + out
                                 + "' is already produced by module '"
                                 + module_[object_[getObjectAddress(out)].module].name
                                 + "' (while pushing module '" + name + "')");
                }
            }
        }
    }
    else
    {
        HADRON_ERROR("module '" + name + "' already exists");
    }
}

unsigned int Environment::getNModule(void) const
{
    return module_.size();
}

void Environment::createModule(const std::string name, const std::string type,
                               XmlReader &reader)
{
    auto &factory = ModuleFactory::getInstance();
    auto pt       = factory.create(type, name);
    
    pt->parseParameters(reader, "options");
    pushModule(pt);
}

ModuleBase * Environment::getModule(const unsigned int address) const
{
    if (hasModule(address))
    {
        return module_[address].data.get();
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
        return module_[address].name;
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
        return typeName(module_[address].type);
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

std::string Environment::getModuleNamespace(const unsigned int address) const
{
    std::string type = getModuleType(address), ns;
    
    auto pos2 = type.rfind("::");
    auto pos1 = type.rfind("::", pos2 - 2);
    
    return type.substr(pos1 + 2, pos2 - pos1 - 2);
}

std::string Environment::getModuleNamespace(const std::string name) const
{
    return getModuleNamespace(getModuleAddress(name));
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
        for (auto &j: module_[i].input)
        {
            moduleGraph.addEdge(object_[j].module, i);
        }
    }
    
    return moduleGraph;
}

void Environment::checkGraph(void) const
{
    for (auto &o: object_)
    {
        if (o.module < 0)
        {
            HADRON_ERROR("object '" + o.name + "' does not have a creator");
        }
    }
}

#define BIG_SEP "==============="
#define SEP     "---------------"
#define MEM_MSG(size)\
sizeString((size)*locVol_) << " (" << sizeString(size)  << "/site)"

Environment::Size
Environment::executeProgram(const std::vector<unsigned int> &p)
{
    Size                                memPeak = 0, sizeBefore, sizeAfter;
    std::vector<std::set<unsigned int>> freeProg;
    bool                                continueCollect, nothingFreed;
    
    // build garbage collection schedule
    LOG(Debug) << "Building garbage collection schedule..." << std::endl;
    freeProg.resize(p.size());
    for (unsigned int i = 0; i < object_.size(); ++i)
    {
        auto pred = [i, this](const unsigned int j)
        {
            auto &in = module_[j].input;
            auto it  = std::find(in.begin(), in.end(), i);
            
            return (it != in.end()) or (j == object_[i].module);
        };
        auto it = std::find_if(p.rbegin(), p.rend(), pred);
        if (it != p.rend())
        {
            freeProg[std::distance(it, p.rend()) - 1].insert(i);
        }
    }

    // program execution
    LOG(Debug) << "Executing program..." << std::endl;
    for (unsigned int i = 0; i < p.size(); ++i)
    {
        // execute module
        if (!isDryRun())
        {
            LOG(Message) << SEP << " Measurement step " << i+1 << "/"
                         << p.size() << " (module '" << module_[p[i]].name
                         << "') " << SEP << std::endl;
        }
        (*module_[p[i]].data)();
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
        // free temporaries
        for (unsigned int i = 0; i < object_.size(); ++i)
        {
            if ((object_[i].storage == Storage::temporary) 
                and hasCreatedObject(i))
            {
                freeObject(i);
            }
        }
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
        if (!isDryRun())
        {
            sizeAfter = getTotalSize();
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

Environment::Size Environment::executeProgram(const std::vector<std::string> &p)
{
    std::vector<unsigned int> pAddress;
    
    for (auto &n: p)
    {
        pAddress.push_back(getModuleAddress(n));
    }
    
    return executeProgram(pAddress);
}

// general memory management ///////////////////////////////////////////////////
void Environment::addObject(const std::string name, const int moduleAddress)
{
    if (!hasObject(name))
    {
        ObjInfo info;
        
        info.name   = name;
        info.module = moduleAddress;
        info.data   = nullptr;
        object_.push_back(std::move(info));
        objectAddress_[name] = static_cast<unsigned int>(object_.size() - 1);
    }
    else
    {
        HADRON_ERROR("object '" + name + "' already exists");
    }
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
        return object_[address].name;
    }
    else
    {
        HADRON_ERROR("no object with address " + std::to_string(address));
    }
}

std::string Environment::getObjectType(const unsigned int address) const
{
    if (hasObject(address))
    {
        if (object_[address].type)
        {
            return typeName(object_[address].type);
        }
        else
        {
            return "<no type>";
        }
    }
    else
    {
        HADRON_ERROR("no object with address " + std::to_string(address));
    }
}

std::string Environment::getObjectType(const std::string name) const
{
    return getObjectType(getObjectAddress(name));
}

Environment::Size Environment::getObjectSize(const unsigned int address) const
{
    if (hasObject(address))
    {
        return object_[address].size;
    }
    else
    {
        HADRON_ERROR("no object with address " + std::to_string(address));
    }
}

Environment::Size Environment::getObjectSize(const std::string name) const
{
    return getObjectSize(getObjectAddress(name));
}

unsigned int Environment::getObjectModule(const unsigned int address) const
{
    if (hasObject(address))
    {
        return object_[address].module;
    }
    else
    {
        HADRON_ERROR("no object with address " + std::to_string(address));
    }
}

unsigned int Environment::getObjectModule(const std::string name) const
{
    return getObjectModule(getObjectAddress(name));
}

unsigned int Environment::getObjectLs(const unsigned int address) const
{
    if (hasObject(address))
    {
        return object_[address].Ls;
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

bool Environment::hasCreatedObject(const unsigned int address) const
{
    if (hasObject(address))
    {
        return (object_[address].data != nullptr);
    }
    else
    {
        return false;
    }
}

bool Environment::hasCreatedObject(const std::string name) const
{
    if (hasObject(name))
    {
        return hasCreatedObject(getObjectAddress(name));
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

Environment::Size Environment::getTotalSize(void) const
{
    Environment::Size size = 0;
    
    for (auto &o: object_)
    {
        size += o.size;
    }
    
    return size;
}

void Environment::addOwnership(const unsigned int owner,
                               const unsigned int property)
{
    if (hasObject(property))
    {
        object_[property].owners.insert(owner);
    }
    else
    {
        HADRON_ERROR("no object with address " + std::to_string(property));
    }
    if (hasObject(owner))
    {
        object_[owner].properties.insert(property);
    }
    else
    {
        HADRON_ERROR("no object with address " + std::to_string(owner));
    }
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
        return (!object_[address].owners.empty());
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
        if (!isDryRun() and hasCreatedObject(address))
        {
            LOG(Message) << "Destroying object '" << object_[address].name
                         << "'" << std::endl;
        }
        for (auto &p: object_[address].properties)
        {
            object_[p].owners.erase(address);
        }
        object_[address].size = 0;
        object_[address].type = nullptr;
        object_[address].owners.clear();
        object_[address].properties.clear();
        object_[address].data.reset(nullptr);
        
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
    for (unsigned int i = 0; i < object_.size(); ++i)
    {
        freeObject(i);
    }
}

void Environment::printContent(void)
{
    LOG(Debug) << "Modules: " << std::endl;
    for (unsigned int i = 0; i < module_.size(); ++i)
    {
        LOG(Debug) << std::setw(4) << i << ": "
                   << getModuleName(i) << std::endl;
    }
    LOG(Debug) << "Objects: " << std::endl;
    for (unsigned int i = 0; i < object_.size(); ++i)
    {
        LOG(Debug) << std::setw(4) << i << ": "
                   << getObjectName(i) << std::endl;
    }
}
