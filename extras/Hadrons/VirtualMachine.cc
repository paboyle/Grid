/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: extras/Hadrons/VirtualMachine.cc

Copyright (C) 2017

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

#include <Grid/Hadrons/VirtualMachine.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>

using namespace Grid;
using namespace QCD;
using namespace Hadrons;

/******************************************************************************
 *                      VirtualMachine implementation                         *
 ******************************************************************************/
// dry run /////////////////////////////////////////////////////////////////////
void VirtualMachine::dryRun(const bool isDry)
{
    dryRun_ = isDry;
}

bool VirtualMachine::isDryRun(void) const
{
    return dryRun_;
}

void VirtualMachine::memoryProfile(const bool doMemoryProfile)
{
    memoryProfile_ = doMemoryProfile;
}

bool VirtualMachine::doMemoryProfile(void) const
{
    return memoryProfile_;
}

// trajectory counter //////////////////////////////////////////////////////////
void VirtualMachine::setTrajectory(const unsigned int traj)
{
    traj_ = traj;
}

unsigned int VirtualMachine::getTrajectory(void) const
{
    return traj_;
}

// module management ///////////////////////////////////////////////////////////
void VirtualMachine::pushModule(VirtualMachine::ModPt &pt)
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
        for (auto &in: m.data->getInput())
        {
            if (!env().hasObject(in))
            {
                env().addObject(in , -1);
            }
            m.input.push_back(env().getObjectAddress(in));
        }
        for (auto &ref: m.data->getReference())
        {
            if (!env().hasObject(ref))
            {
                env().addObject(ref , -1);
            }
            m.input.push_back(env().getObjectAddress(ref));
        }
        module_.push_back(std::move(m));
        address              = static_cast<unsigned int>(module_.size() - 1);
        moduleAddress_[name] = address;
        for (auto &out: getModule(address)->getOutput())
        {
            if (!env().hasObject(out))
            {
                env().addObject(out, address);
            }
            else
            {
                if (env().getObjectModule(env().getObjectAddress(out)) < 0)
                {
                    env().setObjectModule(env().getObjectAddress(out), address);
                }
                else
                {
                    HADRON_ERROR(Definition, "object '" + out
                                 + "' is already produced by module '"
                                 + module_[env().getObjectModule(out)].name
                                 + "' (while pushing module '" + name + "')");
                }
                if (getModule(address)->getReference().size() > 0)
                {
                    auto pred = [this, out](const ModuleInfo &n)
                    {
                        auto &in = n.input;
                        auto it  = std::find(in.begin(), in.end(), env().getObjectAddress(out));
                        
                        return (it != in.end());
                    };
                    auto it = std::find_if(module_.begin(), module_.end(), pred);
                    while (it != module_.end())
                    {
                        for (auto &ref: getModule(address)->getReference())
                        {
                            it->input.push_back(env().getObjectAddress(ref));
                        }
                        it = std::find_if(++it, module_.end(), pred);
                    }   
                }
            }
        }
    }
    else
    {
        HADRON_ERROR(Definition, "module '" + name + "' already exists");
    }
}

unsigned int VirtualMachine::getNModule(void) const
{
    return module_.size();
}

void VirtualMachine::createModule(const std::string name, const std::string type,
                               XmlReader &reader)
{
    auto &factory = ModuleFactory::getInstance();
    auto pt       = factory.create(type, name);
    
    pt->parseParameters(reader, "options");
    pushModule(pt);
}

ModuleBase * VirtualMachine::getModule(const unsigned int address) const
{
    if (hasModule(address))
    {
        return module_[address].data.get();
    }
    else
    {
        HADRON_ERROR(Definition, "no module with address " + std::to_string(address));
    }
}

ModuleBase * VirtualMachine::getModule(const std::string name) const
{
    return getModule(getModuleAddress(name));
}

unsigned int VirtualMachine::getModuleAddress(const std::string name) const
{
    if (hasModule(name))
    {
        return moduleAddress_.at(name);
    }
    else
    {
        HADRON_ERROR(Definition, "no module with name '" + name + "'");
    }
}

std::string VirtualMachine::getModuleName(const unsigned int address) const
{
    if (hasModule(address))
    {
        return module_[address].name;
    }
    else
    {
        HADRON_ERROR(Definition, "no module with address " + std::to_string(address));
    }
}

std::string VirtualMachine::getModuleType(const unsigned int address) const
{
    if (hasModule(address))
    {
        return typeName(module_[address].type);
    }
    else
    {
        HADRON_ERROR(Definition, "no module with address " + std::to_string(address));
    }
}

std::string VirtualMachine::getModuleType(const std::string name) const
{
    return getModuleType(getModuleAddress(name));
}

std::string VirtualMachine::getModuleNamespace(const unsigned int address) const
{
    std::string type = getModuleType(address), ns;
    
    auto pos2 = type.rfind("::");
    auto pos1 = type.rfind("::", pos2 - 2);
    
    return type.substr(pos1 + 2, pos2 - pos1 - 2);
}

std::string VirtualMachine::getModuleNamespace(const std::string name) const
{
    return getModuleNamespace(getModuleAddress(name));
}

bool VirtualMachine::hasModule(const unsigned int address) const
{
    return (address < module_.size());
}

bool VirtualMachine::hasModule(const std::string name) const
{
    return (moduleAddress_.find(name) != moduleAddress_.end());
}

Graph<unsigned int> VirtualMachine::makeModuleGraph(void) const
{
    Graph<unsigned int> moduleGraph;
    
    // create vertices
    for (unsigned int m = 0; m < module_.size(); ++m)
    {
        moduleGraph.addVertex(m);
    }
    // create edges
    for (unsigned int m = 0; m < module_.size(); ++m)
    {
        for (auto &in: module_[m].input)
        {
            moduleGraph.addEdge(env().getObjectModule(in), m);
        }
    }
    
    return moduleGraph;
}

// void VirtualMachine::checkGraph(void) const
// {
//     for (auto &o: object_)
//     {
//         if (o.module < 0)
//         {
//             HADRON_ERROR(Runtime, "object '" + o.name + "' does not have a creator");
//         }
//     }
// }

// general execution ///////////////////////////////////////////////////////////
#define BIG_SEP "==============="
#define SEP     "---------------"
#define MEM_MSG(size) sizeString(size)

VirtualMachine::Size
VirtualMachine::executeProgram(const std::vector<unsigned int> &p)
{
    Size                                memPeak = 0, sizeBefore, sizeAfter;
    std::vector<std::set<unsigned int>> freeProg;
    
    // build garbage collection schedule
    LOG(Debug) << "Building garbage collection schedule..." << std::endl;
    freeProg.resize(p.size());
    for (unsigned int i = 0; i < env().getMaxAddress(); ++i)
    {
        auto pred = [i, this](const unsigned int j)
        {
            auto &in = module_[j].input;
            auto it  = std::find(in.begin(), in.end(), i);
            
            return (it != in.end()) or (j == env().getObjectModule(i));
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
        sizeBefore = env().getTotalSize();
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
        for (auto &j: freeProg[i])
        {
            env().freeObject(j);
        }
        // free temporaries
        for (unsigned int i = 0; i < env().getMaxAddress(); ++i)
        {
            if ((env().getObjectStorage(i) == Environment::Storage::temporary) 
                and env().hasCreatedObject(i))
            {
                env().freeObject(i);
            }
        }
        // print used memory after garbage collection if necessary
        if (!isDryRun())
        {
            sizeAfter = env().getTotalSize();
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

VirtualMachine::Size VirtualMachine::executeProgram(const std::vector<std::string> &p)
{
    std::vector<unsigned int> pAddress;
    
    for (auto &n: p)
    {
        pAddress.push_back(getModuleAddress(n));
    }
    
    return executeProgram(pAddress);
}

// print VM content ////////////////////////////////////////////////////////////
void VirtualMachine::printContent(void) const
{
    LOG(Debug) << "Modules: " << std::endl;
    for (unsigned int i = 0; i < module_.size(); ++i)
    {
        LOG(Debug) << std::setw(4) << i << ": "
                   << getModuleName(i) << std::endl;
    }
}
