/*******************************************************************************
Grid physics library, www.github.com/paboyle/Grid 

Source file: programs/Hadrons/Application.cc

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

#include <Hadrons/Application.hpp>
#include <Hadrons/Graph.hpp>

using namespace Grid;
using namespace Hadrons;

/******************************************************************************
 *                       Application implementation                           *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
Application::Application(const std::string parameterFileName)
: parameterFileName_(parameterFileName)
, env_(Environment::getInstance())
, modFactory_(ModuleFactory::getInstance())
{
    LOG(Message) << "Modules available:" << std::endl;
    auto list = modFactory_.getModuleList();
    for (auto &m: list)
    {
        LOG(Message) << "  " << m << std::endl;
    }
}

// destructor //////////////////////////////////////////////////////////////////
Application::~Application(void)
{}

// execute /////////////////////////////////////////////////////////////////////
void Application::run(void)
{
    parseParameterFile();
    schedule();
    configLoop();
}

// parse parameter file ////////////////////////////////////////////////////////
class ModuleId: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(ModuleId,
                                    std::string, name,
                                    std::string, type);
};

void Application::parseParameterFile(void)
{
    XmlReader reader(parameterFileName_);
    ModuleId  id;
    
    LOG(Message) << "Reading '" << parameterFileName_ << "'..." << std::endl;
    read(reader, "parameters", par_);
    push(reader, "modules");
    push(reader, "module");
    do
    {
        read(reader, "id", id);
        module_[id.name] = modFactory_.create(id.type, id.name);
        module_[id.name]->parseParameters(reader, "options");
        std::vector<std::string> output = module_[id.name]->getOutput();
        for (auto &n: output)
        {
            associatedModule_[n] = id.name;
        }
        input_[id.name] = module_[id.name]->getInput();
    } while (reader.nextElement("module"));
    pop(reader);
    pop(reader);
}

// schedule computation ////////////////////////////////////////////////////////
void Application::schedule(void)
{
    Graph<std::string> moduleGraph;
    
    LOG(Message) << "Scheduling computation..." << std::endl;
    
    // create dependency graph
    for (auto &m: module_)
    {
        std::vector<std::string> input = m.second->getInput();
        for (auto &n: input)
        {
            try
            {
                moduleGraph.addEdge(associatedModule_.at(n), m.first);
            }
            catch (std::out_of_range &)
            {
                HADRON_ERROR("unknown object '" + n + "'");
            }
        }
    }
    
    // topological sort
    unsigned int k = 0;
    
    std::vector<Graph<std::string>> con = moduleGraph.getConnectedComponents();
    LOG(Message) << "Program:" << std::endl;
    for (unsigned int i = 0; i < con.size(); ++i)
    {
        std::vector<std::vector<std::string>> t = con[i].allTopoSort();
        int                                   memPeak, minMemPeak = -1;
        unsigned int                          bestInd;
        bool                                  msg;

        env_.dryRun(true);
        for (unsigned int p = 0; p < t.size(); ++p)
        {
            msg = HadronsLogMessage.isActive();
            HadronsLogMessage.Active(false);
            memPeak = execute(t[p]);
            if ((memPeak < minMemPeak) or (minMemPeak < 0))
            {
                minMemPeak = memPeak;
                bestInd = p;
            }
            HadronsLogMessage.Active(msg);
            env_.freeAll();
        }
        env_.dryRun(false);
        for (unsigned int j = 0; j < t[bestInd].size(); ++j)
        {
            program_.push_back(t[bestInd][j]);
            LOG(Message) << std::setw(4) << std::right << k << ": "
                         << program_[k] << std::endl;
            k++;
        }
    }
}

// program execution ///////////////////////////////////////////////////////////
void Application::configLoop(void)
{
    auto range = par_.configs.range;
    
    for (unsigned int t = range.start; t < range.end; t += range.step)
    {
        LOG(Message) << "Starting measurement for trajectory " << t
                     << std::endl;
        env_.loadUnitGauge();
        execute(program_);
        env_.freeAll();
    }
}

unsigned int Application::execute(const std::vector<std::string> &program)
{
    unsigned int memPeak = 0;
    
    for (unsigned int i = 0; i < program.size(); ++i)
    {
        LOG(Message) << "Measurement step (" << i+1 << "/" << program.size()
                     << ")" << std::endl;
        (*module_[program[i]])(env_);
        LOG(Message) << "allocated propagators: " << env_.nProp() << std::endl;
        if (env_.nProp() > memPeak)
        {
            memPeak = env_.nProp();
        }
        for (auto &n: associatedModule_)
        {
            bool canFree = true;
            
            for (unsigned int j = i + 1; j < program.size(); ++j)
            {
                auto &in = input_[program[j]];
                auto it  = std::find(in.begin(), in.end(), n.first);
                canFree  = canFree and (it == in.end());
            }
            if (canFree and env_.propExists(n.first))
            {
                LOG(Message) << "freeing '" << n.first << "'" << std::endl;
                env_.free(n.first);
            }
        }
    }
    
    return memPeak;
}
