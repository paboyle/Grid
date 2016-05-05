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
using namespace QCD;
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
    auto list = modFactory_.getBuilderList();
    for (auto &m: list)
    {
        LOG(Message) << "  " << m << std::endl;
    }
    auto dim = GridDefaultLatt(), mpi = GridDefaultMpi(), loc(dim);
    locVol_ = 1;
    for (unsigned int d = 0; d < dim.size(); ++d)
    {
        loc[d]  /= mpi[d];
        locVol_ *= loc[d];
    }
    LOG(Message) << "Global lattice: " << dim << std::endl;
    LOG(Message) << "MPI partition : " << mpi << std::endl;
    LOG(Message) << "Local lattice : " << loc << std::endl;
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
class ObjectId: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(ObjectId,
                                    std::string, name,
                                    std::string, type);
};

void Application::parseParameterFile(void)
{
    XmlReader reader(parameterFileName_);
    ObjectId  id;
    
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
    env_.setSeed(strToVec<int>(par_.seed));
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
    
    for (unsigned int i = 0; i < con.size(); ++i)
    {
//        std::vector<std::vector<std::string>> t = con[i].allTopoSort();
//        int                                   memPeak, minMemPeak = -1;
//        unsigned int                          bestInd;
//        bool                                  msg;
//
//        LOG(Message) << "analyzing " << t.size() << " possible programs..."
//                     << std::endl;
//        env_.dryRun(true);
//        for (unsigned int p = 0; p < t.size(); ++p)
//        {
//            msg = HadronsLogMessage.isActive();
//            HadronsLogMessage.Active(false);
//        
//            memPeak = execute(t[p]);
//            if ((memPeak < minMemPeak) or (minMemPeak < 0))
//            {
//                minMemPeak = memPeak;
//                bestInd = p;
//            }
//            HadronsLogMessage.Active(msg);
//            env_.freeAll();
//        }
//        env_.dryRun(false);
        std::vector<std::string> t = con[i].topoSort();
        LOG(Message) << "Program " << i + 1 << ":" << std::endl;
        for (unsigned int j = 0; j < t.size(); ++j)
        {
            program_.push_back(t[j]);
            LOG(Message) << std::setw(4) << std::right << k + 1 << ": "
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
        LOG(Message) << "========== Starting measurement for trajectory " << t
                     << " ==========" << std::endl;
        env_.setTrajectory(t);
        execute(program_);
        env_.freeAll();
    }
}

unsigned int Application::execute(const std::vector<std::string> &program)
{
    unsigned int                       memPeak = 0, size;
    std::vector<std::set<std::string>> freeProg;
    bool                               continueCollect;
    
    // build garbage collection schedule
    freeProg.resize(program.size());
    for (auto &n: associatedModule_)
    {
        auto pred = [&n, this](const std::string &s)
        {
            auto &in = input_[s];
            auto it  = std::find(in.begin(), in.end(), n.first);
            
            return (it != in.end()) or (s == n.second);
        };
        auto it = std::find_if(program.rbegin(), program.rend(), pred);
        if (it != program.rend())
        {
            freeProg[program.rend() - it - 1].insert(n.first);
        }
    }
    // program execution
    for (unsigned int i = 0; i < program.size(); ++i)
    {
        // execute module
        LOG(Message) << "---------- Measurement step " << i+1 << "/"
                     << program.size() << " (module '" << program[i] << "')"
                     << " ----------" << std::endl;
        (*module_[program[i]])();
        size = env_.getTotalSize();
        // print used memory after execution
        LOG(Message) << "Allocated objects: " << sizeString(size*locVol_)
                     << " (" << sizeString(size)  << "/site)" << std::endl;
        if (size > memPeak)
        {
            memPeak = size;
        }
        // garbage collection for step i
        LOG(Message) << "Garbage collection..." << std::endl;
        do
        {
            continueCollect = false;
            auto toFree = freeProg[i];
            for (auto &n: toFree)
            {
                // continue garbage collection while there are still
                // objects without owners
                continueCollect = continueCollect or !env_.hasOwners(n);
                if(env_.free(n))
                {
                    // if an object has been freed, remove it from
                    // the garbage collection schedule
                    freeProg[i].erase(n);
                }
            }
        } while (continueCollect);
        // any remaining objects in step i garbage collection schedule
        // is scheduled for step i + 1
        if (i + 1 < program.size())
        {
            for (auto &n: freeProg[i])
            {
                freeProg[i + 1].insert(n);
            }
        }
        // print used memory after garbage collection
        size = env_.getTotalSize();
        LOG(Message) << "Allocated objects: " << sizeString(size*locVol_)
                     << " (" << sizeString(size)  << "/site)" << std::endl;
    }
    
    return memPeak;
}

// pretty size formatting //////////////////////////////////////////////////////
std::string Application::sizeString(long unsigned int bytes)

{
    constexpr unsigned int bufSize = 256;
    const char             *suffixes[7] = {"", "K", "M", "G", "T", "P", "E"};
    char                   buf[256];
    long unsigned int      s     = 0;
    double                 count = bytes;
    
    while (count >= 1024 && s < 7)
    {
        s++;
        count /= 1024;
    }
    if (count - floor(count) == 0.0)
    {
        snprintf(buf, bufSize, "%d %sB", (int)count, suffixes[s]);
    }
    else
    {
        snprintf(buf, bufSize, "%.1f %sB", count, suffixes[s]);
    }
    
    return std::string(buf);
}
