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
#include <Hadrons/GeneticScheduler.hpp>

using namespace Grid;
using namespace QCD;
using namespace Hadrons;

#define BIG_SEP "==============="
#define SEP     "---------------"

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
#define MEM_MSG(size)\
sizeString((size)*locVol_) << " (" << sizeString(size)  << "/site)"

void Application::schedule(void)
{
    // memory peak function
    auto memPeak = [this](const std::vector<std::string> &program)
    {
        unsigned int memPeak;
        bool         msg;
        
        msg = HadronsLogMessage.isActive();
        HadronsLogMessage.Active(false);
        env_.dryRun(true);
        memPeak = execute(program);
        env_.dryRun(false);
        env_.freeAll();
        HadronsLogMessage.Active(true);
        
        return memPeak;
    };
    
    // create dependency graph
    Graph<std::string> moduleGraph;
    
    LOG(Message) << "Scheduling computation..." << std::endl;
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
    
    // constrained topological sort using a genetic algorithm
    constexpr unsigned int maxGen = 200, maxCstGen = 50;
    unsigned int           k = 0, gen, prevPeak, nCstPeak = 0;
    std::vector<Graph<std::string>> con = moduleGraph.getConnectedComponents();
    GeneticScheduler<std::string>::Parameters par;
    std::random_device rd;

    par.popSize      = 20;
    par.mutationRate = .1;
    par.seed         = rd();
    CartesianCommunicator::BroadcastWorld(0, &(par.seed), sizeof(par.seed));
    for (unsigned int i = 0; i < con.size(); ++i)
    {
        GeneticScheduler<std::string> scheduler(con[i], memPeak, par);
        
        gen = 0;
        do
        {
            scheduler.nextGeneration();
            if (gen != 0)
            {
                if (prevPeak == scheduler.getMinValue())
                {
                    nCstPeak++;
                }
                else
                {
                    nCstPeak = 0;
                }
            }
            prevPeak = scheduler.getMinValue();
            if (gen % 10 == 0)
            {
                LOG(Iterative) << "Generation " << gen << ": "
                               << MEM_MSG(scheduler.getMinValue()) << std::endl;
            }
            gen++;
        } while ((gen < maxGen) and (nCstPeak < maxCstGen));
        auto &t = scheduler.getMinSchedule();
        LOG(Message) << "Program " << i + 1 << " (memory peak: "
                     << MEM_MSG(scheduler.getMinValue()) << "):" << std::endl;
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
        LOG(Message) << BIG_SEP << " Starting measurement for trajectory " << t
                     << " " << BIG_SEP << std::endl;
        env_.setTrajectory(t);
        execute(program_);
        env_.freeAll();
    }
    LOG(Message) << BIG_SEP << " End of measurement " << BIG_SEP << std::endl;
}

unsigned int Application::execute(const std::vector<std::string> &program)
{
    unsigned int                       memPeak = 0, sizeBefore, sizeAfter;
    std::vector<std::set<std::string>> freeProg;
    bool                               continueCollect, nothingFreed;
    
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
        LOG(Message) << SEP << " Measurement step " << i+1 << "/"
                     << program.size() << " (module '" << program[i] << "') "
                     << SEP << std::endl;
        (*module_[program[i]])();
        sizeBefore = env_.getTotalSize();
        // print used memory after execution
        LOG(Message) << "Allocated objects: " << MEM_MSG(sizeBefore)
                     << std::endl;
        if (sizeBefore > memPeak)
        {
            memPeak = sizeBefore;
        }
        // garbage collection for step i
        LOG(Message) << "Garbage collection..." << std::endl;
        nothingFreed = true;
        do
        {
            continueCollect = false;
            auto toFree = freeProg[i];
            for (auto &n: toFree)
            {
                // continue garbage collection while there are still
                // objects without owners
                continueCollect = continueCollect or !env_.hasOwners(n);
                if(env_.freeObject(n))
                {
                    // if an object has been freed, remove it from
                    // the garbage collection schedule
                    freeProg[i].erase(n);
                    nothingFreed = false;
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
        // print used memory after garbage collection if necessary
        sizeAfter = env_.getTotalSize();
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
