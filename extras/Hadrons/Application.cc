/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: extras/Hadrons/Application.cc

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

#include <Grid/Hadrons/Application.hpp>
#include <Grid/Hadrons/GeneticScheduler.hpp>

using namespace Grid;
using namespace QCD;
using namespace Hadrons;

#define BIG_SEP "==============="
#define SEP     "---------------"

/******************************************************************************
 *                       Application implementation                           *
 ******************************************************************************/
// constructors ////////////////////////////////////////////////////////////////
Application::Application(void)
{
    LOG(Message) << "Modules available:" << std::endl;
    auto list = ModuleFactory::getInstance().getBuilderList();
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

Application::Application(const Application::GlobalPar &par)
: Application()
{
    setPar(par);
}

Application::Application(const std::string parameterFileName)
: Application()
{
    parameterFileName_ = parameterFileName;
}

// environment shortcut ////////////////////////////////////////////////////////
Environment & Application::env(void) const
{
    return Environment::getInstance();
}

// access //////////////////////////////////////////////////////////////////////
void Application::setPar(const Application::GlobalPar &par)
{
    par_ = par;
    env().setSeed(strToVec<int>(par_.seed));
}

const Application::GlobalPar & Application::getPar(void)
{
    return par_;
}

// execute /////////////////////////////////////////////////////////////////////
void Application::run(void)
{
    if (!parameterFileName_.empty() and (env().getNModule() == 0))
    {
        parseParameterFile(parameterFileName_);
    }
    if (!scheduled_)
    {
        schedule();
    }
    printSchedule();
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

void Application::parseParameterFile(const std::string parameterFileName)
{
    XmlReader reader(parameterFileName);
    GlobalPar par;
    ObjectId  id;
    
    LOG(Message) << "Building application from '" << parameterFileName << "'..." << std::endl;
    read(reader, "parameters", par);
    setPar(par);
    push(reader, "modules");
    push(reader, "module");
    do
    {
        read(reader, "id", id);
        env().createModule(id.name, id.type, reader);
    } while (reader.nextElement("module"));
    pop(reader);
    pop(reader);
}

void Application::saveParameterFile(const std::string parameterFileName)
{
    XmlWriter          writer(parameterFileName);
    ObjectId           id;
    const unsigned int nMod = env().getNModule();
    
    LOG(Message) << "Saving application to '" << parameterFileName << "'..." << std::endl;
    write(writer, "parameters", getPar());
    push(writer, "modules");
    for (unsigned int i = 0; i < nMod; ++i)
    {
        push(writer, "module");
        id.name = env().getModuleName(i);
        id.type = env().getModule(i)->getRegisteredName();
        write(writer, "id", id);
        env().getModule(i)->saveParameters(writer, "options");
        pop(writer);
    }
    pop(writer);
    pop(writer);
}

// schedule computation ////////////////////////////////////////////////////////
#define MEM_MSG(size)\
sizeString((size)*locVol_) << " (" << sizeString(size)  << "/site)"

#define DEFINE_MEMPEAK \
GeneticScheduler<unsigned int>::ObjFunc memPeak = \
[this](const std::vector<unsigned int> &program)\
{\
    unsigned int memPeak;\
    bool         msg;\
    \
    msg = HadronsLogMessage.isActive();\
    HadronsLogMessage.Active(false);\
    env().dryRun(true);\
    memPeak = env().executeProgram(program);\
    env().dryRun(false);\
    env().freeAll();\
    HadronsLogMessage.Active(true);\
    \
    return memPeak;\
}

void Application::schedule(void)
{
    DEFINE_MEMPEAK;
    
    // build module dependency graph
    LOG(Message) << "Building module graph..." << std::endl;
    auto graph = env().makeModuleGraph();
    auto con = graph.getConnectedComponents();
    
    // constrained topological sort using a genetic algorithm
    LOG(Message) << "Scheduling computation..." << std::endl;
    LOG(Message) << "               #module= " << graph.size() << std::endl;
    LOG(Message) << "       population size= " << par_.genetic.popSize << std::endl;
    LOG(Message) << "       max. generation= " << par_.genetic.maxGen << std::endl;
    LOG(Message) << "  max. cst. generation= " << par_.genetic.maxCstGen << std::endl;
    LOG(Message) << "         mutation rate= " << par_.genetic.mutationRate << std::endl;
    
    unsigned int                               k = 0, gen, prevPeak, nCstPeak = 0;
    std::random_device                         rd;
    GeneticScheduler<unsigned int>::Parameters par;
    
    par.popSize      = par_.genetic.popSize;
    par.mutationRate = par_.genetic.mutationRate;
    par.seed         = rd();
    memPeak_         = 0;
    CartesianCommunicator::BroadcastWorld(0, &(par.seed), sizeof(par.seed));
    for (unsigned int i = 0; i < con.size(); ++i)
    {
        GeneticScheduler<unsigned int> scheduler(con[i], memPeak, par);
        
        gen = 0;
        do
        {
            LOG(Debug) << "Generation " << gen << ":" << std::endl;
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
        } while ((gen < par_.genetic.maxGen)
                 and (nCstPeak < par_.genetic.maxCstGen));
        auto &t = scheduler.getMinSchedule();
        if (scheduler.getMinValue() > memPeak_)
        {
            memPeak_ = scheduler.getMinValue();
        }
        for (unsigned int j = 0; j < t.size(); ++j)
        {
            program_.push_back(t[j]);
        }
    }
    scheduled_ = true;
}

void Application::saveSchedule(const std::string filename)
{
    TextWriter               writer(filename);
    std::vector<std::string> program;
    
    if (!scheduled_)
    {
        HADRON_ERROR("Computation not scheduled");
    }
    LOG(Message) << "Saving current schedule to '" << filename << "'..."
                 << std::endl;
    for (auto address: program_)
    {
        program.push_back(env().getModuleName(address));
    }
    write(writer, "schedule", program);
}

void Application::loadSchedule(const std::string filename)
{
    DEFINE_MEMPEAK;
    
    TextReader               reader(filename);
    std::vector<std::string> program;
    
    LOG(Message) << "Loading schedule from '" << filename << "'..."
                 << std::endl;
    read(reader, "schedule", program);
    program_.clear();
    for (auto &name: program)
    {
        program_.push_back(env().getModuleAddress(name));
    }
    scheduled_ = true;
    memPeak_   = memPeak(program_);
}

void Application::printSchedule(void)
{
    if (!scheduled_)
    {
        HADRON_ERROR("Computation not scheduled");
    }
    LOG(Message) << "Schedule (memory peak: " << MEM_MSG(memPeak_) << "):"
                 << std::endl;
    for (unsigned int i = 0; i < program_.size(); ++i)
    {
        LOG(Message) << std::setw(4) << i + 1 << ": "
                     << env().getModuleName(program_[i]) << std::endl;
    }
}

// loop on configurations //////////////////////////////////////////////////////
void Application::configLoop(void)
{
    auto range = par_.trajCounter;
    
    for (unsigned int t = range.start; t < range.end; t += range.step)
    {
        LOG(Message) << BIG_SEP << " Starting measurement for trajectory " << t
                     << " " << BIG_SEP << std::endl;
        env().setTrajectory(t);
        env().executeProgram(program_);
    }
    LOG(Message) << BIG_SEP << " End of measurement " << BIG_SEP << std::endl;
    env().freeAll();
}
