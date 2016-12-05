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
: env_(Environment::getInstance())
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

// access //////////////////////////////////////////////////////////////////////
void Application::setPar(const Application::GlobalPar &par)
{
    par_ = par;
    env_.setSeed(strToVec<int>(par_.seed));
}

// execute /////////////////////////////////////////////////////////////////////
void Application::run(void)
{
    if (!parameterFileName_.empty())
    {
        parseParameterFile();
    }
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
    GlobalPar par;
    ObjectId  id;
    
    LOG(Message) << "Reading '" << parameterFileName_ << "'..." << std::endl;
    read(reader, "parameters", par);
    setPar(par);
    push(reader, "modules");
    push(reader, "module");
    do
    {
        read(reader, "id", id);
        env_.createModule(id.name, id.type, reader);
    } while (reader.nextElement("module"));
    pop(reader);
    pop(reader);
    env_.setSeed(strToVec<int>(par_.seed));
    env_.printContent();
}

// schedule computation ////////////////////////////////////////////////////////
#define MEM_MSG(size)\
sizeString((size)*locVol_) << " (" << sizeString(size)  << "/site)"

void Application::schedule(void)
{
    // memory peak function
    auto memPeak = [this](const std::vector<unsigned int> &program)
    {
        unsigned int memPeak;
        bool         msg;
        
        msg = HadronsLogMessage.isActive();
        HadronsLogMessage.Active(false);
        env_.dryRun(true);
        memPeak = env_.executeProgram(program);
        env_.dryRun(false);
        env_.freeAll();
        HadronsLogMessage.Active(true);
        
        return memPeak;
    };
    
    // constrained topological sort using a genetic algorithm
    LOG(Message) << "Scheduling computation..." << std::endl;
    constexpr unsigned int maxGen = 200, maxCstGen = 50;
    unsigned int           k = 0, gen, prevPeak, nCstPeak = 0;
    auto                   graph = env_.makeModuleGraph();
    auto                   con = graph.getConnectedComponents();
    std::random_device     rd;
    GeneticScheduler<unsigned int>::Parameters par;

    par.popSize      = 10;
    par.mutationRate = .1;
    par.seed         = rd();
    CartesianCommunicator::BroadcastWorld(0, &(par.seed), sizeof(par.seed));
    for (unsigned int i = 0; i < con.size(); ++i)
    {
        GeneticScheduler<unsigned int> scheduler(con[i], memPeak, par);
        
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
            LOG(Debug) << "generation " << gen << ":\n" << scheduler;
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
            LOG(Message) << std::setw(4) << k + 1 << ": "
                         << env_.getModuleName(program_[k]) << std::endl;
            k++;
        }
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
        env_.setTrajectory(t);
        env_.executeProgram(program_);
        env_.freeAll();
    }
    LOG(Message) << BIG_SEP << " End of measurement " << BIG_SEP << std::endl;
}
