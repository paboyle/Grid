/*
 * Application.cc, part of Grid
 *
 * Copyright (C) 2015 Antonin Portelli
 *
 * Grid is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Grid is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Grid.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <Hadrons/Application.hpp>
#include <Hadrons/Graph.hpp>

using namespace std;
using namespace Grid;
using namespace Hadrons;

/******************************************************************************
 *                       Application implementation                           *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
Application::Application(int argc, char *argv[])
: env_(Environment::getInstance())
, modFactory_(ModuleFactory::getInstance())
{
    if (argc < 2)
    {
        cerr << "usage: " << argv[0] << " <parameter file> [Grid options]";
        cerr << endl;
        exit(EXIT_FAILURE);
    }
    parameterFileName_ = argv[1];
    Grid_init(&argc, &argv);
    HadronsLogError.Active(GridLogError.isActive());
    HadronsLogWarning.Active(GridLogWarning.isActive());
    HadronsLogMessage.Active(GridLogMessage.isActive());
    HadronsLogDebug.Active(GridLogDebug.isActive());
    LOG(Message) << "Grid initialized" << endl;
    LOG(Message) << "Modules available:" << endl;
    auto list = modFactory_.getModuleList();
    for (auto &m: list)
    {
        LOG(Message) << "  " << m << endl;
    }
}

// destructor //////////////////////////////////////////////////////////////////
Application::~Application(void)
{
    LOG(Message) << "Grid is finalizing now" << endl;
    Grid_finalize();
}

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
    
    LOG(Message) << "Reading '" << parameterFileName_ << "'..." << endl;
    read(reader, "parameters", par_);
    push(reader, "modules");
    push(reader, "module");
    do
    {
        read(reader, "id", id);
        module_[id.name] = modFactory_.create(id.type, id.name);
        module_[id.name]->parseParameters(reader, "options");
        vector<string> output = module_[id.name]->getOutput();
        for (auto &n: output)
        {
            associatedModule_[n] = id.name;
        }
    } while (reader.nextElement("module"));
    pop(reader);
    pop(reader);
}

// schedule computation ////////////////////////////////////////////////////////
void Application::schedule(void)
{
    Graph<string> moduleGraph;
    
    LOG(Message) << "Scheduling computation..." << endl;
    
    // create dependency graph
    for (auto &m: module_)
    {
        vector<string> input = m.second->getInput();
        for (auto &n: input)
        {
            try
            {
                moduleGraph.addEdge(associatedModule_.at(n), m.first);
            }
            catch (out_of_range &)
            {
                HADRON_ERROR("unknown object '" + n + "'");
            }
        }
    }
    
    // topological sort
    map<string, map<string, bool>> m;
    unsigned int                   k = 0;
    
    vector<Graph<string>> con = moduleGraph.getConnectedComponents();
    LOG(Message) << "Program:" << endl;
    for (unsigned int i = 0; i < con.size(); ++i)
    {
        vector<vector<string>> t = con[i].allTopoSort();

        m = makeDependencyMatrix(t);
        for (unsigned int j = 0; j < t[0].size(); ++j)
        {
            program_.push_back(t[0][j]);
            LOG(Message) << setw(4) << right << k << ": "
                         << program_[k] << endl;
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
        LOG(Message) << "Starting measurement for trajectory " << t << endl;
        execute();
    }
}

void Application::execute(void)
{
    for (unsigned int i = 0; i < program_.size(); ++i)
    {
        LOG(Message) << "Measurement step (" << i+1 << "/" << program_.size()
                     << ")" << endl;
        (*module_[program_[i]])(env_);
    }
}
