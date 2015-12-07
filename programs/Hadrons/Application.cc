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
}

// parse parameter file ////////////////////////////////////////////////////////
void Application::parseParameterFile(void)
{
    XmlReader reader(parameterFileName_);
    
    LOG(Message) << "Reading '" << parameterFileName_ << "'..." << endl;
    read(reader, "parameters", parameters_);
}

// schedule computation ////////////////////////////////////////////////////////
void Application::schedule(void)
{
    Graph<string> moduleGraph;
    
    LOG(Message) << "Scheduling computation..." << endl;
    for (auto &m: parameters_.modules)
    {
        for (auto &p: m.in)
        {
            moduleGraph.addEdge(p, m.name);
        }
    }
    
    vector<Graph<string>> con = moduleGraph.getConnectedComponents();
    
    LOG(Message) << "Program:" << endl;
    LOG(Message) << "  #segments: " << con.size() << endl;
    for (unsigned int i = 0; i < con.size(); ++i)
    {
        vector<vector<string>> t = con[i].allTopoSort();
        auto                   m = makeDependencyMatrix(t);
        
        for (auto &v: t[0])
        {
            cout << v << " ";
        }
        cout << endl;
        for (auto &v1: t[0])
        {
            for (auto &v2: t[0])
            {
                cout << m[v1][v2] << " ";
            }
            cout << endl;
        }
        
        LOG(Message) << "  segment " << i << ":" << endl;
    }
}
