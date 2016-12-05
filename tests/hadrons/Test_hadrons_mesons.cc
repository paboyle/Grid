/*******************************************************************************
 Grid physics library, www.github.com/paboyle/Grid
 
 Source file: programs/Hadrons/Hadrons.cc
 
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

using namespace Grid;
using namespace Hadrons;

int main(int argc, char *argv[])
{
    // initialization //////////////////////////////////////////////////////////
    Grid_init(&argc, &argv);
    HadronsLogError.Active(GridLogError.isActive());
    HadronsLogWarning.Active(GridLogWarning.isActive());
    HadronsLogMessage.Active(GridLogMessage.isActive());
    HadronsLogIterative.Active(GridLogIterative.isActive());
    HadronsLogDebug.Active(GridLogDebug.isActive());
    LOG(Message) << "Grid initialized" << std::endl;
    
    // run setup ///////////////////////////////////////////////////////////////
    Application                   application;

    std::vector<std::string> flavour = {"l", "s", "c1", "c2", "c3"};
    std::vector<double>      mass    = {.01, .04, .2  , .25 , .3  };
    
    // global parameters
    Application::GlobalPar globalPar;
    globalPar.trajCounter.start = 1500;
    globalPar.trajCounter.end   = 1520;
    globalPar.trajCounter.step  = 20;
    globalPar.seed              = "1 2 3 4";
    application.setPar(globalPar);
    // gauge field
    application.createModule<MGauge::Unit>("gauge");
    // sources
    MSource::Z2::Par z2Par;
    z2Par.tA = 0;
    z2Par.tB = 0;
    application.createModule<MSource::Z2>("z2", z2Par);
    MSource::Point::Par ptPar;
    ptPar.position = "0 0 0 0";
    application.createModule<MSource::Point>("pt", ptPar);
    for (unsigned int i = 0; i < flavour.size(); ++i)
    {
        // actions
        MAction::DWF::Par actionPar;
        actionPar.gauge = "gauge";
        actionPar.Ls    = 12;
        actionPar.M5    = 1.8;
        actionPar.mass  = mass[i];
        application.createModule<MAction::DWF>("DWF_" + flavour[i], actionPar);
        
        // solvers
        MSolver::RBPrecCG::Par solverPar;
        solverPar.action   = "DWF_" + flavour[i];
        solverPar.residual = 1.0e-8;
        application.createModule<MSolver::RBPrecCG>("CG_" + flavour[i],
                                                    solverPar);
        
        // propagators
        MQuark::Par quarkPar;
        quarkPar.solver = "CG_" + flavour[i];
        quarkPar.source = (flavour[i][0] == 'c') ? "z2" : "pt";
        application.createModule<MQuark>("Q_" + flavour[i], quarkPar);
    }
    for (unsigned int i = 0;     i < flavour.size(); ++i)
    for (unsigned int j = i + 1; j < flavour.size(); ++j)
    {
        MContraction::Meson::Par mesPar;
        
        mesPar.output = "mesons/" + flavour[i] + flavour[j];
        mesPar.q1     = "Q_" + flavour[i];
        mesPar.q2     = "Q_" + flavour[j];
        application.createModule<MContraction::Meson>("meson_"
                                                      + flavour[i] + flavour[j],
                                                      mesPar);
    }
    
    // execution
    application.run();
    
    // epilogue
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();
    
    return EXIT_SUCCESS;
}
