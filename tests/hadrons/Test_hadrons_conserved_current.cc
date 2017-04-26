/*******************************************************************************
 Grid physics library, www.github.com/paboyle/Grid
 
 Source file: tests/hadrons/Test_hadrons_conserved_current.cc
 
 Copyright (C) 2017
 
 Author: Andrew Lawson <andrew.lawson1991@gmail.com>
 
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

#include "Test_hadrons.hpp"

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
    Application  application;
    unsigned int nt = GridDefaultLatt()[Tp];
    double       mass = 0.04;

    // global parameters
    Application::GlobalPar globalPar;
    globalPar.trajCounter.start    = 1500;
    globalPar.trajCounter.end      = 1520;
    globalPar.trajCounter.step     = 20;
    globalPar.seed                 = "1 2 3 4";
    globalPar.genetic.maxGen       = 1000;
    globalPar.genetic.maxCstGen    = 200;
    globalPar.genetic.popSize      = 20;
    globalPar.genetic.mutationRate = .1;
    application.setPar(globalPar);

    // gauge field
    application.createModule<MGauge::Unit>("gauge");

    // action
    std::string actionName = "DWF";
    MAction::DWF::Par actionPar;
    actionPar.gauge = "gauge";
    actionPar.Ls    = 12;
    actionPar.M5    = 1.8;
    actionPar.mass  = mass;
    application.createModule<MAction::DWF>(actionName, actionPar);

    // solver
    std::string solverName = "CG";
    MSolver::RBPrecCG::Par solverPar;
    solverPar.action   = actionName;
    solverPar.residual = 1.0e-8;
    application.createModule<MSolver::RBPrecCG>(solverName,
                                                solverPar);

    // Conserved current sink contractions: use a single point propagator.
    std::string pointProp = "q_0";
    std::string pos = "0 0 0 0";
    std::string modName = "Ward Identity Test";
    MAKE_POINT_PROP(pos, pointProp, solverName);
    if (!(Environment::getInstance().hasModule(modName)))
    {
        MContraction::WardIdentity::Par wiPar;
        wiPar.q      = pointProp + "_5d";
        wiPar.q4d    = pointProp;
        wiPar.action = actionName;
        wiPar.mass   = mass;
        application.createModule<MContraction::WardIdentity>(modName, wiPar);
    }

    // Conserved current contractions with sequential insertion of vector
    // current.
    std::string q_x = "q_x";
    std::string q_y = "q_y";
    std::string q_z = "q_z";
    std::string q_t = "q_t";
    std::string mom = ZERO_MOM;
    modName         = "Sequential Ward Identity Test";
    MAKE_SEQUENTIAL_PROP(nt/2, pointProp, mom, q_x, solverName);
    MAKE_SEQUENTIAL_PROP(nt/2, pointProp, mom, q_y, solverName);
    MAKE_SEQUENTIAL_PROP(nt/2, pointProp, mom, q_z, solverName);
    MAKE_SEQUENTIAL_PROP(nt/2, pointProp, mom, q_t, solverName);
    if (!(Environment::getInstance().hasModule(modName)))
    {
        MContraction::WardIdentitySeq::Par wiPar;
        wiPar.q_x = q_x;
        wiPar.q_y = q_y;
        wiPar.q_z = q_z;
        wiPar.q_t = q_t;
        application.createModule<MContraction::WardIdentitySeq>(modName, wiPar);
    }

    // execution
    application.saveParameterFile("ConservedCurrentTest.xml");
    application.run();

    // epilogue
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();

    return EXIT_SUCCESS;
}