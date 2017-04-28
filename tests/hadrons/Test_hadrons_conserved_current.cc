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
    unsigned int Ls = 12;

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
    actionPar.Ls    = Ls;
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
    makeWITest(application, modName, pointProp, actionName, mass, Ls);

    // Conserved current contractions with sequential insertion of vector/axial
    // current.
    std::string mom      = ZERO_MOM;
    unsigned int t_J     = nt/2;
    std::string seqPropA = ADD_INDEX(pointProp + "_seq_A", t_J);
    std::string seqPropV = ADD_INDEX(pointProp + "_seq_V", t_J);
    std::string seqSrcA  = seqPropA + "_src";
    std::string seqSrcV  = seqPropV + "_src";
    std::string point5d  = LABEL_5D(pointProp);
    makeConservedSequentialSource(application, seqSrcA, point5d, 
                                  actionName, t_J, Current::Axial, Tp, mom);
    makePropagator(application, seqPropA, seqSrcA, solverName);
    makeConservedSequentialSource(application, seqSrcV, point5d, 
                                  actionName, t_J, Current::Vector, Tp, mom);
    makePropagator(application, seqPropV, seqSrcV, solverName);

    std::string modNameA = "Axial Sequential Test";
    std::string modNameV = "Vector Sequential Test";
    makeSeqTest(application, modNameA, pointProp, seqPropA, 
                actionName, pos, t_J, Tp, Current::Axial, Ls);
    makeSeqTest(application, modNameV, pointProp, seqPropV, 
                actionName, pos, t_J, Tp, Current::Vector, Ls);

    // execution
    application.saveParameterFile("ConservedCurrentTest.xml");
    application.run();

    // epilogue
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();

    return EXIT_SUCCESS;
}