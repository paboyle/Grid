/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Tests/Hadrons/Test_hadrons_meson_conserved_3pt.cc

Copyright (C) 2015-2018

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

See the full license in the file "LICENSE" in the top level distribution directory
*************************************************************************************/
/*  END LEGAL */

#include "Test_hadrons.hpp"

using namespace Grid;
using namespace Hadrons;

int main(int argc, char *argv[])
{
    // initialization //////////////////////////////////////////////////////////
    HADRONS_DEFAULT_INIT;
    
    // run setup ///////////////////////////////////////////////////////////////
    Application              application;
    
    // actions parameters
    double        mass = 0.04;
    unsigned int  Ls = 16;
    double        M5 = 1.8;
    
    // kinematics
    unsigned int  nt    = GridDefaultLatt()[Tp];
    unsigned int  tSrc  = 0;
    unsigned int  tJ    = nt / 4;
    std::string   kmom  = "0. 0. 0. 0.";
    std::string   pmom  = "1. 0. 0. 0.";

    // Global parameters.
    HADRONS_DEFAULT_GLOBALS(application);

    // Unit gauge field.
    std::string gaugeField = "Unit gauge";
    application.createModule<MGauge::Unit>(gaugeField);

    // DWF action
    std::string actionName = "DWF";
    makeDWFAction(application, actionName, gaugeField, mass, M5, Ls);

    // Solver
    std::string solver = "CG";
    makeRBPrecCGSolver(application, solver, actionName);
    
    // main test body //////////////////////////////////////////////////////////
    // Point sink modules.
    std::string sink_0 = "sink_0";
    std::string sink_p = "sink_p";
    MSink::Point::Par sinkPar;
    sinkPar.mom = kmom;
    application.createModule<MSink::ScalarPoint>(sink_0, sinkPar);
    sinkPar.mom = pmom;
    application.createModule<MSink::ScalarPoint>(sink_p, sinkPar);

    // 2pt pion contraction, zero momentum.
    std::string q_0 = "Q_0";
    MAKE_WALL_PROP(tSrc, q_0, solver);
    std::string modName = INIT_INDEX("2pt_pion_WP", tSrc);
    std::string output  = "2pt/pion_WP_0";
    mesonContraction(application, modName, output, q_0, q_0, sink_0);

    // 2pt pion contraction, with momentum p.
    std::string q_p = "Q_p";
    MAKE_3MOM_WALL_PROP(tSrc, pmom, q_p, solver);
    modName = INIT_INDEX("2pt_pion_WP_p", tSrc);
    output  = "2pt/pion_WP_p";
    mesonContraction(application, modName, output, q_0, q_p, sink_p);

    // 3pt pion(0) -> pion(p), with sequentially inserted vector current in 
    // time direction.
    std::string qSeq    = q_0 + INIT_INDEX("_seq_Vc3", tJ);
    std::string q5d     = LABEL_5D(q_0); // Need 5D prop for DWF conserved current.
    std::string srcName = qSeq + "_src";
    modName = LABEL_3PT("3pt_pion_Vc3", tSrc, tJ);
    output  = "3pt/pion_Vc3_p";
    makeConservedSequentialSource(application, srcName, q5d, actionName,
                                  tJ, Current::Vector, Tp, pmom);
    makePropagator(application, qSeq, srcName, solver);
    mesonContraction(application, modName, output, q_0, qSeq, sink_p);

    std::string par_file_name = "conserved_3pt.xml";
    application.saveParameterFile(par_file_name);
    application.run();

    // epilogue
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();

    return EXIT_SUCCESS;
}
    
    
