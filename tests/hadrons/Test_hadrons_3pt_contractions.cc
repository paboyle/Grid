/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Tests/Hadrons/Test_hadrons_3pt_contractions.cc

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
    Application  application;
    double       mass = 0.04;
    double       M5   = 1.8;
    unsigned int Ls   = 12;
    unsigned int nt   = GridDefaultLatt()[Tp];
    unsigned int t_i  = 0;
    unsigned int t_f  = nt / 2;
    std::string  mom  = "1. 0. 0. 0.";

    // global parameters
    HADRONS_DEFAULT_GLOBALS(application);

    // gauge field
    std::string gaugeField = "gauge";
    application.createModule<MGauge::Unit>(gaugeField);

    // Action & solver setup.
    std::string action = "DWF";
    std::string solver = "CG";
    makeDWFAction(application, action, gaugeField, mass, M5, Ls);
    makeRBPrecCGSolver(application, solver, action);

    /***************************************************************************
     * Weak Contraction test: Non-Eye class.
     **************************************************************************/
    // Make wall source propagators for each leg of 4-quark vertex.
    std::string q_i_0 = "q_i_0";
    std::string q_i_p = "q_i_p";
    std::string q_f_0 = "q_f_0";
    std::string q_f_p = "q_f_p";
    MAKE_WALL_PROP(t_i, q_i_0, solver);
    MAKE_WALL_PROP(t_f, q_f_0, solver);
    MAKE_3MOM_WALL_PROP(t_i, mom, q_i_p, solver);
    MAKE_3MOM_WALL_PROP(t_f, mom, q_f_p, solver);

    // Perform contractions, zero and non-zero momentum.
    std::string HW_CW_0 = LABEL_3PT("HW_CW_0", t_i, t_f);
    std::string HW_CW_p = LABEL_3PT("HW_CW_p", t_i, t_f);
    weakContractionNonEye(application, 3, q_i_0, q_i_0, q_f_0, q_f_0, HW_CW_0);
    weakContractionNonEye(application, 3, q_i_0, q_i_p, q_f_p, q_f_0, HW_CW_p);

    /***************************************************************************
     * Weak Contraction test: Eye-class.
     **************************************************************************/
    // Create random propagator for loop.
    std::string eta = "noise_source";
    makeNoiseSource(application, eta, 0, nt - 1);
    std::string loopProp = "loop";
    std::string loopRes  = loopProp + "_res";
    makePropagator(application, loopRes, eta, solver);
    makeLoop(application, loopProp, eta, loopRes);

    // Wall sink smear the propagator directly connecting the source & sink.
    // (i.e. make point sink but smear before the contraction)
    std::string wallSink = "wall_sink";
    std::string qWall    = "q_wall";
    makePointSink(application, wallSink);
    sinkSmear(application, wallSink, q_i_0, qWall);

    // Perform contractions, zero and non-zero momentum.
    std::string HW_SE_0 = LABEL_3PT("HW_SE_0", t_i, t_f);
    std::string HW_SE_p = LABEL_3PT("HW_SE_p", t_i, t_f);
    weakContractionEye(application, 3, qWall, q_i_0, q_f_p, loopProp, HW_SE_0, t_f);
    weakContractionEye(application, 3, qWall, q_i_p, q_f_p, loopProp, HW_SE_p, t_f);

    /***************************************************************************
     * Gamma insertion test.
     **************************************************************************/
    Gamma::Algebra gamma = Gamma::Algebra::GammaT;
    std::string sd_0 = LABEL_3PT("sd_0", t_i, t_f);
    std::string sd_p = LABEL_3PT("sd_p", t_i, t_f);
    gamma3ptContraction(application, 3, qWall, q_i_0, q_f_0, sd_0, t_f, gamma);
    gamma3ptContraction(application, 3, qWall, q_i_p, q_f_p, sd_p, t_f, gamma);

    // execution
    application.saveParameterFile("ContractionTest3pt.xml");
    application.run();

    // epilogue
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();

    return EXIT_SUCCESS;
}
