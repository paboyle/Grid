/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Tests/Hadrons/Test_hadrons_seq_gamma.cc

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
using namespace QCD;
using namespace Hadrons;

/*******************************************************************************
 * Consistency test for sequential gamma insertion.
 ******************************************************************************/

int main(int argc, char *argv[])
{
    // initialization //////////////////////////////////////////////////////////
    HADRONS_DEFAULT_INIT;

    // run setup ///////////////////////////////////////////////////////////////
    Application  application;
    unsigned int nt   = GridDefaultLatt()[Tp];
    unsigned int tS   = nt / 2;
    unsigned int Ls   = 12;
    double       mass = 0.04;
    double       M5   = 1.8;

    // global parameters
    HADRONS_DEFAULT_GLOBALS(application);

    // gauge field
    std::string gaugeField = "gauge";
    application.createModule<MGauge::Unit>(gaugeField);

    // action
    std::string actionName = "DWF";
    makeDWFAction(application, actionName, gaugeField, mass, M5, Ls);

    // solver
    std::string solverName = "CG";
    makeRBPrecCGSolver(application, solverName, actionName);

    // test sequential propagator, with g5 insertion.
    Gamma::Algebra g = Gamma::Algebra::Gamma5;
    std::string pointProp = "q_0";
    std::string point5d   = LABEL_5D(pointProp);
    std::string origin    = "0 0 0 0";
    MAKE_POINT_PROP(origin, pointProp, solverName);

    std::string seqProp = ADD_INDEX(pointProp + "_seqg5", tS);
    std::string seqSrc  = seqProp + "_src";
    MAKE_SEQUENTIAL_PROP(tS, pointProp, ZERO_MOM, seqProp, solverName, g);

    std::string modName = "Test g5 sequential insertion";
    makeSeqGamComparison(application, modName, pointProp, seqProp, origin, g, tS);

    // execution
    application.saveParameterFile("SeqGamma5Test.xml");
    application.run();

    // epilogue
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();

    return EXIT_SUCCESS;
}
