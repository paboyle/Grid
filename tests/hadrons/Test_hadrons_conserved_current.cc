/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Tests/Hadrons/Test_hadrons_conserved_current.cc

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

inline void setupSeqCurrTests(Application &application, std::string modStem, 
                              std::string &pointProp, std::string &seqStem,
                              std::string &actionName, std::string &solverName,
                              std::string &origin, Current curr, 
                              unsigned int t_J, unsigned int mu,
                              unsigned int Ls = 1)
{
    std::string modName = ADD_INDEX(modStem, mu);
    std::string seqProp = ADD_INDEX(seqStem, mu);
    std::string seqSrc  = seqProp + "_src";

    // 5D actions require 5D propagator as input for conserved current
    // insertions.
    std::string propIn;
    if (Ls > 1) 
    {
        propIn = LABEL_5D(pointProp); 
    }
    else
    { 
        propIn = pointProp; 
    }

    makeConservedSequentialSource(application, seqSrc, propIn,
                                  actionName, t_J, curr, mu);
    makePropagator(application, seqProp, seqSrc, solverName);
    makeSeqCurrComparison(application, modName, propIn, seqProp, 
                          actionName, origin, t_J, mu, curr);
}

inline void setupWardIdentityTests(Application &application,
                                   std::string &actionName, 
                                   double mass,
                                   unsigned int Ls = 1, 
                                   bool perform_axial_tests = false)
{
    // solver
    std::string solverName = actionName + "_CG";
    makeRBPrecCGSolver(application, solverName, actionName);

    unsigned int nt = GridDefaultLatt()[Tp];
    unsigned int t_J = nt/2;

    /***************************************************************************
     * Conserved current sink contractions: use a single point propagator for
     * the Ward Identity test.
     **************************************************************************/
    std::string pointProp = actionName + "_q_0";
    std::string origin    = "0 0 0 0";
    std::string modName   = actionName + " Ward Identity Test";
    MAKE_POINT_PROP(origin, pointProp, solverName);
    makeWITest(application, modName, pointProp, actionName, mass, Ls,
               perform_axial_tests);

    /***************************************************************************
     * Conserved current tests with sequential insertion of vector/axial 
     * current. If above Ward Identity passes, sufficient to test sequential
     * insertion of conserved current agrees with contracted version.
     **************************************************************************/
    // Compare sequential insertion to contraction. Should be enough to perform 
    // for time and one space component.
    std::string seqStem = ADD_INDEX(pointProp + "seq_V", t_J);
    std::string modStem = actionName + " Vector Sequential Test mu";
    setupSeqCurrTests(application, modStem, pointProp, seqStem, actionName, 
                      solverName, origin, Current::Vector, t_J, Tp, Ls);
    setupSeqCurrTests(application, modStem, pointProp, seqStem, actionName, 
                      solverName, origin, Current::Vector, t_J, Xp, Ls);

    // Perform axial tests only if partially-conserved axial current exists for
    // the action.
    if (perform_axial_tests)
    {
        seqStem = ADD_INDEX(pointProp + "seq_A", t_J);
        modStem = actionName + " Axial Sequential Test mu";
        setupSeqCurrTests(application, modStem, pointProp, seqStem, actionName, 
                          solverName, origin, Current::Axial, t_J, Tp, Ls);
        setupSeqCurrTests(application, modStem, pointProp, seqStem, actionName, 
                          solverName, origin, Current::Axial, t_J, Xp, Ls);
    }
}

int main(int argc, char *argv[])
{
    // initialization //////////////////////////////////////////////////////////
    HADRONS_DEFAULT_INIT;

    // run setup ///////////////////////////////////////////////////////////////
    Application  application;
    double       mass = 0.04;
    double       M5 = 1.8;
    unsigned int Ls = 12;

    // global parameters
    HADRONS_DEFAULT_GLOBALS(application);

    // gauge field
    std::string gaugeField = "gauge";
    application.createModule<MGauge::Unit>(gaugeField);

    // Setup each action and the conserved current tests relevant to it.
    std::string actionName = "DWF";
    makeDWFAction(application, actionName, gaugeField, mass, M5, Ls);
    setupWardIdentityTests(application, actionName, mass, Ls, true);

    actionName = "Wilson";
    makeWilsonAction(application, actionName, gaugeField, mass);
    setupWardIdentityTests(application, actionName, mass);

    // execution
    application.saveParameterFile("ConservedCurrentTest.xml");
    application.run();

    // epilogue
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();

    return EXIT_SUCCESS;
}
