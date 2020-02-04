/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Tests/Hadrons/Test_hadrons_spectrum.cc

Copyright (C) 2015-2018

 Author: Antonin Portelli <antonin.portelli@me.com>
 Author: Felix Erben <felix.erben@ed.ac.uk>

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

#include <Hadrons/Application.hpp>
#include <Hadrons/Modules.hpp>

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
    Application              application;
    std::vector<std::string> flavour = {"l", "s", "c"};
    std::vector<double>      mass    = {.01, .04, .2 };
    
    // global parameters
    Application::GlobalPar globalPar;
    globalPar.trajCounter.start = 1500;
    globalPar.trajCounter.end   = 1520;
    globalPar.trajCounter.step  = 20;
    globalPar.runId             = "test";
    application.setPar(globalPar);
    // gauge field
    application.createModule<MGauge::Unit>("gauge");
    // sources
    MSource::Point::Par ptPar;
    ptPar.position = "0 0 0 0";
    application.createModule<MSource::Point>("pt_0", ptPar);
    ptPar.position = "0 0 0 4";
    application.createModule<MSource::Point>("pt_4", ptPar);
    // sink
    MSink::Point::Par sinkPar;
    sinkPar.mom = "0 0 0";
    application.createModule<MSink::ScalarPoint>("sink", sinkPar);
    application.createModule<MSink::Point>("sink_spec", sinkPar);
    
    // set fermion boundary conditions to be periodic space, antiperiodic time.
    std::string boundary = "1 1 1 -1";
    std::string twist = "0. 0. 0. 0.";

    for (unsigned int i = 0; i < flavour.size(); ++i)
    {
        // actions
        MAction::DWF::Par actionPar;
        actionPar.gauge = "gauge";
        actionPar.Ls    = 12;
        actionPar.M5    = 1.8;
        actionPar.mass  = mass[i];
        actionPar.boundary = boundary;
        actionPar.twist = twist;
        application.createModule<MAction::DWF>("DWF_" + flavour[i], actionPar);
        
        // solvers
        MSolver::RBPrecCG::Par solverPar;
        solverPar.action       = "DWF_" + flavour[i];
        solverPar.residual     = 1.0e-8;
        solverPar.maxIteration = 10000;
        application.createModule<MSolver::RBPrecCG>("CG_" + flavour[i],
                                                    solverPar);
        
    }

    // propagators
    MFermion::GaugeProp::Par quarkPar;
    quarkPar.solver = "CG_l";
    quarkPar.source = "pt_0";
    application.createModule<MFermion::GaugeProp>("Qpt_l_0", quarkPar);
    quarkPar.source = "pt_4";
    application.createModule<MFermion::GaugeProp>("Qpt_l_4", quarkPar);
    quarkPar.solver = "CG_s";
    quarkPar.source = "pt_0";
    application.createModule<MFermion::GaugeProp>("Qpt_s_0", quarkPar);
    //This should be a loop - how do I make this?
    quarkPar.solver = "CG_c";
    quarkPar.source = "pt_0";
    application.createModule<MFermion::GaugeProp>("Qpt_c_loop", quarkPar);
    quarkPar.solver = "CG_l";
    quarkPar.source = "pt_0";
    application.createModule<MFermion::GaugeProp>("Qpt_l_loop", quarkPar);

    MSink::Smear::Par smearPar;
    smearPar.q="Qpt_l_0";
    smearPar.sink = "sink_spec";
    application.createModule<MSink::Smear>("Qpt_d_spec",smearPar);
    smearPar.q="Qpt_s_0";
    smearPar.sink = "sink_spec";
    application.createModule<MSink::Smear>("Qpt_s_spec",smearPar);

    MContraction::XiToSigmaEye::Par EyePar;
    EyePar.output  = "XiToSigma/Eye_u";
    EyePar.qqLoop = "Qpt_l_loop";
    EyePar.qdSpec = "Qpt_d_spec";
    EyePar.qdTf   = "Qpt_l_4";
    EyePar.qsSpec = "Qpt_s_spec";
    EyePar.qsTi   = "Qpt_s_0";
    EyePar.tf    = 4;
    EyePar.sink    = "sink";
    application.createModule<MContraction::XiToSigmaEye>("XiToSigmaEye_u", EyePar);
    EyePar.output  = "XiToSigma/Eye_c";
    EyePar.qqLoop = "Qpt_c_loop";
    application.createModule<MContraction::XiToSigmaEye>("XiToSigmaEye_c", EyePar);

    // execution
    application.saveParameterFile("xts.xml");
    application.run();
    
    // epilogue
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();
    
    return EXIT_SUCCESS;
}
