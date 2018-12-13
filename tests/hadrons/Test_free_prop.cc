/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Tests/Hadrons/Test_free_prop.cc

Copyright (C) 2015-2018

 Author: Antonin Portelli <antonin.portelli@me.com>
 Author: Vera Guelpers    <v.m.guelpers@soton.ac.uk>

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
    std::vector<std::string> flavour = {"h"}; //{"l", "s", "c1", "c2", "c3"};
    std::vector<double>      mass    = {.2}; //{.01, .04, .2  , .25 , .3  };
    std::vector<std::string> lepton_flavour    = {"mu"};
    std::vector<double>      lepton_mass    = {.2};

    unsigned int  nt    = GridDefaultLatt()[Tp];
    
    // global parameters
    Application::GlobalPar globalPar;
    globalPar.trajCounter.start = 1500;
    globalPar.trajCounter.end   = 1520;
    globalPar.trajCounter.step  = 20;
    globalPar.runId             = "test";
    application.setPar(globalPar);
    // gauge field
    application.createModule<MGauge::Unit>("gauge");
    // unit gauge field for lepton 
    application.createModule<MGauge::Unit>("free_gauge");
    // pt source
    MSource::Point::Par ptPar;
    ptPar.position = "0 0 0 0";
    application.createModule<MSource::Point>("pt", ptPar);
    // sink
    MSink::Point::Par sinkPar;
    sinkPar.mom = "0 0 0";
    application.createModule<MSink::ScalarPoint>("sink", sinkPar);
    
    // set fermion boundary conditions to be periodic space, antiperiodic time.
    std::string boundary = "1 1 1 -1";


    //Propagators from FFT and Feynman rules
    for (unsigned int i = 0; i < lepton_mass.size(); ++i)
    {
        //DWF actions
        MAction::DWF::Par actionPar_lep;
        actionPar_lep.gauge = "free_gauge";
        actionPar_lep.Ls    = 8;
        actionPar_lep.M5    = 1.8;
        actionPar_lep.mass  = lepton_mass[i];
        actionPar_lep.boundary = boundary;
        application.createModule<MAction::DWF>("free_DWF_" + lepton_flavour[i], actionPar_lep);

        //DWF free propagators
        MFermion::FreeProp::Par freePar;
        freePar.source = "pt";
	freePar.action = "free_DWF_" + lepton_flavour[i];
	freePar.twist = "0 0 0 0.5";
        freePar.mass = lepton_mass[i];
        application.createModule<MFermion::FreeProp>("Lpt_" + lepton_flavour[i],
							 freePar);

        //Wilson actions
        MAction::Wilson::Par actionPar_lep_W;
        actionPar_lep_W.gauge = "free_gauge";
        actionPar_lep_W.mass  = lepton_mass[i];
        actionPar_lep_W.boundary = boundary;
        application.createModule<MAction::Wilson>("free_W_" + lepton_flavour[i], actionPar_lep_W);

        //Wilson free propagators
        MFermion::FreeProp::Par freePar_W;
        freePar_W.source = "pt";
	freePar_W.action = "free_W_" + lepton_flavour[i];
	freePar_W.twist = "0 0 0 0.5";
        freePar_W.mass = lepton_mass[i];
        application.createModule<MFermion::FreeProp>("W_Lpt_" + lepton_flavour[i],
							 freePar_W);


    }



    //Propagators from inversion
    for (unsigned int i = 0; i < flavour.size(); ++i)
    {
        //DWF actions
        MAction::DWF::Par actionPar;
        actionPar.gauge = "gauge";
        actionPar.Ls    = 8;
        actionPar.M5    = 1.8;
        actionPar.mass  = mass[i];
        actionPar.boundary = boundary;
        application.createModule<MAction::DWF>("DWF_" + flavour[i], actionPar);
        
        // solvers
        MSolver::RBPrecCG::Par solverPar;
        solverPar.action       = "DWF_" + flavour[i];
        solverPar.residual     = 1.0e-8;
        solverPar.maxIteration = 10000;
        application.createModule<MSolver::RBPrecCG>("CG_" + flavour[i],
                                                    solverPar);
        
        //DWF propagators
        MFermion::GaugeProp::Par quarkPar;
        quarkPar.solver = "CG_" + flavour[i];
        quarkPar.source = "pt";
        application.createModule<MFermion::GaugeProp>("Qpt_" + flavour[i],
							 quarkPar);



        //Wilson actions
        MAction::Wilson::Par actionPar_W;
        actionPar_W.gauge = "gauge";
        actionPar_W.mass  = mass[i];
        actionPar_W.boundary = boundary;
        application.createModule<MAction::Wilson>("W_" + flavour[i], actionPar_W);

        
        // solvers
        MSolver::RBPrecCG::Par solverPar_W;
        solverPar_W.action       = "W_" + flavour[i];
        solverPar_W.residual     = 1.0e-8;
        solverPar_W.maxIteration = 10000;
        application.createModule<MSolver::RBPrecCG>("W_CG_" + flavour[i],
                                                    solverPar_W);
        
        //Wilson propagators
        MFermion::GaugeProp::Par quarkPar_W;
        quarkPar_W.solver = "W_CG_" + flavour[i];
        quarkPar_W.source = "pt";
        application.createModule<MFermion::GaugeProp>("W_Qpt_" + flavour[i],
							 quarkPar_W);

    }


    //2pt contraction for Propagators from FFT and Feynman rules
    for (unsigned int i = 0; i < lepton_flavour.size(); ++i)
    for (unsigned int j = i; j < lepton_flavour.size(); ++j)
    {
        //2pt function contraction DWF
        MContraction::Meson::Par freemesPar;
        freemesPar.output  = "2pt_free/DWF_L_pt_" + lepton_flavour[i] + lepton_flavour[j];
        freemesPar.q1      = "Lpt_" + lepton_flavour[i];
        freemesPar.q2      = "Lpt_" + lepton_flavour[j];
        freemesPar.gammas  = "(Gamma5 Gamma5)";
        freemesPar.sink    = "sink";
        application.createModule<MContraction::Meson>("meson_L_pt_"
                                                      + lepton_flavour[i] + lepton_flavour[j],
                                                      freemesPar);

        //2pt function contraction Wilson
        MContraction::Meson::Par freemesPar_W;
        freemesPar_W.output  = "2pt_free/W_L_pt_" + lepton_flavour[i] + lepton_flavour[j];
        freemesPar_W.q1      = "W_Lpt_" + lepton_flavour[i];
        freemesPar_W.q2      = "W_Lpt_" + lepton_flavour[j];
        freemesPar_W.gammas  = "(Gamma5 Gamma5)";
        freemesPar_W.sink    = "sink";
        application.createModule<MContraction::Meson>("W_meson_L_pt_"
                                                      + lepton_flavour[i] + lepton_flavour[j],
                                                      freemesPar_W);

    }

    //2pt contraction for Propagators from inverion
    for (unsigned int i = 0; i < flavour.size(); ++i)
    for (unsigned int j = i; j < flavour.size(); ++j)
    {
        //2pt function contraction DWF
        MContraction::Meson::Par mesPar;
        mesPar.output  = "2pt_free/DWF_pt_" + flavour[i] + flavour[j];
        mesPar.q1      = "Qpt_" + flavour[i];
        mesPar.q2      = "Qpt_" + flavour[j];
        mesPar.gammas  = "(Gamma5 Gamma5)";
        mesPar.sink    = "sink";
        application.createModule<MContraction::Meson>("meson_pt_"
                                                      + flavour[i] + flavour[j],
                                                      mesPar);


        //2pt function contraction Wilson
        MContraction::Meson::Par mesPar_W;
        mesPar_W.output  = "2pt_free/W_pt_" + flavour[i] + flavour[j];
        mesPar_W.q1      = "W_Qpt_" + flavour[i];
        mesPar_W.q2      = "W_Qpt_" + flavour[j];
        mesPar_W.gammas  = "(Gamma5 Gamma5)";
        mesPar_W.sink    = "sink";
        application.createModule<MContraction::Meson>("W_meson_pt_"
                                                      + flavour[i] + flavour[j],
                                                      mesPar_W);

    }



    // execution
    application.saveParameterFile("free_prop.xml");
    application.run();
    
    // epilogue
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();
    
    return EXIT_SUCCESS;
}
