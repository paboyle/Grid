/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Tests/Hadrons/Test_hadrons_meson_3pt.cc

Copyright (C) 2015-2018

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
    std::vector<std::string> flavour = {"l", "s", "c1", "c2", "c3"};
    std::vector<double>      mass    = {.01, .04, .2  , .25 , .3  };
    unsigned int             nt      = GridDefaultLatt()[Tp];
    
    // global parameters
    Application::GlobalPar globalPar;
    globalPar.trajCounter.start    = 1500;
    globalPar.trajCounter.end      = 1520;
    globalPar.trajCounter.step     = 20;
    globalPar.runId                = "test";
    globalPar.genetic.maxGen       = 1000;
    globalPar.genetic.maxCstGen    = 200;
    globalPar.genetic.popSize      = 20;
    globalPar.genetic.mutationRate = .1;
    application.setPar(globalPar);
    
    // gauge field
    application.createModule<MGauge::Unit>("gauge");
    
    // set fermion boundary conditions to be periodic space, antiperiodic time.
    std::string boundary = "1 1 1 -1";
    std::string twist = "0. 0. 0. 0.";

    // sink
    MSink::Point::Par sinkPar;
    sinkPar.mom = "0 0 0";
    application.createModule<MSink::ScalarPoint>("sink", sinkPar);
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
    for (unsigned int t = 0; t < nt; t += 1)
    {
        std::string                           srcName;
        std::vector<std::string>              qName;
        std::vector<std::vector<std::string>> seqName;
        
        // Z2 source
        MSource::Z2::Par z2Par;
        z2Par.tA = t;
        z2Par.tB = t;
        srcName  = "z2_" + std::to_string(t);
        application.createModule<MSource::Z2>(srcName, z2Par);
        for (unsigned int i = 0; i < flavour.size(); ++i)
        {
            // sequential sources
            MSource::SeqGamma::Par seqPar;
            qName.push_back("QZ2_" + flavour[i] + "_" + std::to_string(t));
            seqPar.q   = qName[i];
            seqPar.tA  = (t + nt/4) % nt;
            seqPar.tB  = (t + nt/4) % nt;
            seqPar.mom = "1. 0. 0. 0.";
            seqName.push_back(std::vector<std::string>(Nd));
            for (unsigned int mu = 0; mu < Nd; ++mu)
            {
                seqPar.gamma   = 0x1 << mu;
                seqName[i][mu] = "G" + std::to_string(seqPar.gamma)
                                 + "_" + std::to_string(seqPar.tA) + "-"
                                 + qName[i];
                application.createModule<MSource::SeqGamma>(seqName[i][mu], seqPar);
            }
            
            // propagators
            MFermion::GaugeProp::Par quarkPar;
            quarkPar.solver = "CG_" + flavour[i];
            quarkPar.source = srcName;
            application.createModule<MFermion::GaugeProp>(qName[i], quarkPar);
            for (unsigned int mu = 0; mu < Nd; ++mu)
            {
                quarkPar.source = seqName[i][mu];
                seqName[i][mu]  = "Q_" + flavour[i] + "-" + seqName[i][mu];
                application.createModule<MFermion::GaugeProp>(seqName[i][mu], quarkPar);
            }
        }
        
        // contractions
        MContraction::Meson::Par mesPar;
        for (unsigned int i = 0; i < flavour.size(); ++i)
        for (unsigned int j = i; j < flavour.size(); ++j)
        {
            mesPar.output = "mesons/Z2_" + flavour[i] + flavour[j];
            mesPar.q1     = qName[i];
            mesPar.q2     = qName[j];
            mesPar.gammas = "all";
            mesPar.sink   = "sink";
            application.createModule<MContraction::Meson>("meson_Z2_"
                                                          + std::to_string(t)
                                                          + "_"
                                                          + flavour[i]
                                                          + flavour[j],
                                                          mesPar);
        }
        for (unsigned int i = 0; i < flavour.size(); ++i)
        for (unsigned int j = 0; j < flavour.size(); ++j)
        for (unsigned int mu = 0; mu < Nd; ++mu)
        {
            MContraction::Meson::Par mesPar;
            
            mesPar.output = "3pt/Z2_" + flavour[i] + flavour[j] + "_"
                            + std::to_string(mu);
            mesPar.q1     = qName[i];
            mesPar.q2     = seqName[j][mu];
            mesPar.gammas = "all";
            mesPar.sink   = "sink";
            application.createModule<MContraction::Meson>("3pt_Z2_"
                                                          + std::to_string(t)
                                                          + "_"
                                                          + flavour[i]
                                                          + flavour[j]
                                                          + "_"
                                                          + std::to_string(mu),
                                                          mesPar);
        }
    }
    
    // execution
    application.saveParameterFile("meson3pt.xml");
    application.run();
    
    // epilogue
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();
    
    return EXIT_SUCCESS;
}
