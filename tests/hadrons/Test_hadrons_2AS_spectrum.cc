/*******************************************************************************
 Grid physics library, www.github.com/paboyle/Grid
 
 Source file: tests/hadrons/Test_hadrons_spectrum.cc
 
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


    BEGIN_HADRONS_NAMESPACE
    BEGIN_MODULE_NAMESPACE(MFermion)
    MODULE_REGISTER_NS(GaugeProp2AS, TGaugeProp<WilsonTwoIndexAntiSymmetricImplR>, MFermion);
    END_MODULE_NAMESPACE
    BEGIN_MODULE_NAMESPACE(MSource)
    MODULE_REGISTER_NS(Point2AS, TPoint<WilsonTwoIndexAntiSymmetricImplR>, MSource);
    END_MODULE_NAMESPACE
    BEGIN_MODULE_NAMESPACE(MContraction)
    MODULE_REGISTER_NS(Meson2AS, ARG(TMeson<WilsonTwoIndexAntiSymmetricImplR, WilsonTwoIndexAntiSymmetricImplR>), MContraction);
//    MODULE_REGISTER_NS(BaryonMultirep, ARG(TBaryon<FIMPL, FIMPL, FIMPL>), MContraction);
    END_MODULE_NAMESPACE
    BEGIN_MODULE_NAMESPACE(MSink)
    MODULE_REGISTER_NS(ScalarPoint2AS, TPoint<WilsonTwoIndexAntiSymmetricImplR>, MSink);
    END_MODULE_NAMESPACE
    BEGIN_MODULE_NAMESPACE(MSolver)
    MODULE_REGISTER_NS(RBPrecCG2AS, TRBPrecCG<WilsonTwoIndexAntiSymmetricImplR>, MSolver);
    END_MODULE_NAMESPACE
    BEGIN_MODULE_NAMESPACE(MAction)
    MODULE_REGISTER_NS(WilsonClover2AS, TWilsonClover<WilsonTwoIndexAntiSymmetricImplR>, MAction);
    END_MODULE_NAMESPACE
    END_HADRONS_NAMESPACE


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
    std::vector<std::string> flavour = {"l", "s"};
    std::vector<double>      mass    = {-0.01, -0.04};
    double                   csw     = 1.0;
    // global parameters
    Application::GlobalPar globalPar;
    globalPar.trajCounter.start = 1500;
    globalPar.trajCounter.end   = 1520;
    globalPar.trajCounter.step  = 20;
    globalPar.seed              = "1 2 3 4";
    application.setPar(globalPar);
    // gauge field
    application.createModule<MGauge::Unit>("gauge");
    MSource::Point2AS::Par ptPar;
    ptPar.position = "0 0 0 0";
    application.createModule<MSource::Point2AS>("pt", ptPar);
    // sink
    MSink::ScalarPoint2AS::Par sinkPar;
    sinkPar.mom = "0 0 0";
    application.createModule<MSink::ScalarPoint2AS>("sink", sinkPar);
    
    // set fermion boundary conditions to be periodic space, antiperiodic time.
    std::string boundary = "1 1 1 -1";

    for (unsigned int i = 0; i < flavour.size(); ++i)
    {
        // actions
        MAction::WilsonClover2AS::Par actionPar;
        actionPar.gauge = "gauge";
        actionPar.mass  = mass[i];
        actionPar.csw_r = csw;
	    actionPar.csw_t = csw;
        actionPar.clover_anisotropy.isAnisotropic= false;
        actionPar.clover_anisotropy.t_direction  = Nd-1 ;
        actionPar.clover_anisotropy.xi_0         = 1.0  ;
        actionPar.clover_anisotropy.nu           = 1.0  ;
	    actionPar.boundary = boundary;
        application.createModule<MAction::WilsonClover2AS>("WilsonClover2AS_" + flavour[i], actionPar);
        
        // solvers
        MSolver::RBPrecCG2AS::Par solverPar;
        solverPar.action       = "WilsonClover2AS_" + flavour[i];
        solverPar.residual     = 1.0e-8;
        solverPar.maxIteration = 10000;
        application.createModule<MSolver::RBPrecCG2AS>("CG_" + flavour[i],
                                                    solverPar);
        
        // propagators
        MFermion::GaugeProp2AS::Par quarkPar;
        quarkPar.solver = "CG_" + flavour[i];
        quarkPar.source = "pt";
        application.createModule<MFermion::GaugeProp2AS>("Qpt_" + flavour[i], quarkPar);
        quarkPar.source = "z2";
        application.createModule<MFermion::GaugeProp2AS>("QZ2_" + flavour[i], quarkPar);
    }
    for (unsigned int i = 0; i < flavour.size(); ++i)
    for (unsigned int j = i; j < flavour.size(); ++j)
    {
        MContraction::Meson2AS::Par mesPar;
        
        mesPar.output  = "mesons2AS/pt_" + flavour[i] + flavour[j];
        mesPar.q1      = "Qpt_" + flavour[i];
        mesPar.q2      = "Qpt_" + flavour[j];
        mesPar.gammas  = "all";
        mesPar.sink    = "sink";
        application.createModule<MContraction::Meson2AS>("meson_pt_"
                                                      + flavour[i] + flavour[j],
                                                      mesPar);
     
     //   mesPar.output  = "mesons2AS/Z2_" + flavour[i] + flavour[j];
     //   mesPar.q1      = "QZ2_" + flavour[i];
     //   mesPar.q2      = "QZ2_" + flavour[j];
     //   mesPar.gammas  = "all";
     //   mesPar.sink    = "sink";
     //   application.createModule<MContraction::Meson2AS>("meson_Z2_"
     //                                                 + flavour[i] + flavour[j],
     //                                                 mesPar);
    }
    for (unsigned int i = 0; i < flavour.size(); ++i)
    for (unsigned int j = i; j < flavour.size(); ++j)
    for (unsigned int k = j; k < flavour.size(); ++k)
    {
        MContraction::Baryon::Par barPar;
        
        barPar.output = "baryons/pt_" + flavour[i] + flavour[j] + flavour[k];
        barPar.q1     = "Qpt_" + flavour[i];
        barPar.q2     = "Qpt_" + flavour[j];
        barPar.q3     = "Qpt_" + flavour[k];
        application.createModule<MContraction::Baryon>(
            "baryon_pt_" + flavour[i] + flavour[j] + flavour[k], barPar);
    }
    
    // execution
    application.saveParameterFile("spectrum.xml");
    application.run();
    
    // epilogue
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();
    
    return EXIT_SUCCESS;
}
