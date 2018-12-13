/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Tests/Hadrons/Test_hadrons.hpp

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

#include <Hadrons/Application.hpp>
#include <Hadrons/Modules.hpp>

using namespace Grid;
using namespace Hadrons;

/*******************************************************************************
 * Macros to reduce code duplication.
 ******************************************************************************/
// Common initialisation
#define HADRONS_DEFAULT_INIT \
    Grid_init(&argc, &argv); \
    HadronsLogError.Active(GridLogError.isActive()); \
    HadronsLogWarning.Active(GridLogWarning.isActive()); \
    HadronsLogMessage.Active(GridLogMessage.isActive()); \
    HadronsLogIterative.Active(GridLogIterative.isActive()); \
    HadronsLogDebug.Active(GridLogDebug.isActive()); \
    LOG(Message) << "Grid initialized" << std::endl;

#define HADRONS_DEFAULT_GLOBALS(application) \
{ \
    Application::GlobalPar globalPar;           \
    globalPar.trajCounter.start    = 1500;      \
    globalPar.trajCounter.end      = 1520;      \
    globalPar.trajCounter.step     = 20;        \
    globalPar.runId                = "test";    \
    globalPar.genetic.maxGen       = 1000;      \
    globalPar.genetic.maxCstGen    = 200;       \
    globalPar.genetic.popSize      = 20;        \
    globalPar.genetic.mutationRate = .1;        \
    application.setPar(globalPar);              \
}

// Useful definitions
#define ZERO_MOM "0. 0. 0. 0."
#define INIT_INDEX(s, n) (std::string(s) + "_" + std::to_string(n))
#define ADD_INDEX(s, n) (s + "_" + std::to_string(n))
#define LABEL_3PT(s, t1, t2) ADD_INDEX(INIT_INDEX(s, t1), t2)
#define LABEL_4PT(s, t1, t2, t3) ADD_INDEX(ADD_INDEX(INIT_INDEX(s, t1), t2), t3)
#define LABEL_4PT_NOISE(s, t1, t2, t3, nn) ADD_INDEX(ADD_INDEX(ADD_INDEX(INIT_INDEX(s, t1), t2), t3), nn)
#define LABEL_5D(s) s + "_5d";

// Wall source/sink macros
#define NAME_3MOM_WALL_SOURCE(t, mom) ("wall_" + std::to_string(t) + "_" + mom)
#define NAME_WALL_SOURCE(t) NAME_3MOM_WALL_SOURCE(t, ZERO_MOM)
#define NAME_POINT_SOURCE(pos) ("point_" + pos)

// Meson module "gammas" special values
#define ALL_GAMMAS "all"

#define MAKE_3MOM_WALL_PROP(tW, mom, propName, solver)\
{\
    std::string srcName = NAME_3MOM_WALL_SOURCE(tW, mom);\
    makeWallSource(application, srcName, tW, mom);\
    makePropagator(application, propName, srcName, solver);\
}

#define MAKE_WALL_PROP(tW, propName, solver)\
        MAKE_3MOM_WALL_PROP(tW, ZERO_MOM, propName, solver)

// Sequential source macros
#define MAKE_SEQUENTIAL_PROP(tS, qSrc, mom, seqPropName, solver, gamma)\
{\
    std::string srcName = seqPropName + "_src";\
    makeSequentialSource(application, srcName, qSrc, tS, gamma, mom);\
    makePropagator(application, seqPropName, srcName, solver);\
}

// Point source macros
#define MAKE_POINT_PROP(pos, propName, solver)\
{\
    std::string srcName = NAME_POINT_SOURCE(pos);\
    makePointSource(application, srcName, pos);\
    makePropagator(application, propName, srcName, solver);\
}

/*******************************************************************************
 * Action setups.
 ******************************************************************************/

/*******************************************************************************
 * Name: makeWilsonAction
 * Parameters: application - main application that stores modules.
 *             actionName  - name of action module to create.
 *             gaugeField  - gauge field module.     
 *             mass        - quark mass.
 *             boundary    - fermion boundary conditions (default to periodic
 *                           space, antiperiodic time).
 * Returns: None.
 ******************************************************************************/
inline void makeWilsonAction(Application &application, std::string actionName,
                             std::string &gaugeField, double mass,
                             std::string boundary = "1 1 1 -1")
{
    if (!(VirtualMachine::getInstance().hasModule(actionName)))
    {
        MAction::Wilson::Par actionPar;
        actionPar.gauge = gaugeField;
        actionPar.mass  = mass;
        actionPar.boundary = boundary;
        actionPar.twist = "0. 0. 0. 0.";
        application.createModule<MAction::Wilson>(actionName, actionPar);
    }
}

/*******************************************************************************
 * Name: makeDWFAction
 * Parameters: application - main application that stores modules.
 *             actionName  - name of action module to create.
 *             gaugeField  - gauge field module.     
 *             mass        - quark mass.
 *             M5          - domain wall height.
 *             Ls          - fifth dimension extent.
 *             boundary    - fermion boundary conditions (default to periodic
 *                           space, antiperiodic time).
 * Returns: None.
 ******************************************************************************/
inline void makeDWFAction(Application &application, std::string actionName,
                          std::string &gaugeField, double mass, double M5,
                          unsigned int Ls, std::string boundary = "1 1 1 -1")
{
    if (!(VirtualMachine::getInstance().hasModule(actionName)))
    {
        MAction::DWF::Par actionPar;
        actionPar.gauge = gaugeField;
        actionPar.Ls    = Ls;
        actionPar.M5    = M5;
        actionPar.mass  = mass;
        actionPar.boundary = boundary;
        actionPar.twist = "0. 0. 0. 0.";
        application.createModule<MAction::DWF>(actionName, actionPar);
    }
}

/*******************************************************************************
 * Functions for propagator construction.
 ******************************************************************************/
 
/*******************************************************************************
 * Name: makeRBPrecCGSolver
 * Purpose: Make RBPrecCG solver module for specified action.
 * Parameters: application - main application that stores modules.
 *             solverName  - name of solver module to create.
 *             actionName  - action module corresponding to propagators to be
 *                           computed.
 *             residual    - CG target residual.
 * Returns: None.
 ******************************************************************************/
inline void makeRBPrecCGSolver(Application &application, std::string &solverName,
                               std::string &actionName, double residual = 1e-8)
{
    if (!(VirtualMachine::getInstance().hasModule(solverName)))
    {
        MSolver::RBPrecCG::Par solverPar;
        solverPar.action       = actionName;
        solverPar.residual     = residual;
        solverPar.maxIteration = 10000;
        application.createModule<MSolver::RBPrecCG>(solverName,
                                                    solverPar);
    }
}

/*******************************************************************************
 * Name: makePointSource
 * Purpose: Construct point source and add to application module.
 * Parameters: application - main application that stores modules.
 *             srcName     - name of source module to create.
 *             pos         - Position of point source.
 * Returns: None.
 ******************************************************************************/
inline void makePointSource(Application &application, std::string srcName,
                            std::string pos)
{
    // If the source already exists, don't make the module again.
    if (!(VirtualMachine::getInstance().hasModule(srcName)))
    {
        MSource::Point::Par pointPar;
        pointPar.position = pos;
        application.createModule<MSource::Point>(srcName, pointPar);
    }
}

/*******************************************************************************
 * Name: makeSequentialSource
 * Purpose: Construct sequential source and add to application module.
 * Parameters: application - main application that stores modules.
 *             srcName     - name of source module to create.
 *             qSrc        - Input quark for sequential inversion.
 *             tS          - sequential source timeslice.
 *             mom         - momentum insertion (default is zero).
 * Returns: None.
 ******************************************************************************/
inline void makeSequentialSource(Application &application, std::string srcName,
                                 std::string qSrc, unsigned int tS,
                                 Gamma::Algebra gamma = Gamma::Algebra::GammaT,
                                 std::string mom = ZERO_MOM)
{
    // If the source already exists, don't make the module again.
    if (!(VirtualMachine::getInstance().hasModule(srcName)))
    {
        MSource::SeqGamma::Par seqPar;
        seqPar.q   = qSrc;
        seqPar.tA  = tS;
        seqPar.tB  = tS;
        seqPar.mom = mom;
        seqPar.gamma = gamma;
        application.createModule<MSource::SeqGamma>(srcName, seqPar);
    }
}

/*******************************************************************************
 * Name: makeConservedSequentialSource
 * Purpose: Construct sequential source with conserved current insertion and 
 *          add to application module.
 * Parameters: application - main application that stores modules.
 *             srcName     - name of source module to create.
 *             qSrc        - Input quark for sequential inversion.
 *             actionName  - action corresponding to quark.
 *             tS          - sequential source timeslice.
 *             curr        - conserved current type to insert.
 *             mu          - Lorentz index of current to insert.
 *             mom         - momentum insertion (default is zero).
 * Returns: None.
 ******************************************************************************/
inline void makeConservedSequentialSource(Application &application,
                                          std::string &srcName,
                                          std::string &qSrc, 
                                          std::string &actionName,
                                          unsigned int tS,
                                          Current curr,
                                          unsigned int mu,
                                          std::string mom = ZERO_MOM)
{
    // If the source already exists, don't make the module again.
    if (!(VirtualMachine::getInstance().hasModule(srcName)))
    {
        MSource::SeqConserved::Par seqPar;
        seqPar.q         = qSrc;
        seqPar.action    = actionName;
        seqPar.tA        = tS;
        seqPar.tB        = tS;
        seqPar.curr_type = curr;
        seqPar.mu_min    = mu;
        seqPar.mu_min    = mu;
        seqPar.mom       = mom;
        application.createModule<MSource::SeqConserved>(srcName, seqPar);
    }
}

/*******************************************************************************
 * Name: makeNoiseSource
 * Parameters: application - main application that stores modules.
 *             srcName     - name of source module to create.
 *             tA          - lower source timeslice limit.
 *             tB          - upper source timeslice limit.
 * Returns: None.
 ******************************************************************************/
inline void makeNoiseSource(Application &application, std::string &srcName,
                            unsigned int tA, unsigned int tB)
{
    if (!(VirtualMachine::getInstance().hasModule(srcName)))
    {
        MSource::Z2::Par noisePar;
        noisePar.tA = tA;
        noisePar.tB = tB;
        application.createModule<MSource::Z2>(srcName, noisePar);
    }
 }

/*******************************************************************************
 * Name: makeWallSource
 * Purpose: Construct wall source and add to application module.
 * Parameters: application - main application that stores modules.
 *             srcName     - name of source module to create.
 *             tW          - wall source timeslice.
 *             mom         - momentum insertion (default is zero).
 * Returns: None.
 ******************************************************************************/
inline void makeWallSource(Application &application, std::string &srcName,
                           unsigned int tW, std::string mom = ZERO_MOM)
{
    // If the source already exists, don't make the module again.
    if (!(VirtualMachine::getInstance().hasModule(srcName)))
    {
        MSource::Wall::Par wallPar;
        wallPar.tW  = tW;
        wallPar.mom = mom;
        application.createModule<MSource::Wall>(srcName, wallPar);
    }
}

/*******************************************************************************
 * Name: makePointSink
 * Purpose: Create function for point sink smearing of a propagator.
 * Parameters: application - main application that stores modules.
 *             propName    - name of input propagator.
 *             sinkFnct    - name of output sink smearing module.
 *             mom         - momentum insertion (default is zero).
 * Returns: None.
 ******************************************************************************/
inline void makePointSink(Application &application, std::string &sinkFnct,
                          std::string mom = ZERO_MOM)
{
    // If the sink function already exists, don't make it again.
    if (!(VirtualMachine::getInstance().hasModule(sinkFnct)))
    {
        MSink::Point::Par pointPar;
        pointPar.mom = mom;
        application.createModule<MSink::Point>(sinkFnct, pointPar);
    }
}

/*******************************************************************************
 * Name: sinkSmear
 * Purpose: Perform sink smearing of a propagator.
 * Parameters: application - main application that stores modules.
 *             sinkFnct    - sink smearing module.
 *             propName    - propagator to smear.
 *             smearedProp - name of output smeared propagator.
 * Returns: None.
 ******************************************************************************/
inline void sinkSmear(Application &application, std::string &sinkFnct,
                      std::string &propName, std::string &smearedProp)
{
    // If the propagator has already been smeared, don't smear it again.
    if (!(VirtualMachine::getInstance().hasModule(smearedProp)))
    {
        MSink::Smear::Par smearPar;
        smearPar.q    = propName;
        smearPar.sink = sinkFnct;
        application.createModule<MSink::Smear>(smearedProp, smearPar);
    }
}

/*******************************************************************************
 * Name: makePropagator
 * Purpose: Construct source and propagator then add to application module.
 * Parameters: application - main application that stores modules.
 *             propName    - name of propagator module to create.
 *             srcName     - name of source module to use.
 *             solver      - solver to use (default is CG).
 * Returns: None.
 ******************************************************************************/
inline void makePropagator(Application &application, std::string &propName,
                           std::string &srcName, std::string &solver)
{
    // If the propagator already exists, don't make the module again.
    if (!(VirtualMachine::getInstance().hasModule(propName)))
    {
        MFermion::GaugeProp::Par quarkPar;
        quarkPar.source = srcName;
        quarkPar.solver = solver;
        application.createModule<MFermion::GaugeProp>(propName, quarkPar);
    }
}

/*******************************************************************************
 * Name: makeLoop
 * Purpose: Use noise source and inversion result to make loop propagator, then 
 *          add to application module.
 * Parameters: application - main application that stores modules.
 *             propName    - name of propagator module to create.
 *             srcName     - name of noise source module to use.
 *             resName     - name of inversion result on given noise source.
 * Returns: None.
 ******************************************************************************/
inline void makeLoop(Application &application, std::string &propName,
                     std::string &srcName, std::string &resName)
{
    // If the loop propagator already exists, don't make the module again.
    if (!(VirtualMachine::getInstance().hasModule(propName)))
    {
        MLoop::NoiseLoop::Par loopPar;
        loopPar.q   = resName;
        loopPar.eta = srcName;
        application.createModule<MLoop::NoiseLoop>(propName, loopPar);
    }
}

/*******************************************************************************
 * Contraction module creation.
 ******************************************************************************/

/*******************************************************************************
 * Name: mesonContraction
 * Purpose: Create meson contraction module and add to application module.
 * Parameters: application - main application that stores modules.
 *             modName     - unique module name.
 *             output      - name of output files.
 *             q1          - quark propagator 1.
 *             q2          - quark propagator 2.
 *             sink        - sink smearing module.
 *             gammas      - gamma insertions at source and sink.
 * Returns: None.
 ******************************************************************************/
inline void mesonContraction(Application &application, 
                             std::string &modName, std::string &output, 
                             std::string &q1, std::string &q2,
                             std::string &sink,
                             std::string gammas = "<Gamma5 Gamma5>")
{
    if (!(VirtualMachine::getInstance().hasModule(modName)))
    {
        MContraction::Meson::Par mesPar;
        mesPar.output = output;
        mesPar.q1 = q1;
        mesPar.q2 = q2;
        mesPar.sink = sink;
        mesPar.gammas = gammas;
        application.createModule<MContraction::Meson>(modName, mesPar);
    }
 }

/*******************************************************************************
 * Name: gamma3ptContraction
 * Purpose: Create gamma3pt contraction module and add to application module.
 * Parameters: application - main application that stores modules.
 *             npt         - specify n-point correlator (for labelling).
 *             q1          - quark propagator 1, sink smeared.
 *             q2          - quark propagator 2.
 *             q3          - quark propagator 3.
 *             label       - unique label to construct module name.
 *             tSnk        - sink position of sink for q1.
 *             gamma       - gamma insertions between q2 and q3.
 * Returns: None.
 ******************************************************************************/
inline void gamma3ptContraction(Application &application, unsigned int npt, 
                                std::string &q1, std::string &q2,
                                std::string &q3, std::string &label,
                                unsigned int tSnk = 0,
                                Gamma::Algebra gamma = Gamma::Algebra::Identity)
{
    std::string modName = std::to_string(npt) + "pt_" + label;
    if (!(VirtualMachine::getInstance().hasModule(modName)))
    {
        MContraction::Gamma3pt::Par gamma3ptPar;
        gamma3ptPar.output = std::to_string(npt) + "pt/" + label;
        gamma3ptPar.q1 = q1;
        gamma3ptPar.q2 = q2;
        gamma3ptPar.q3 = q3;
        gamma3ptPar.tSnk = tSnk;
        gamma3ptPar.gamma = gamma;
        application.createModule<MContraction::Gamma3pt>(modName, gamma3ptPar);
    }
 }

/*******************************************************************************
 * Name: weakContraction[Eye,NonEye]
 * Purpose: Create Weak Hamiltonian contraction module for Eye/NonEye topology
 *          and add to application module.
 * Parameters: application - main application that stores modules.
 *             npt         - specify n-point correlator (for labelling).
 *             q1          - quark propagator 1.
 *             q2          - quark propagator 2.
 *             q3          - quark propagator 3.
 *             q4          - quark propagator 4.
 *             label       - unique label to construct module name.
 *             tSnk        - time position of sink (for sink smearing).
 * Returns: None.
 ******************************************************************************/
#define HW_CONTRACTION(top) \
inline void weakContraction##top(Application &application, unsigned int npt,\
                                 std::string &q1, std::string &q2, \
                                 std::string &q3, std::string &q4, \
                                 std::string &label, unsigned int tSnk = 0)\
{\
    std::string modName = std::to_string(npt) + "pt_" + label;\
    if (!(VirtualMachine::getInstance().hasModule(modName)))\
    {\
        MContraction::WeakHamiltonian##top::Par weakPar;\
        weakPar.output = std::to_string(npt) + "pt/" + label;\
        weakPar.q1 = q1;\
        weakPar.q2 = q2;\
        weakPar.q3 = q3;\
        weakPar.q4 = q4;\
        weakPar.tSnk = tSnk;\
        application.createModule<MContraction::WeakHamiltonian##top>(modName, weakPar);\
    }\
}
HW_CONTRACTION(Eye)    // weakContractionEye
HW_CONTRACTION(NonEye) // weakContractionNonEye

/*******************************************************************************
 * Name: disc0Contraction
 * Purpose: Create contraction module for 4pt Weak Hamiltonian + current
 *          disconnected topology for neutral mesons and add to application 
 *          module.
 * Parameters: application - main application that stores modules.
 *             q1          - quark propagator 1.
 *             q2          - quark propagator 2.
 *             q3          - quark propagator 3.
 *             q4          - quark propagator 4.
 *             label       - unique label to construct module name.
 * Returns: None.
 ******************************************************************************/
inline void disc0Contraction(Application &application, 
                             std::string &q1, std::string &q2,
                             std::string &q3, std::string &q4,
                             std::string &label)
{
    std::string modName = "4pt_" + label;
    if (!(VirtualMachine::getInstance().hasModule(modName)))
    {
        MContraction::WeakNeutral4ptDisc::Par disc0Par;
        disc0Par.output = "4pt/" + label;
        disc0Par.q1 = q1;
        disc0Par.q2 = q2;
        disc0Par.q3 = q3;
        disc0Par.q4 = q4;
        application.createModule<MContraction::WeakNeutral4ptDisc>(modName, disc0Par);
    }
 }

/*******************************************************************************
 * Name: discLoopContraction
 * Purpose: Create contraction module for disconnected loop and add to
 *          application module.
 * Parameters: application - main application that stores modules.
 *             q_loop      - loop quark propagator.
 *             modName     - unique module name.
 *             gamma       - gamma matrix to use in contraction.
 * Returns: None.
 ******************************************************************************/
inline void discLoopContraction(Application &application,
                                std::string &q_loop, std::string &modName,
                                Gamma::Algebra gamma = Gamma::Algebra::Identity)
{
    if (!(VirtualMachine::getInstance().hasModule(modName)))
    {
        MContraction::DiscLoop::Par discPar;
        discPar.output = "disc/" + modName;
        discPar.q_loop = q_loop;
        discPar.gamma  = gamma;
        application.createModule<MContraction::DiscLoop>(modName, discPar);
    }
}

/*******************************************************************************
 * Name: makeWITest
 * Purpose: Create module to test Ward Identities for conserved current
 *          contractions and add to application module.
 * Parameters: application - main application that stores modules.
 *             modName     - name of module to create.
 *             propName    - 4D quark propagator.
 *             actionName  - action used to compute quark propagator.
 *             mass        - mass of quark.
 *             Ls          - length of 5th dimension (default = 1).
 *             test_axial  - whether or not to check PCAC relation.
 * Returns: None.
 ******************************************************************************/
inline void makeWITest(Application &application, std::string &modName,
                       std::string &propName, std::string &actionName, 
                       double mass, unsigned int Ls = 1, bool test_axial = false)
{
    if (!(VirtualMachine::getInstance().hasModule(modName)))
    {
        MContraction::WardIdentity::Par wiPar;
        if (Ls > 1)
        {
            wiPar.q = LABEL_5D(propName);
        }
        else
        {
            wiPar.q = propName;
        }
        wiPar.action     = actionName;
        wiPar.mass       = mass;
        wiPar.test_axial = test_axial;
        application.createModule<MContraction::WardIdentity>(modName, wiPar);
    }
}

/*******************************************************************************
 * Name: makeSeqCurrComparison
 * Purpose: Create module to compare sequential insertion of conserved current
 *          against sink contraction and add to application module.
 * Parameters: application - main application that stores modules.
 *             modName     - name of module to create.
 *             propName    - quark propagator (point source), 5D if available.
 *             seqName     - 4D quark propagator with sequential insertion of
 *                           conserved current.
 *             actionName  - action used to compute quark propagators.
 *             origin      - origin of point source propagator.
 *             t_J         - time at which sequential current is inserted.
 *             mu          - Lorentz index of sequential current.
 *             curr        - type of conserved current inserted.
 * Returns: None.
 ******************************************************************************/
inline void makeSeqCurrComparison(Application &application, std::string &modName,
                                 std::string &propName, std::string &seqName,
                                 std::string &actionName, std::string &origin,
                                 unsigned int t_J, unsigned int mu, Current curr)
{
    if (!(VirtualMachine::getInstance().hasModule(modName)))
    {
        MUtilities::TestSeqConserved::Par seqPar;
        seqPar.q      = propName;
        seqPar.qSeq   = seqName;
        seqPar.action = actionName;
        seqPar.origin = origin;
        seqPar.t_J    = t_J;
        seqPar.mu     = mu;
        seqPar.curr   = curr;
        application.createModule<MUtilities::TestSeqConserved>(modName, seqPar);
    }
}

/*******************************************************************************
 * Name: makeSeqGamComparison
 * Purpose: Create module to compare sequential insertion of gamma matrix
 *          against sink contraction and add to application module.
 * Parameters: application - main application that stores modules.
 *             modName     - name of module to create.
 *             propName    - 4D quark propagator.
 *             seqProp     - 4D quark propagator with sequential insertion of
 *                           gamma matrix.
 *             gamma       - Inserted gamma matrix.
 *             t_g         - time at which gamma matrix is inserted 
 *                           sequentially.
 * Returns: None.
 ******************************************************************************/
inline void makeSeqGamComparison(Application &application, std::string &modName,
                                 std::string &propName, std::string &seqProp,
                                 std::string &origin, Gamma::Algebra gamma, 
                                 unsigned int t_g)
{
    if (!(VirtualMachine::getInstance().hasModule(modName)))
    {
        MUtilities::TestSeqGamma::Par seqPar;
        seqPar.q      = propName;
        seqPar.qSeq   = seqProp;
        seqPar.origin = origin;
        seqPar.t_g    = t_g;
        seqPar.gamma  = gamma;
        application.createModule<MUtilities::TestSeqGamma>(modName, seqPar);
    }
}
