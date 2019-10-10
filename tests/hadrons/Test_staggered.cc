/*******************************************************************************

This program computes staggered two-, three-, and four-point functions using
sequential propagators. It is supposed to be a check for correctness of other
methods, e.g., A2A methods using meson fields. What we compute concretely are 
the following correlation functions:

Two-point correlator <P(t)P(=0)>

    t=0                     t={0,1,2,3}
    |P----------->-----------P|
    |------------<------------|

Three-point correlator <P(t)P(t=1)P(t=0)>

    t=0     t=1             t={0,1,2,3}
    |P--->---P--->-----------P|
    |------------<------------|

Four-point correlator <P(t)P(t=2)P(t=1)P(t=0)>

    t=0     t=1     t=2     t={0,1,2,3}
    |P--->---P--->---P--->---P|
    |------------<------------|

In these expressions P is a pseudoscalar current.

The calculation requires three propagators:
 * quark_0t: "normal" propagator from the origin to timeslice t
 * quark_01t: 0 to 1 to t via a sequential propagator
 * quark_012t: 0 to 1 to 2 to t via a sequential propagator

*******************************************************************************/

// We need to preempt the defaul fermion stuff,
// which gets loaded from Global.hpp by Application.hpp
// #ifndef FIMPLBASE
// #define FIMPLBASE WilsonImpl
// #endif

// #ifndef FIMPLBASE
// #define FIMPLBASE StaggeredImpl
// #endif

// #ifndef ZFIMPLBASE
// #define ZFIMPLBASE ZWilsonImpl
// #endif

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
    unsigned int             nt      = GridDefaultLatt()[Tp];
    
    // global parameters
    Application::GlobalPar globalPar;
    globalPar.trajCounter.start    = 20;
    globalPar.trajCounter.end      = 40;
    globalPar.trajCounter.step     = 20;
    globalPar.runId                = "test";
    globalPar.genetic.maxGen       = 1000;
    globalPar.genetic.maxCstGen    = 200;
    globalPar.genetic.popSize      = 20;
    globalPar.genetic.mutationRate = 0.1;
    application.setPar(globalPar);
    
    // gauge field
    // application.createModule<MGauge::Unit>("gauge");
    // Load a stored gauge configuration -- NERSC style
    MIO::LoadNerscPar ldNerscPar;
    ldNerscPar.file = "config/ckpoint_lat";
    application.createModule<MIO::LoadNersc>("gauge", ldNerscPar);

    // set fermion boundary conditions to be periodic space, antiperiodic time.
    std::string boundary = "1 1 1 -1";
    std::string twist = "0. 0. 0. 0.";

    // actions
    MAction::ImprovedStaggered::Par actionPar;
    actionPar.gauge = "gauge";
    actionPar.mass = 0.1;
    actionPar.c1 = 1.0;
    actionPar.c2 = 0.0;
    actionPar.tad = 1.0;
    actionPar.boundary = boundary;
    actionPar.twist = twist;
    std::string actionName = "stag";
    application.createModule<MAction::ImprovedStaggered>(actionName, actionPar);

    // solvers
    MSolver::RBPrecCG::Par solverPar;
    solverPar.action       = "stag";
    solverPar.residual     = 1.0e-8;
    solverPar.maxIteration = 10000;
    std::string solverName = "cg_stag"; 
    application.createModule<MSolver::StagRBPrecCG>(solverName, solverPar);    

    std::string srcName;
    std::string quarkName;
    MFermion::GaugeProp::Par quarkPar;

    // Point source at the origin
    MSource::Point::Par pointPar;
    pointPar.position = "0 0 0 0";
    srcName = "point_origin";
    application.createModule<MSource::StagPoint>(srcName, pointPar);

    // scalar sink, no momentum injected
    MSink::Point::Par sinkPar;
    sinkPar.mom = "0 0 0";
    std::string sinkName = "sink";
    application.createModule<MSink::ScalarPoint>(sinkName, sinkPar);

    // Propagator from origin to (x,y,z,t) via point source at origin
    quarkName       = "quark_0t";
    quarkPar.solver = "cg_stag";
    quarkPar.source = "point_origin";
    application.createModule<MFermion::StagGaugeProp>(quarkName, quarkPar);

    // First sequential propagator (origin) --> (t=1) --> (x,y,z,t)
    // Sequential source
    std::string seqSrcName;
    MSource::SeqGamma::Par seqPar;
    seqPar.q     = "quark_0t";
    seqPar.tA    = 1; // source has support on timeslice t=1 only, so tA=tB=1
    seqPar.tB    = 1; // source has support on timeslice t=1 only, so tA=tB=1
    seqPar.mom   = "0. 0. 0. 0.";
    seqPar.gamma = 1; // i.e., Gamma5, see Gamma.h
    seqSrcName   = "seq_1";
    application.createModule<MSource::StagSeqGamma>(seqSrcName, seqPar);
    // Sequential propagator
    quarkName       = "quark_01t";
    quarkPar.solver = "cg_stag";
    quarkPar.source = "seq_1"; 
    application.createModule<MFermion::StagGaugeProp>(quarkName, quarkPar);

    // Second sequential propagator (origin) --> (t=1) --> (t=2) --> (x,y,z,t)
    seqPar.q     = "quark_01t";
    seqPar.tA    = 2; // source has support on timeslice t=2 only, so tA=tB=1
    seqPar.tB    = 2; // source has support on timeslice t=1 only, so tA=tB=1
    seqPar.mom   = "0. 0. 0. 0.";
    seqPar.gamma = 1; // i.e., Gamma5, see Gamma.h
    seqSrcName   = "seq_2";
    application.createModule<MSource::StagSeqGamma>(seqSrcName, seqPar);
    // Sequential propagator
    quarkName       = "quark_012t";
    quarkPar.solver = "cg_stag";
    quarkPar.source = "seq_2"; 
    application.createModule<MFermion::StagGaugeProp>(quarkName, quarkPar);

    // Two-point function
    std::string mesName;
    MContraction::Meson::Par mesPar;
    mesPar.output = "2pt/point_0t";
    mesPar.q1     = "quark_0t";
    mesPar.q2     = "quark_0t";
    mesPar.gammas = "Gamma5"; // diagonal
    mesPar.sink   = "sink";
    mesName       = "2pt_point";
    application.createModule<MContraction::StagMeson>(mesName, mesPar);

    // Three-point function
    mesPar.output = "3pt/point_01t";
    mesPar.q1     = "quark_0t";
    mesPar.q2     = "quark_01t";
    mesPar.gammas = "Gamma5"; // diagonal
    mesPar.sink   = "sink";
    mesName       = "3pt_point";
    application.createModule<MContraction::StagMeson>(mesName, mesPar);

    // Four-point function
    mesPar.output = "4pt/point_012t";
    mesPar.q1     = "quark_0t";
    mesPar.q2     = "quark_012t";
    mesPar.gammas = "Gamma5"; // diagonal
    mesPar.sink   = "sink";
    mesName       = "4pt_point";
    application.createModule<MContraction::StagMeson>(mesName, mesPar);

    // execution
    application.saveParameterFile("meson4ptV2.xml");
    application.run();
    
    // epilogue
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();
    
    return EXIT_SUCCESS;
}