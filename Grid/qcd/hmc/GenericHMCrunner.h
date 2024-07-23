/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/hmc/GenericHmcRunner.h

Copyright (C) 2015

Author: paboyle <paboyle@ph.ed.ac.uk>
Author: Guido Cossu <guido.cossu@ed.ac.uk>

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
  directory
*************************************************************************************/
			   /*  END LEGAL */
#ifndef GRID_GENERIC_HMC_RUNNER
#define GRID_GENERIC_HMC_RUNNER

#include <unordered_map>

NAMESPACE_BEGIN(Grid);

// very ugly here but possibly resolved if we had a base Reader class
template < class ReaderClass >
class HMCRunnerBase {
public:
  virtual void Run() = 0;
  virtual void initialize(ReaderClass& ) = 0;
};

template <class Implementation,
          template <typename, typename, typename> class Integrator,
          class RepresentationsPolicy = NoHirep, class ReaderClass = XmlReader>
class HMCWrapperTemplate: public HMCRunnerBase<ReaderClass> {
public:
  INHERIT_FIELD_TYPES(Implementation);
  typedef Implementation ImplPolicy;  // visible from outside
  template <typename S = NoSmearing<Implementation> >
  using IntegratorType = Integrator<Implementation, S, RepresentationsPolicy>;

  HMCparameters Parameters;
  std::string ParameterFile;
  HMCResourceManager<Implementation> Resources;

  // The set of actions (keep here for lower level users, for now)
  ActionSet<Field, RepresentationsPolicy> TheAction;

  HMCWrapperTemplate() = default;

  HMCWrapperTemplate(HMCparameters Par){
    Parameters = Par;
  }

  void initialize(ReaderClass & TheReader){
    std::cout  << "Initialization of the HMC" << std::endl;
    Resources.initialize(TheReader);

    // eventually add smearing

    Resources.GetActionSet(TheAction);    
  }


  void ReadCommandLine(int argc, char **argv) {
    std::string arg;

    if (GridCmdOptionExists(argv, argv + argc, "--StartingType")) {
      arg = GridCmdOptionPayload(argv, argv + argc, "--StartingType");

      if (arg != "HotStart" && arg != "ColdStart" && arg != "TepidStart" &&
          arg != "CheckpointStart") {
        std::cout << GridLogError << "Unrecognized option in --StartingType\n";
        std::cout
	  << GridLogError
	  << "Valid [HotStart, ColdStart, TepidStart, CheckpointStart]\n";
        exit(1);
      }
      Parameters.StartingType = arg;
      std::cout <<GridLogMessage << " GenericHMCrunner --StartingType "<<arg<<std::endl;
    }

    if (GridCmdOptionExists(argv, argv + argc, "--StartingTrajectory")) {
      arg = GridCmdOptionPayload(argv, argv + argc, "--StartingTrajectory");
      std::vector<int> ivec(0);
      GridCmdOptionIntVector(arg, ivec);
      Parameters.StartTrajectory = ivec[0];
      std::cout <<GridLogMessage << " GenericHMCrunner --StartingTrajectory "<<ivec[0]<<std::endl;
    }

    if (GridCmdOptionExists(argv, argv + argc, "--Trajectories")) {
      arg = GridCmdOptionPayload(argv, argv + argc, "--Trajectories");
      std::vector<int> ivec(0);
      GridCmdOptionIntVector(arg, ivec);
      Parameters.Trajectories = ivec[0];
      std::cout << GridLogMessage<<" GenericHMCrunner Command Line --Trajectories "<<ivec[0]<<std::endl;
    }

    if (GridCmdOptionExists(argv, argv + argc, "--Thermalizations")) {
      arg = GridCmdOptionPayload(argv, argv + argc, "--Thermalizations");
      std::vector<int> ivec(0);
      GridCmdOptionIntVector(arg, ivec);
      Parameters.NoMetropolisUntil = ivec[0];
      std::cout << GridLogMessage<<" GenericHMCrunner --Thermalizations "<<ivec[0]<<std::endl;
    }
    if (GridCmdOptionExists(argv, argv + argc, "--ParameterFile")) {
      arg = GridCmdOptionPayload(argv, argv + argc, "--ParameterFile");
      ParameterFile = arg;
    }
  }


  template <class SmearingPolicy>
  void Run(SmearingPolicy &S) {
    Runner(S);
  }

  void Run(){
    NoSmearing<Implementation> S;
    Runner(S);
  }

  //Use the checkpointer to initialize the RNGs and the gauge field, writing the resulting gauge field into U.
  //This is called automatically by Run but may be useful elsewhere, e.g. for integrator tuning experiments
  void initializeGaugeFieldAndRNGs(Field &U){
    if(!Resources.haveRNGs()) Resources.AddRNGs();

    if (Parameters.StartingType == "HotStart") {
      // Hot start
      Resources.SeedFixedIntegers();
      Implementation::HotConfiguration(Resources.GetParallelRNG(), U);
    } else if (Parameters.StartingType == "ColdStart") {
      // Cold start
      Resources.SeedFixedIntegers();
      Implementation::ColdConfiguration(Resources.GetParallelRNG(), U);
    } else if (Parameters.StartingType == "TepidStart") {
      // Tepid start
      Resources.SeedFixedIntegers();
      Implementation::TepidConfiguration(Resources.GetParallelRNG(), U);
    } else if (Parameters.StartingType == "CheckpointStart") {
      // CheckpointRestart
      Resources.GetCheckPointer()->CheckpointRestore(Parameters.StartTrajectory, U,
						     Resources.GetSerialRNG(),
						     Resources.GetParallelRNG());
    } else if (Parameters.StartingType == "CheckpointStartReseed") {
      // Same as CheckpointRestart but reseed the RNGs using the fixed integer seeding used for ColdStart and HotStart
      // Useful for creating new evolution streams from an existing stream
      
      // WARNING: Unfortunately because the checkpointer doesn't presently allow us to separately restore the RNG and gauge fields we have to load
      // an existing RNG checkpoint first; make sure one is available and named correctly
      Resources.GetCheckPointer()->CheckpointRestore(Parameters.StartTrajectory, U,
						     Resources.GetSerialRNG(),
						     Resources.GetParallelRNG());
      Resources.SeedFixedIntegers();      
    } else {
      // others
      std::cout << GridLogError << "Unrecognized StartingType\n";
      std::cout
	<< GridLogError
	<< "Valid [HotStart, ColdStart, TepidStart, CheckpointStart, CheckpointStartReseed]\n";
      exit(1);
    }
  }



  //////////////////////////////////////////////////////////////////

private:
  template <class SmearingPolicy>
  void Runner(SmearingPolicy &Smearing) {
    auto UGrid = Resources.GetCartesian();
    Field U(UGrid);

    initializeGaugeFieldAndRNGs(U);

    typedef IntegratorType<SmearingPolicy> TheIntegrator;
    TheIntegrator MDynamics(UGrid, Parameters.MD, TheAction, Smearing);

    // Sets the momentum filter
    MDynamics.setMomentumFilter(*(Resources.GetMomentumFilter()));

    Smearing.set_Field(U);

    HybridMonteCarlo<TheIntegrator> HMC(Parameters, MDynamics,
                                        Resources.GetSerialRNG(),
                                        Resources.GetParallelRNG(), 
                                        Resources.GetObservables(), U);

    // Run it
    HMC.evolve();
  }
};

// These are for gauge fields, default integrator MinimumNorm2
template <template <typename, typename, typename> class Integrator>
using GenericHMCRunner = HMCWrapperTemplate<PeriodicGimplR, Integrator>;
template <template <typename, typename, typename> class Integrator>
using GenericHMCRunnerF = HMCWrapperTemplate<PeriodicGimplF, Integrator>;
template <template <typename, typename, typename> class Integrator>
using GenericHMCRunnerD = HMCWrapperTemplate<PeriodicGimplD, Integrator>;


// These are for gauge fields, default integrator MinimumNorm2
template <template <typename, typename, typename> class Integrator>
using ConjugateHMCRunner = HMCWrapperTemplate<ConjugateGimplR, Integrator>;
template <template <typename, typename, typename> class Integrator>
using ConjugateHMCRunnerF = HMCWrapperTemplate<ConjugateGimplF, Integrator>;
template <template <typename, typename, typename> class Integrator>
using ConjugateHMCRunnerD = HMCWrapperTemplate<ConjugateGimplD, Integrator>;



template <class RepresentationsPolicy,
          template <typename, typename, typename> class Integrator>
using GenericHMCRunnerHirep =
				     HMCWrapperTemplate<PeriodicGimplR, Integrator, RepresentationsPolicy>;

// sp2n

template <template <typename, typename, typename> class Integrator>
using GenericSpHMCRunner = HMCWrapperTemplate<SpPeriodicGimplR, Integrator>;

template <class RepresentationsPolicy,
          template <typename, typename, typename> class Integrator>
using GenericSpHMCRunnerHirep =
                     HMCWrapperTemplate<SpPeriodicGimplR, Integrator, RepresentationsPolicy>;



template <class Implementation, class RepresentationsPolicy, 
          template <typename, typename, typename> class Integrator>
using GenericHMCRunnerTemplate = HMCWrapperTemplate<Implementation, Integrator, RepresentationsPolicy>;

typedef HMCWrapperTemplate<ScalarImplR, MinimumNorm2, ScalarFields>
ScalarGenericHMCRunner;

typedef HMCWrapperTemplate<ScalarAdjImplR, MinimumNorm2, ScalarMatrixFields>
ScalarAdjGenericHMCRunner;

template <int Colours> 
using ScalarNxNAdjGenericHMCRunner = HMCWrapperTemplate < ScalarNxNAdjImplR<Colours>, ForceGradient, ScalarNxNMatrixFields<Colours> >;

NAMESPACE_END(Grid);

#endif  // GRID_GENERIC_HMC_RUNNER
