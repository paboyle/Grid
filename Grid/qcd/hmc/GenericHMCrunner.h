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

namespace Grid {
namespace QCD {


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
    }

    if (GridCmdOptionExists(argv, argv + argc, "--StartingTrajectory")) {
      arg = GridCmdOptionPayload(argv, argv + argc, "--StartingTrajectory");
      std::vector<int> ivec(0);
      GridCmdOptionIntVector(arg, ivec);
      Parameters.StartTrajectory = ivec[0];
    }

    if (GridCmdOptionExists(argv, argv + argc, "--Trajectories")) {
      arg = GridCmdOptionPayload(argv, argv + argc, "--Trajectories");
      std::vector<int> ivec(0);
      GridCmdOptionIntVector(arg, ivec);
      Parameters.Trajectories = ivec[0];
    }

    if (GridCmdOptionExists(argv, argv + argc, "--Thermalizations")) {
      arg = GridCmdOptionPayload(argv, argv + argc, "--Thermalizations");
      std::vector<int> ivec(0);
      GridCmdOptionIntVector(arg, ivec);
      Parameters.NoMetropolisUntil = ivec[0];
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

  //////////////////////////////////////////////////////////////////

 private:
  template <class SmearingPolicy>
  void Runner(SmearingPolicy &Smearing) {
    auto UGrid = Resources.GetCartesian();
    Resources.AddRNGs();
    Field U(UGrid);

    // Can move this outside?
    typedef IntegratorType<SmearingPolicy> TheIntegrator;
    TheIntegrator MDynamics(UGrid, Parameters.MD, TheAction, Smearing);

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
    }

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

template <class Implementation, class RepresentationsPolicy, 
          template <typename, typename, typename> class Integrator>
using GenericHMCRunnerTemplate = HMCWrapperTemplate<Implementation, Integrator, RepresentationsPolicy>;

typedef HMCWrapperTemplate<ScalarImplR, MinimumNorm2, ScalarFields>
    ScalarGenericHMCRunner;

typedef HMCWrapperTemplate<ScalarAdjImplR, MinimumNorm2, ScalarMatrixFields>
    ScalarAdjGenericHMCRunner;

template <int Colours> 
using ScalarNxNAdjGenericHMCRunner = HMCWrapperTemplate < ScalarNxNAdjImplR<Colours>, ForceGradient, ScalarNxNMatrixFields<Colours> >;

}  // namespace QCD
}  // namespace Grid

#endif  // GRID_GENERIC_HMC_RUNNER
