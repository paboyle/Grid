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

template <class Implementation, 
          template < typename, typename, typename > class Integrator, 
          class RepresentationsPolicy = NoHirep >
class BinaryHmcRunnerTemplate {
public:
 INHERIT_FIELD_TYPES(Implementation);
	    typedef Implementation ImplPolicy; // visible from outside
      template < typename S = NoSmearing<Implementation> > 
      using IntegratorType = Integrator<Implementation,S,RepresentationsPolicy>;

      enum StartType_t 
      { 
        ColdStart,
        HotStart,
        TepidStart,
        CheckpointStart,
        FilenameStart 
      };

      struct HMCPayload
      {
        StartType_t StartType;
        HMCparameters Parameters;

        HMCPayload() { StartType = HotStart; }
      };

     // These can be rationalised, some private
	  HMCPayload Payload; // Parameters 
	  HMCResourceManager Resources;
	  IntegratorParameters MDparameters;

	  ActionSet<Field, RepresentationsPolicy> TheAction;

    // A vector of HmcObservable that can be injected from outside
	  std::vector<HmcObservable<typename Implementation::Field> *> ObservablesList;

	  //GridCartesian *        UGrid;

    // These two are unnecessary, eliminate
	  // GridRedBlackCartesian *UrbGrid;
	  // GridCartesian *        FGrid;
	  // GridRedBlackCartesian *FrbGrid;

	  void ReadCommandLine(int argc, char ** argv) {
	  	std::string arg;

	  	if (GridCmdOptionExists(argv, argv + argc, "--StartType")) {
	  		arg = GridCmdOptionPayload(argv, argv + argc, "--StartType");
	  		if (arg == "HotStart") {
	  			Payload.StartType = HotStart;
	  		} else if (arg == "ColdStart") {
	  			Payload.StartType = ColdStart;
	  		} else if (arg == "TepidStart") {
	  			Payload.StartType = TepidStart;
	  		} else if (arg == "CheckpointStart") {
	  			Payload.StartType = CheckpointStart;
	  		} else {
	  			std::cout << GridLogError << "Unrecognized option in --StartType\n";
	  			std::cout
	  			<< GridLogError
	  			<< "Valid [HotStart, ColdStart, TepidStart, CheckpointStart]\n";
	  			assert(0);
	  		}
	  	}

	  	if (GridCmdOptionExists(argv, argv + argc, "--StartTrajectory")) {
	  		arg = GridCmdOptionPayload(argv, argv + argc, "--StartTrajectory");
	  		std::vector<int> ivec(0);
	  		GridCmdOptionIntVector(arg, ivec);
	  		Payload.Parameters.StartTrajectory = ivec[0];
	  	}

	  	if (GridCmdOptionExists(argv, argv + argc, "--Trajectories")) {
	  		arg = GridCmdOptionPayload(argv, argv + argc, "--Trajectories");
	  		std::vector<int> ivec(0);
	  		GridCmdOptionIntVector(arg, ivec);
	  		Payload.Parameters.Trajectories = ivec[0];
	  	}

	  	if (GridCmdOptionExists(argv, argv + argc, "--Thermalizations")) {
	  		arg = GridCmdOptionPayload(argv, argv + argc, "--Thermalizations");
	  		std::vector<int> ivec(0);
	  		GridCmdOptionIntVector(arg, ivec);
	  		Payload.Parameters.NoMetropolisUntil = ivec[0];
	  	}

	  }

  // A couple of wrapper functions
  template <class IOCheckpointer> void Run(IOCheckpointer &CP)  {
    NoSmearing<Implementation> S;
    Runner(CP, S);
  }

  template <class IOCheckpointer, class SmearingPolicy> void Run(IOCheckpointer &CP, SmearingPolicy &S) {
    Runner(CP, S);
  }

  //////////////////////////////////////////////////////////////////

private:
  template <class SmearingPolicy, class IOCheckpointer>
    void Runner(IOCheckpointer &Checkpoint,	SmearingPolicy &Smearing) {
  auto UGrid = Resources.GetCartesian();
  Resources.AddRNGs();
  Field U(UGrid);

  typedef IntegratorType<SmearingPolicy> TheIntegrator;
  TheIntegrator MDynamics(UGrid, MDparameters, TheAction, Smearing);

  if (Payload.StartType == HotStart) {
        // Hot start
   Payload.Parameters.MetropolisTest = true;
   Resources.SeedFixedIntegers();
   Implementation::HotConfiguration(Resources.GetParallelRNG(), U);
 } else if (Payload.StartType == ColdStart) {
        // Cold start
   Payload.Parameters.MetropolisTest = true;
   Resources.SeedFixedIntegers();
   Implementation::ColdConfiguration(Resources.GetParallelRNG(), U);
 } else if (Payload.StartType == TepidStart) {
        // Tepid start
   Payload.Parameters.MetropolisTest = true;
   Resources.SeedFixedIntegers();
   Implementation::TepidConfiguration(Resources.GetParallelRNG(), U);
 } else if (Payload.StartType == CheckpointStart) {
   Payload.Parameters.MetropolisTest = true;
        // CheckpointRestart
   Checkpoint.CheckpointRestore(Payload.Parameters.StartTrajectory, U, Resources.GetSerialRNG(), Resources.GetParallelRNG());
 }

 Smearing.set_Field(U);

 HybridMonteCarlo<TheIntegrator> HMC(Payload.Parameters, MDynamics, Resources.GetSerialRNG(), Resources.GetParallelRNG(), U);

 for (int obs = 0; obs < ObservablesList.size(); obs++)
   HMC.AddObservable(ObservablesList[obs]);

      // Run it
 HMC.evolve();
}
};

// These are for gauge fields, default integrator MinimumNorm2
template <template <typename, typename, typename> class Integrator >  using BinaryHmcRunner = BinaryHmcRunnerTemplate<PeriodicGimplR, Integrator > ;
template <template <typename, typename, typename> class Integrator >  using BinaryHmcRunnerF = BinaryHmcRunnerTemplate<PeriodicGimplF, Integrator > ;
template <template <typename, typename, typename> class Integrator >  using BinaryHmcRunnerD = BinaryHmcRunnerTemplate<PeriodicGimplD, Integrator > ;

template <class RepresentationsPolicy, template <typename, typename, typename> class Integrator >
using BinaryHmcRunnerTemplateHirep = BinaryHmcRunnerTemplate<PeriodicGimplR, Integrator, RepresentationsPolicy>;

typedef BinaryHmcRunnerTemplate<ScalarImplR, MinimumNorm2, ScalarFields> ScalarBinaryHmcRunner;

}  // namespace QCD
}  // namespace Grid

#endif  // GRID_GENERIC_HMC_RUNNER
