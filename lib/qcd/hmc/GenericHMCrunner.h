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
#ifndef GENERIC_HMC_RUNNER
#define GENERIC_HMC_RUNNER

namespace Grid {
namespace QCD {

  // Virtual Class for HMC specific for gauge theories
  // implement a specific theory by defining the BuildTheAction
  template <class Implementation, class RepresentationsPolicy = NoHirep>
  class BinaryHmcRunnerTemplate {
  public:
    INHERIT_FIELD_TYPES(Implementation);
    typedef Implementation ImplPolicy;

    enum StartType_t { ColdStart,
                       HotStart,
                       TepidStart,
                       CheckpointStart };

    ActionSet<Field, RepresentationsPolicy> TheAction;

    // A vector of HmcObservable
    // that can be injected from outside
    std::vector<HmcObservable<typename Implementation::Field> *>
        ObservablesList;

    IntegratorParameters MDparameters;

    GridCartesian *        UGrid;
    GridCartesian *        FGrid;
    GridRedBlackCartesian *UrbGrid;
    GridRedBlackCartesian *FrbGrid;

    std::vector<int> SerialSeed;
    std::vector<int> ParallelSeed;

    void RNGSeeds(std::vector<int> S, std::vector<int> P) {
      SerialSeed   = S;
      ParallelSeed = P;
    }

    virtual void BuildTheAction(int argc, char **argv) = 0; // necessary?

    // A couple of wrapper classes
    template <class IOCheckpointer>
    void Run(int argc, char **argv, IOCheckpointer &Checkpoint) {
      NoSmearing<Implementation> S;
      Runner(argc, argv, Checkpoint, S);
    }

    template <class IOCheckpointer, class SmearingPolicy>
    void Run(int argc, char **argv, IOCheckpointer &CP, SmearingPolicy &S) {
      Runner(argc, argv, CP, S);
    }
    //////////////////////////////

    template <class SmearingPolicy, class IOCheckpointer>
    void Runner(int             argc,
                char **         argv,
                IOCheckpointer &Checkpoint,
                SmearingPolicy &Smearing) {
      StartType_t StartType = HotStart;

      std::string arg;

      if (GridCmdOptionExists(argv, argv + argc, "--StartType")) {
        arg = GridCmdOptionPayload(argv, argv + argc, "--StartType");
        if (arg == "HotStart") {
          StartType = HotStart;
        } else if (arg == "ColdStart") {
          StartType = ColdStart;
        } else if (arg == "TepidStart") {
          StartType = TepidStart;
        } else if (arg == "CheckpointStart") {
          StartType = CheckpointStart;
        } else {
          std::cout << GridLogError << "Unrecognized option in --StartType\n";
          std::cout
              << GridLogError
              << "Valid [HotStart, ColdStart, TepidStart, CheckpointStart]\n";
          assert(0);
        }
      }

      int StartTraj = 0;
      if (GridCmdOptionExists(argv, argv + argc, "--StartTrajectory")) {
        arg = GridCmdOptionPayload(argv, argv + argc, "--StartTrajectory");
        std::vector<int> ivec(0);
        GridCmdOptionIntVector(arg, ivec);
        StartTraj = ivec[0];
      }

      int NumTraj = 1;
      if (GridCmdOptionExists(argv, argv + argc, "--Trajectories")) {
        arg = GridCmdOptionPayload(argv, argv + argc, "--Trajectories");
        std::vector<int> ivec(0);
        GridCmdOptionIntVector(arg, ivec);
        NumTraj = ivec[0];
      }

      int NumThermalizations = 10;
      if (GridCmdOptionExists(argv, argv + argc, "--Thermalizations")) {
        arg = GridCmdOptionPayload(argv, argv + argc, "--Thermalizations");
        std::vector<int> ivec(0);
        GridCmdOptionIntVector(arg, ivec);
        NumThermalizations = ivec[0];
      }

      GridSerialRNG   sRNG;
      GridParallelRNG pRNG(UGrid);
      Field           U(UGrid);

      typedef MinimumNorm2<Implementation,
                           SmearingPolicy,
                           RepresentationsPolicy>
          IntegratorType; // change here to change the algorithm
      IntegratorType MDynamics(UGrid, MDparameters, TheAction, Smearing);

      HMCparameters HMCpar;
      HMCpar.StartTrajectory   = StartTraj;
      HMCpar.Trajectories      = NumTraj;
      HMCpar.NoMetropolisUntil = NumThermalizations;

      if (StartType == HotStart) {
        // Hot start
        HMCpar.MetropolisTest = true;
        sRNG.SeedFixedIntegers(SerialSeed);
        pRNG.SeedFixedIntegers(ParallelSeed);
        Implementation::HotConfiguration(pRNG, U);
      } else if (StartType == ColdStart) {
        // Cold start
        HMCpar.MetropolisTest = true;
        sRNG.SeedFixedIntegers(SerialSeed);
        pRNG.SeedFixedIntegers(ParallelSeed);
        Implementation::ColdConfiguration(pRNG, U);
      } else if (StartType == TepidStart) {
        // Tepid start
        HMCpar.MetropolisTest = true;
        sRNG.SeedFixedIntegers(SerialSeed);
        pRNG.SeedFixedIntegers(ParallelSeed);
        Implementation::TepidConfiguration(pRNG, U);
      } else if (StartType == CheckpointStart) {
        HMCpar.MetropolisTest = true;
        // CheckpointRestart
        Checkpoint.CheckpointRestore(StartTraj, U, sRNG, pRNG);
      }

      Smearing.set_Field(U);

      HybridMonteCarlo<IntegratorType> HMC(HMCpar, MDynamics, sRNG, pRNG, U);

      for (int obs = 0; obs < ObservablesList.size(); obs++)
        HMC.AddObservable(ObservablesList[obs]);

      // Run it
      HMC.evolve();
    }
  };

  // These are for gauge fields
  typedef BinaryHmcRunnerTemplate<PeriodicGimplR> BinaryHmcRunner;
  typedef BinaryHmcRunnerTemplate<PeriodicGimplF> BinaryHmcRunnerF;
  typedef BinaryHmcRunnerTemplate<PeriodicGimplD> BinaryHmcRunnerD;

  template <class RepresentationsPolicy>
  using BinaryHmcRunnerTemplateHirep = BinaryHmcRunnerTemplate<PeriodicGimplR, RepresentationsPolicy>;

  typedef BinaryHmcRunnerTemplate<ScalarImplR, ScalarFields>
      ScalarBinaryHmcRunner;
}
}
#endif
