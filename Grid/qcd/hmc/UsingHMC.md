# Using HMC in Grid

These are the instructions to use the Generalised HMC on Grid as of commit `749b802`.
Disclaimer: Grid is still under active development so any information here can be changed in future releases.


## Command line options

(relevant file `GenericHMCrunner.h`)
The initial configuration can be changed at the command line using 
`--StartingType STARTING_TYPE`, where `STARTING_TYPE` is one of
`HotStart`, `ColdStart`, `TepidStart`, and `CheckpointStart`.
Default: `--StartingType HotStart`

Example:
```
./My_hmc_exec  --StartingType HotStart
```

The `CheckpointStart` option uses the prefix for the configurations and rng seed files defined in your executable and the initial configuration is specified by
`--StartingTrajectory STARTING_TRAJECTORY`, where `STARTING_TRAJECTORY` is an integer.
Default: `--StartingTrajectory 0`

The number of trajectories for a specific run are specified at command line by
`--Trajectories TRAJECTORIES`, where `TRAJECTORIES` is an integer.
Default: `--Trajectories 1`

The number of thermalization steps (i.e. steps when the Metropolis acceptance check is turned off) is specified by
`--Thermalizations THERMALIZATIONS`, where `THERMALIZATIONS` is an integer.
Default: `--Thermalizations 10`

Any other parameter is defined in the source for the executable.

## HMC controls

The lines 

```
  std::vector<int> SerSeed({1, 2, 3, 4, 5});
  std::vector<int> ParSeed({6, 7, 8, 9, 10});
```

define the seeds for the serial and the parallel RNG.

The line 

```
  TheHMC.MDparameters.set(20, 1.0);// MDsteps, traj length
```

declares the number of molecular dynamics steps and the total trajectory length.


## Actions

Action names are defined in the directory `Grid/qcd/action`.

Gauge actions list (from `Grid/qcd/action/gauge/Gauge.h`):

```
WilsonGaugeActionR;
WilsonGaugeActionF;
WilsonGaugeActionD;
PlaqPlusRectangleActionR;
PlaqPlusRectangleActionF;
PlaqPlusRectangleActionD;
IwasakiGaugeActionR;
IwasakiGaugeActionF;
IwasakiGaugeActionD;
SymanzikGaugeActionR;
SymanzikGaugeActionF;
SymanzikGaugeActionD;
```

```
ConjugateWilsonGaugeActionR;
ConjugateWilsonGaugeActionF;
ConjugateWilsonGaugeActionD;
ConjugatePlaqPlusRectangleActionR;
ConjugatePlaqPlusRectangleActionF;
ConjugatePlaqPlusRectangleActionD;
ConjugateIwasakiGaugeActionR;
ConjugateIwasakiGaugeActionF;
ConjugateIwasakiGaugeActionD;
ConjugateSymanzikGaugeActionR;
ConjugateSymanzikGaugeActionF;
ConjugateSymanzikGaugeActionD;
```

Each of these action accepts one single parameter at creation time (beta).
Example for creating a Symanzik action with beta=4.0

```
  SymanzikGaugeActionR(4.0)
```

Scalar actions list (from `Grid/qcd/action/scalar/Scalar.h`):

```
ScalarActionR;
ScalarActionF;
ScalarActionD;
```

The suffixes `R`, `F`, `D` in the action names refer to the `Real`
(the precision is defined at compile time by the `--enable-precision` flag in the configure),
`Float` and `Double`, that force the precision of the action to be 32, 64 bit respectively.
