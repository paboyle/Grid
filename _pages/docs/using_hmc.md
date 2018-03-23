---
title : "Documentation"
author_profile: false
excerpt: "Using HMC"
header:
  overlay_color: "#5DADE2"
permalink: /docs/theories/hmc/
sidebar:
  nav : docs
---

Using HMC in Grid version 0.7.0

These are the instructions to use the Generalised HMC on Grid version 0.7.0.
Disclaimer: GRID is still under active development so any information provided here can be changed in future releases.

Introduction
=======

The general problem is to generate a Markov Chain distributed according to the action $$S(\psi)$$ in order to compute observables expectaction values.

$$ \langle O \rangle = \frac{1}{Z} \int O e^{-S(\psi)} D\psi $$ 

The Hybrid Monte Carlo approach is to introduce ficticious random momenta to construct an Hamiltonian $$H(\psi)$$ and generate 
new configurations by integrating the corresponding Hamilton equations.

$$H(\psi) = \frac{1}{2} P^2 + S(\psi)$$


Command line options
====================

(relevant file `GenericHMCrunner.h`)
List of command line options, specific for HMC

* ```--StartingType  <string>```  

Choices: HotStart, ColdStart, TepidStart, CheckpointStart

* ```--StartingTrajectory <integer>```

Only for CheckpointStart, ignored otherwise. 
Expected format of the filename is ```<config prefix>.<integer>``` for the configuration and ```<rng file prefix>.<integer>```

* ```--Trajectories <integer>```

Number of trajectories in this run, excluding the thermalization steps. Default: 1.

* ```--Thermalizations <integers>```

Default: 10

* ```--ParameterFile <string>```

The filename for the input parameters deserialisation.

All of them, except the starting trajectory, can be overridden by the input file (but this behaviour can be easily changed by the user writing the source file).

Actions
======

Action names are defined in the file
lib/qcd/Actions.h

Gauge action names list:

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

Scalar action names list

```
ScalarActionR;
ScalarActionF;
ScalarActionD;
```

each of these action accept one single parameter at creation time (beta).
Example for creating a Symanzik action with beta=4.0

	SymanzikGaugeActionR(4.0)

The suffixes R,F,D in the action names refer to the Real
(the precision is defined at compile time by the --enable-precision flag in the configure),
Float and Double, that force the precision of the action to be 32, 64 bit respectively.




